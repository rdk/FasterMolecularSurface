# Performance Optimization Lessons

A field guide distilled from optimizing the numerical solvent-accessible-surface (SAS) kernel in
this library. It records what worked, what surprised us, and what is explicitly not worth doing, so
future work does not relearn it the hard way.

The kernel: for each atom, tessellate a unit sphere into points; a point is *buried* if any
neighboring atom occludes it (`diff . p > thresh`); surviving points form the surface and their count
gives the per-atom area. The work is `O(atoms x tessPoints x neighbors)` with an early exit.

> [!IMPORTANT]
> **Bit-for-bit is the contract.** Every optimization below is verified to reproduce CDK's reference
> `NumericalSurface` exactly (golden baseline + `assertBitForBit` equivalence over the 10-structure
> corpus at multiple tessellation levels and solvent radii). Optimizations that would change the
> floating-point result ship only as separate, opt-in variants with their own tolerance tests, never
> as the default. The test harness is what made aggressive optimization safe: a wrong idea fails
> loudly instead of silently corrupting a surface.

---

## The optimization ladder

Each step is a separate class; the prior ones are left untouched as a measurable baseline. Mean
speedup vs CDK `NumericalSurface` (GraalVM 25, this machine: 16 physical cores, AVX-512):

| variant | technique added | tess 4 | tess 3 | tess 2 (p2rank) |
|---|---|---|---|---|
| `FasterNumericalSurface` | reference optimized impl | 1.22x | 1.25x | 1.32x |
| `SoaNumericalSurface` | structure-of-arrays scratch | 1.47x | 1.50x | 1.71x |
| `GridSoaNumericalSurface` | flat cell-grid neighbor index | 1.46x | 1.55x | 1.87x |
| `OrderedGridSoa…` | sort neighbors by threshold | 4.62x | 3.50x | 2.31x |
| `HintedGridSoa…` | last-occluder hint (no sort) | 6.92x | 4.87x | 3.23x |
| `SymmetricHintedGridSoa…` | compute each neighbor pair once | 7.13x | 5.34x | 3.67x |
| `LowAllocSymmetricHintedGridSoa…` | bufferless two-pass build (GC tradeoff) | 6.89x | 5.10x | 3.42x |
| `VectorizedSymmetricHintedGridSoa…` | SIMD occlusion scan (Vector API) | 10.68x | 7.24x | 4.61x |
| `PrunedVectorizedSymmetricHintedGridSoa…` | per-pair occlusion cutoff + 256-bit SIMD | 13.16x | 8.85x | 5.72x |
| `DedupVectorizedSymmetricHintedGridSoa…` | scan distinct tessellation directions | **18.61x** | **14.76x** | **10.25x** |

Net at p2rank's operating point (tess 2): **10.25x CDK**, about **2.2x** over the first vectorized
variant, all bit-for-bit identical.

---

## Hard-won lessons

### SIMD and the Java Vector API

1. **The `VectorSpecies` must be a `static final` constant.** Passing the species as an instance
   field or method parameter defeats the JIT's vector intrinsics and collapses throughput by roughly
   **10x** (measured: the variant fell *below* CDK). To support two lane widths we duplicated the
   scan into two classes, each baking in its own `static final` species, rather than parameterizing
   one. The ~45 lines of duplication are the price of intrinsification.

2. **Widest is not fastest: AVX-512 frequency downclock is real.** `SPECIES_PREFERRED` is 512-bit on
   this box, but once neighbor lists are short (after the cutoff below) the 512-bit bursts are too
   brief to amortize AVX-512's license downclock and the scalar<->512-bit transition penalties, and
   the kernel ran *slower than scalar*. Pinning to **256-bit** fixed it. Rule of thumb: for short,
   bursty SIMD loops on AVX-512 hardware, measure 256-bit against `SPECIES_PREFERRED`; do not assume
   wider wins. The unpruned variant (longer lists) is fine at 512-bit, so the right width is
   workload-dependent.

3. **Vectors do not escape if they stay local to the scanning method.** We feared the strategy-object
   indirection (the scan is reached through an interface field, megamorphic across variants) would
   prevent inlining and force `DoubleVector` boxes onto the heap. JFR allocation profiling refuted it:
   ~6 vector-box samples across 800 builds of a 4779-atom protein. The `DoubleVector`s never escape
   `collect()`, so they are scalar-replaced into registers regardless of whether `collect()` itself
   inlines into the caller. Megamorphic dispatch to a vectorized method is fine.

4. **Bit-exact SIMD means no FMA and a fixed reduction order.** The dot product uses lane-wise
   `mul` then `add` in the same order as the scalar code (`(dx*px + dy*py) + dz*pz`); a fused
   multiply-add rounds once instead of twice and diverges from the reference. The buried predicate is
   reduced with `anyTrue` (an OR), which is order-independent, so neighbor scan order is free to change.

### Algorithmic wins

5. **Deduplicate tessellation directions (the biggest single win).** CDK's icosahedral tessellation
   emits each direction about **5.7x** (240 points but only 42 distinct directions at tess 2; 960/162
   at tess 3; 3840/642 at tess 4), because it returns all three vertices of every shared-vertex
   triangle. The buried verdict depends only on direction, so evaluate each distinct direction once
   and re-expand surviving directions into the original point order and multiplicity. Same multiset,
   same order, same areas, ~5.7x fewer scans. This was the single largest speedup (+79% at tess 2).
   The mapping is a pure function of the tessellation, so memoize it per build.

6. **The neighbor cutoff was over-inflated; tighten it per pair.** A neighbor `j` can bury a point of
   atom `i` only if their expanded spheres overlap, `d(i,j) < R_i + R_j`. When `d >= R_i + R_j`,
   `diff . p <= |diff| <= thresh` for every `p`, so `j` never buries anything and can be dropped with
   zero effect on the result. The cell grid's global cutoff (`2*(maxRadius+solvent)`) is much looser
   and admits many such never-firing neighbors. Filtering to the exact per-pair bound roughly halved
   the neighbor count fed to the (dominant) scan. Bit-exact, with a one-line geometric proof.

7. **A cheap hint beat an expensive sort.** Sorting each atom's neighbors by occlusion strength
   (`OrderedGridSoa`) helped, but profiling showed the sort had become ~42% of runtime. Replacing it
   with a last-occluder hint (remember the neighbor that buried the previous point; test it first)
   captured the same early-exit benefit at near-zero cost, because consecutive tessellation points are
   spatially coherent. The hint variant beat the sort variant at every tessellation level.

8. **Exploit symmetry in the neighbor build.** The neighbor relation is symmetric, so a half-stencil
   that computes each unordered pair's distance once (recording both directions) halves the distance
   math versus querying every atom independently.

### Numerics and what does or does not break bit-exactness

9. **Safe (container/order only), keep as default:** changing the result container (e.g. `List<Point3d>`
   to a flat `double[]`), reordering neighbors, reordering the occlusion scan, dropping neighbors that
   provably never bury, deduplicating directions. None touch the arithmetic that produces a surviving
   point.

10. **Unsafe (changes the floating-point result), opt-in variant only:** FMA, `float` instead of
    `double`, and reassociating sums. Note the subtle one: **`x / c` is not bit-identical to
    `x * (1/c)`** (the reciprocal rounds, then the multiply rounds again). We rejected a
    reciprocal-multiply micro-opt for the threshold division for exactly this reason; it can flip a
    `>` comparison in the last ULP on a boundary point.

### Parallelism and scaling

11. **16-thread scaling is bound by memory bandwidth and turbo downclock, not GC.** At 16 physical
    cores throughput scaled ~13.2x (about 82% efficiency); the loss is shared memory bandwidth (the
    scan streams contiguous arrays) plus the all-core clock being lower than single-core boost. GC was
    only ~2-4% of wall time. Part of the 18% loss (the downclock half) is a hardware limit and not
    recoverable in software.

12. **Reducing allocation does not help throughput when GC is not the bottleneck.** The
    low-allocation variant cut allocation ~63% and GC ~6-7x, but was ~5% *slower* at 16 threads,
    because its extra compute (a second distance pass) cost more than the GC it saved. Lower
    allocation only wins under genuine GC/heap pressure (small heap, high core count, throughput
    collector). On a RAM-rich box it loses. Always confirm what the bottleneck actually is before
    optimizing for it.

13. **Intra-surface parallelism conflicts with across-item parallelism.** p2rank already runs one
    protein per core. Parallelizing a single surface across threads would oversubscribe and worsen the
    bandwidth/downclock contention in that case; it only helps the single-protein / idle-core regime,
    and only with adaptive gating.

### Methodology

14. **Measure, do not theorize.** Three times a confident hypothesis was overturned by a five-minute
    measurement: (a) the vector-escape fear (refuted by JFR allocation counts), (b) the
    low-allocation-wins-at-scale fear (refuted by a throughput run), (c) the AVX-512 regression
    (discovered only because the benchmark caught the vectorized variant dropping below CDK). The
    benchmark and the profiler decided the design, not intuition.

15. **Lean on the equivalence harness.** Because every variant is checked bit-for-bit against CDK over
    the corpus, we could refactor the hottest code aggressively and trust the tests to catch any
    divergence (including subtle ones like point ordering after the dedup). Build the oracle first;
    optimize fearlessly after.

16. **Benchmark hygiene.**
    - Compare **ratios vs CDK on the same JVM**, not absolute milliseconds across runs (machine load,
      turbo, and JIT warmup make cross-run absolutes unreliable).
    - The harness is coarse (median of a few runs after warmup), good for relative comparison, not for
      reporting absolute numbers as if from JMH.
    - To benchmark on the deployment JVM (GraalVM 25) while the build pins a Java 17 toolchain,
      temporarily disable the toolchain block so the launcher JVM both compiles and runs, then revert.
    - For allocation questions, `com.sun.management.ThreadMXBean.getThreadAllocatedBytes` gives a clean
      per-thread byte delta without JFR overhead; for hot-spot questions, JFR execution samples
      aggregated by leaf method; for escape questions, JFR `ObjectAllocationSample` filtered by class.

---

## Explicitly not worth doing (closed, with reason)

- **Cache-blocking the tessellation x neighbor loops.** The per-atom working set (neighbor scratch +
  tessellation arrays) already fits L1; there is nothing to block.
- **Optimizing `firstTrue()` in the SIMD scan.** It runs only on a burying block, which the scalar
  hint already short-circuits; it is cold.
- **NUMA / first-touch tuning.** This box is single-socket; the read-only tessellation arrays are tiny
  and live in shared LLC. Revisit only on multi-socket hardware.
- **Caching the neighbor grid across builds.** Coordinates differ per protein, so the grid genuinely
  cannot be reused; it is correctly per-build.
- **Reciprocal-multiply for the threshold division.** Breaks bit-exactness (see lesson 10).
- **Parameterizing the vector species.** Breaks intrinsification (see lesson 1).

---

## Open ideas (not yet built)

- **Flat `double[]` surface-point output** instead of one `Point3d` per point (and the per-atom
  `ArrayList` copy). Bit-exact (container only); cuts the dominant allocation and a downstream copy in
  p2rank, which consumes coordinates in bulk and discards the `Point3d` objects. Mostly an
  at-scale / GC-pressure win.
- **Global tessellation/template cache keyed by tessellation level**, including the dedup mapping, so
  thousands of builds share one immutable tessellation instead of rebuilding it each construction.
  Helps small proteins and batch runs most.
- **`float32` opt-in variant.** Halves memory traffic and doubles SIMD lanes, the most direct lever on
  the bandwidth ceiling; single precision is well within the tessellation discretization error for
  SAS. Breaks bit-exactness, so it ships separately with tolerance tests (keep the coordinate
  subtraction in `double`, narrow to `float` for the scan, to avoid cancellation).
