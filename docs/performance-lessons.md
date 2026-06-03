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

Each step is a separate class; the prior ones are left untouched as a measurable baseline. The
shipped production impl is `FasterNumericalSurface` (used by p2rank); every step we developed is an
experimental class named `DevSurfaceV<n><Increment>` (prefix + version + the increment it adds), so
the ladder is ordered and self-documenting without 9-qualifier class names. Mean speedup vs CDK
`NumericalSurface` (GraalVM 25 single-thread, this machine: 16 physical cores, AVX-512):

| # | class | increment | tess 4 | tess 3 | tess 2 (p2rank) |
|---|---|---|---|---|---|
| - | `FasterNumericalSurface` | reference optimized impl (production) | 1.22x | 1.24x | 1.34x |
| 1 | `DevSurfaceV1Soa` | structure-of-arrays scratch (also the shared engine) | 1.47x | 1.50x | 1.71x |
| 2 | `DevSurfaceV2Grid` | flat cell-grid neighbor index | 1.46x | 1.55x | 1.87x |
| 3 | `DevSurfaceV3Sorted` | sort neighbors by threshold | 4.62x | 3.50x | 2.31x |
| 4 | `DevSurfaceV4Hinted` | last-occluder hint (drops the sort) | 6.92x | 4.87x | 3.23x |
| 5 | `DevSurfaceV5Symmetric` | compute each neighbor pair once | 7.13x | 5.34x | 3.67x |
| 6 | `DevSurfaceV6LowAlloc` | bufferless two-pass build (GC side-branch off V5) | 6.89x | 5.10x | 3.42x |
| 7 | `DevSurfaceV7Simd` | SIMD occlusion scan (Vector API) + scalar fallback | 10.67x | 7.24x | 4.67x |
| 8 | `DevSurfaceV8Pruned` | per-pair occlusion cutoff + 256-bit lanes | 12.19x | 8.86x | 5.78x |
| 9 | `DevSurfaceV9Dedup` | scan distinct tessellation directions | 19.06x | 14.90x | 10.46x |
| 10 | `DevSurfaceV10CachedMap` | cache the direction mapping process-wide | 20.10x | 14.22x | 10.53x |
| 11 | `DevSurfaceV11CachedTess` | cache the tessellation arrays | 20.31x | 15.38x | 10.59x |
| 12 | `DevSurfaceV12Flat` | flat `double[]` point output (B) + cached VdW radii (D) | 18.30x | 14.91x | 10.64x |
| 13 | `DevSurfaceV13Arena` | per-thread scratch arena (A) on V11 | 21.83x | 15.47x | 10.67x |
| 14 | `DevSurfaceV14ArenaFlat` | per-thread scratch arena (A) on V12 | 19.72x | 15.03x | 10.69x |
| 16 | `DevSurfaceV16DirectNbr` | copy-free CSR neighbor access (#2) + cached VdW, on V11 | see below | | |
| 17 | `DevSurfaceV17PackedNbr` | packed single-`int[]` edge buffer, on V16 | see below | | |
| 18 | `DevSurfaceV18SortedCoords` | cell-sorted candidate coordinates (sequential read, no gather), on V17 **(compute champion)** | see below | | |
| 19 | `DevSurfaceV19FlatStore` | V18 + flat `double[]` point store / zero-copy `surfacePointsXYZ()` (storage only, compute-neutral) | see below | | |

Net at p2rank's operating point (tess 2): **~10.6x CDK** (V11), about **2.3x** over the first
vectorized step (V7), all bit-for-bit identical. V12-V14 (flat output, arena) are allocation/GC
reductions: bit-exact and lower-allocation, but ~flat on single-thread speed here (the kernel is
bandwidth/turbo-bound, not GC-bound). On HotSpot 25 (C2) the same ladder tops out around 11-12x at
tess 4 / ~9x at tess 2; Graal's JIT optimizes these SoA/SIMD/dedup loops markedly better, and the gap
widens with tessellation.

The table's column is the per-structure median-wall-time harness (`./gradlew benchmark`); V16 was added
later and measured with the steady-state aggregate-throughput harness instead, which reports higher
absolute multiples (e.g. V11 at ~11.5x/26x CDK for tess 2/4 single-thread there), so V16's figures are
quoted relative to V11 measured on the *same* harness rather than dropped into the column above (the two
harnesses are not interchangeable; regenerate the whole column with one harness if you need V16 in it).

**Rung 16, `DevSurfaceV16DirectNbr`.** Two changes on V11, both bit-for-bit identical and both "less
work, no added compute": (#2) the pruned neighbor source already stores neighbors in one CSR array, so
read that array directly per atom instead of copying each slice into a reused `IntArrayList` per query
(the copy, `IntArrayList.add`, was ~15% of single-thread CPU, ~2/3 of it this per-query copy); and (#3)
resolve VdW radii through the process-wide cache on the engine side too. Measured vs V11 (GraalVM 25,
aggregate-throughput harness, median of 7): **1.03x** at tess 2 single-thread, **1.015x** at tess 2 / 16
threads, **1.01x** at tess 4 single-thread, **1.05x** at tess 4 / 16 threads - a small but consistent win
everywhere, no regression. Enabling the copy-free path is opt-in (`directNeighbors` engine flag), so V11
and the others keep the copy path and stay byte-identical baselines; the V11 row re-measured at 8.11 vs
8.02 Matoms/s confirms the shared-engine change is neutral.

**Rung 17, `DevSurfaceV17PackedNbr`.** Profiling V16 showed that once the per-query copy was
gone, the construction's *own* edge buffering became the top non-scan hot method: the pruned build
records each kept edge with `edgeI.add(i); edgeJ.add(j)` into two HPPC `IntArrayList`s, which the V16
profile put at ~25% of CPU single-thread and **~34% at 16 threads** (and most of the 63% `int[]`
allocation). V17 packs each kept edge as an interleaved `(i, j)` pair into one cursor-managed `int[]`
grown by doubling: one sequential write-stream instead of two, one capacity check per edge instead of
two, no per-add wrapper overhead. It keeps the **single** distance pass, so unlike V15 it adds no
arithmetic. Measured vs V16 (same harness, median of 7): **1.057x** at tess 2 single-thread, **1.074x**
at tess 2 / 16 threads, **1.005x** at tess 4 single-thread, **1.017x** at tess 4 / 16 threads (vs V11:
1.076x / 1.090x / 1.018x / 1.044x). The win is largest at tess 2 / 16 threads - exactly where the
neighbor build is the biggest fraction and the bandwidth-bound regime rewards a single write-stream -
and shrinks at tess 4 where the SAS scan dominates. No regression anywhere. This is the V16 profile's
predicted #1 lever (make the edge buffer cheap, don't remove it), and it landed.

**Champion: `DevSurfaceV18SortedCoords` (rung 18).** Profiling V17 showed the build was still the top
hot method at 16 threads (~47%), with ~87% of its samples on the inner cell-scan loop reading each
candidate's coordinates `ax[j], ay[j], az[j], expandedR[j]` at the candidate's *atom* index
`j = cellAtoms[p]` - a 4-way **random gather**, since `j` runs in cell order, not atom order. The build
had become memory-latency-bound on that gather. V18 builds, once, a cell-sorted interleaved coordinate
array (`coords[4p..4p+3] = (x,y,z,expandedR)` of `cellAtoms[p]`) with a single gather pass, then the
distance pass reads it **sequentially** as it scans a cell - one streaming read per candidate instead of
four random loads. Each candidate's coordinates are read many times (once per appearance across neighbor
cells) but the sorted copy is built once, so the trade favors it. Measured vs V17 (same harness, median
of 7): **1.031x** at tess 2 single-thread, **1.056x** at tess 2 / 16 threads (vs V11: 1.131x / 1.169x).
Largest at 16 threads, the memory-bandwidth-bound regime that most rewards sequential access; smaller at
tess 4 where the SAS scan dominates. No regression. `CellGrid` is untouched (the source builds its own
sorted copy), so the other variants stay byte-identical baselines. This is the V17 profile's predicted
lever (turn the candidate gather into a sequential stream), and the last random-access hot path in the
build: at p2rank's tess 2 the ladder now reaches ~13x CDK single-thread / ~13.2x at 16 threads.

**Rung 19, `DevSurfaceV19FlatStore`.** The first rung that changes point *storage* rather than the
compute. The neighbor build, scan, and tessellation are V18's, bit-for-bit; the only change is swapping
the `Point3d`-per-point list store for `FlatSurfacePointStore` (one flat `double[]`, `Point3d`
materialized lazily). The point is the zero-copy `PackedSurfaceAccess` delivery path: a caller that
consumes raw coordinates in bulk (p2rank turns each point into its own three-double `Point`) gets them
via `surfacePointsXYZ()` with **no copy** and never allocates a `Point3d` or the `Point3d[]`. Through the
`Point3d`-valued `MolecularSurface` accessors it behaves exactly like V18 (lazy materialization), so on
this box it is **throughput-neutral** single-thread (bandwidth-bound, not GC-bound); its payoff is the
lower allocation footprint and the zero-copy bulk path, an at-scale / GC-pressure win rather than a raw
single-thread speedup. Output is bit-for-bit identical to `FasterNumericalSurface`. This builds out the
"flat `double[]` surface-point output" idea formerly under Open ideas. V18 remains the *compute*
champion; V19 is the recommended shape when bulk coordinate delivery matters.

**Side-branch (regression, kept as a documented negative result):** `DevSurfaceV15LeanNbr` (rung 15)
applies the V6 bufferless two-pass build to V11's pruned neighbor source (#1) plus engine-side cached
VdW radii (#3). It is bit-for-bit identical but **slower than V11**: ~0.86x at tess 2 single-thread,
~0.88x at tess 2 / 16 threads, ~0.94x at tess 4 single-thread, ~1.00x (break-even) at tess 4 / 16
threads. It is a side-branch off V11, not a rung: see lesson 13 for why a profile that flagged the
neighbor build as the top hot method and top allocator still did not make allocation-reduction pay. V16
is the same neighbor build attacked the right way (remove a redundant copy instead of trading compute
for allocation), and it wins where V15 lost.

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

5. **A `float` verdict pays off on GraalVM, not HotSpot.** Narrowing only the occlusion *comparison*
   to single precision (`FloatNumericalSurface`: 256-bit `FloatVector`, 8 lanes vs 4, and half the
   neighbor-scratch memory traffic, while point positions and areas stay `double` — the coordinate
   subtraction stays in `double` to avoid cancellation, and `emitWeighted` uses the `double`
   directions) is **~1.05–1.14x** faster than the `double` distinct surface on **GraalVM 25** (Graal
   JIT) — uniformly across 1 and 16 threads and tess 2–4. On **HotSpot 25** (Oracle C2) the *same*
   bytecode is **neutral to slightly slower** (0.93–1.01x): C2's `double` Vector path is already good
   and the per-atom `double→float` narrowing is not amortized, worst at low tess where neighbor lists
   are short and the 8-wide loop falls to its scalar tail. So the float win is JIT-dependent, worth
   shipping only because the deployment JVM is GraalVM — and even there it is modest, because the
   (unchanged) neighbor build is a large fixed share of the time. Not bit-exact: a tessellation point
   within `float` epsilon of the occlusion boundary may flip survival, but over the whole corpus
   (incl. the degenerate solvent=0 surface) total-area error stayed ≤ 1.4e-5 and the distinct
   point-set symmetric difference ≤ 1.4e-5 of points — well within SAS discretization error. (An
   earlier 1.45x figure came from a 10-build-per-round micro-bench and was a warmup artifact; the
   200-build median put it at ~1.12x. Corollary to lesson 15.)

### Algorithmic wins

6. **Deduplicate tessellation directions (the biggest single win).** CDK's icosahedral tessellation
   emits each direction about **5.7x** (240 points but only 42 distinct directions at tess 2; 960/162
   at tess 3; 3840/642 at tess 4), because it returns all three vertices of every shared-vertex
   triangle. The buried verdict depends only on direction, so evaluate each distinct direction once
   and re-expand surviving directions into the original point order and multiplicity. Same multiset,
   same order, same areas, ~5.7x fewer scans. This was the single largest speedup (+79% at tess 2).
   The mapping is a pure function of the tessellation, so memoize it per build.

7. **The neighbor cutoff was over-inflated; tighten it per pair.** A neighbor `j` can bury a point of
   atom `i` only if their expanded spheres overlap, `d(i,j) < R_i + R_j`. When `d >= R_i + R_j`,
   `diff . p <= |diff| <= thresh` for every `p`, so `j` never buries anything and can be dropped with
   zero effect on the result. The cell grid's global cutoff (`2*(maxRadius+solvent)`) is much looser
   and admits many such never-firing neighbors. Filtering to the exact per-pair bound roughly halved
   the neighbor count fed to the (dominant) scan. Bit-exact, with a one-line geometric proof.

8. **A cheap hint beat an expensive sort.** Sorting each atom's neighbors by occlusion strength
   (`DevSurfaceV3Sorted`) helped, but profiling showed the sort had become ~42% of runtime. Replacing it
   with a last-occluder hint (remember the neighbor that buried the previous point; test it first)
   captured the same early-exit benefit at near-zero cost, because consecutive tessellation points are
   spatially coherent. The hint variant beat the sort variant at every tessellation level.

9. **Exploit symmetry in the neighbor build.** The neighbor relation is symmetric, so a half-stencil
   that computes each unordered pair's distance once (recording both directions) halves the distance
   math versus querying every atom independently.

### Numerics and what does or does not break bit-exactness

10. **Safe (container/order only), keep as default:** changing the result container (e.g. `List<Point3d>`
   to a flat `double[]`), reordering neighbors, reordering the occlusion scan, dropping neighbors that
   provably never bury, deduplicating directions. None touch the arithmetic that produces a surviving
   point.

11. **Unsafe (changes the floating-point result), opt-in variant only:** FMA, `float` instead of
    `double`, and reassociating sums. Note the subtle one: **`x / c` is not bit-identical to
    `x * (1/c)`** (the reciprocal rounds, then the multiply rounds again). We rejected a
    reciprocal-multiply micro-opt for the threshold division for exactly this reason; it can flip a
    `>` comparison in the last ULP on a boundary point. The `float`-verdict variant (lesson 5) is the
    one place this opt-in is actually shipped, behind tolerance tests rather than the bit-exact harness.

### Parallelism and scaling

12. **16-thread scaling is bound by memory bandwidth and turbo downclock, not GC.** At 16 physical
    cores throughput scaled ~13.2x (about 82% efficiency); the loss is shared memory bandwidth (the
    scan streams contiguous arrays) plus the all-core clock being lower than single-core boost. GC was
    only ~2-4% of wall time. Part of the 18% loss (the downclock half) is a hardware limit and not
    recoverable in software.

13. **Reducing allocation does not help throughput when GC is not the bottleneck.** The
    low-allocation variant cut allocation ~63% and GC ~6-7x, but was ~5% *slower* at 16 threads,
    because its extra compute (a second distance pass) cost more than the GC it saved. Lower
    allocation only wins under genuine GC/heap pressure (small heap, high core count, throughput
    collector). On a RAM-rich box it loses. Always confirm what the bottleneck actually is before
    optimizing for it.

    **Confirmed a second time, against a profile that pointed the other way (`DevSurfaceV15LeanNbr`).**
    Profiling V11 at tess 2 showed the *neighbor build* had become the dominant cost - 25% of CPU
    single-thread, **44% (the hottest method) at 16 threads**, and 51-61% of all allocation (its `int[]`
    edge buffers). The obvious read is "allocation-bound, so go bufferless." We applied the V6 two-pass
    trick to the pruned source: it removed the edge buffers but doubled the per-pair distance pass, and
    it **regressed 12-14% at tess 2** (both thread modes), breaking even only at tess 4 / 16 threads
    where the occlusion scan dominates and GC pressure is highest. The lesson: a hot method that is also
    the top allocator is *not* therefore allocation-bound. The neighbor build is CPU-bound on its
    distance arithmetic; trading that arithmetic for less allocation loses even at 16 threads, because
    (per lesson 12) the box is bandwidth/turbo-bound, not GC-bound. The right lever for that hot method
    is to cut *work* with no added compute (e.g. eliminate the redundant per-query CSR-to-`IntArrayList`
    copy), not to trade compute for allocation.

    **The predicted lever worked (`DevSurfaceV16DirectNbr`).** That same CSR source was already holding
    every atom's neighbors contiguously, but the `NeighborSource` contract forced a per-query copy of
    each slice into a reused `IntArrayList`. Exposing the array directly (a `DirectNeighborSource`
    capability + an opt-in engine flag) removed the copy with zero added arithmetic and won at every
    point measured (1.01-1.05x V11), most at tess 4 / 16 threads. The contrast with V15 is the whole
    lesson: attacking the same hot method by *removing redundant work* wins where attacking it by
    *trading compute for allocation* lost. Diagnose what a hot method is bound by before optimizing it.

    **And again, on the edge buffer itself (`DevSurfaceV17PackedNbr`).** Removing the copy exposed that
    the build's own `edgeI.add(i); edgeJ.add(j)` into two `IntArrayList`s was now the top non-scan hot
    method (~34% at 16 threads). V15 had tried to *remove* that buffer (two-pass, recompute distances)
    and lost. V17 instead makes it *cheap* - one packed `int[]` with a manual cursor, single distance
    pass - and wins (up to 1.074x V16 at tess 2 / 16 threads). Same hot method, third time: the buffer
    was never the problem, the wrapper overhead and the second write-stream were. Make the necessary
    work cheap; don't pay to avoid it with more compute.

14. **Intra-surface parallelism conflicts with across-item parallelism.** p2rank already runs one
    protein per core. Parallelizing a single surface across threads would oversubscribe and worsen the
    bandwidth/downclock contention in that case; it only helps the single-protein / idle-core regime,
    and only with adaptive gating.

### Methodology

15. **Measure, do not theorize.** Three times a confident hypothesis was overturned by a five-minute
    measurement: (a) the vector-escape fear (refuted by JFR allocation counts), (b) the
    low-allocation-wins-at-scale fear (refuted by a throughput run), (c) the AVX-512 regression
    (discovered only because the benchmark caught the vectorized variant dropping below CDK). The
    benchmark and the profiler decided the design, not intuition.

16. **Lean on the equivalence harness.** Because every variant is checked bit-for-bit against CDK over
    the corpus, we could refactor the hottest code aggressively and trust the tests to catch any
    divergence (including subtle ones like point ordering after the dedup). Build the oracle first;
    optimize fearlessly after.

17. **Benchmark hygiene.**
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
- **Reciprocal-multiply for the threshold division.** Breaks bit-exactness (see lesson 11).
- **Parameterizing the vector species.** Breaks intrinsification (see lesson 1).

---

## Open ideas (not yet built)

- ~~**Flat `double[]` surface-point output**~~ **Built** as `DevSurfaceV19FlatStore` (and shipped in the
  production `Packed`/`DistinctPacked` surfaces via `FlatSurfacePointStore` + `PackedSurfaceAccess`):
  one flat `double[]` instead of one `Point3d` per point, with a zero-copy `surfacePointsXYZ()` bulk path
  for consumers like p2rank that discard the `Point3d` objects. Bit-exact (container only). As predicted,
  it is throughput-neutral single-thread here (bandwidth-bound, not GC-bound) — the win is allocation /
  at-scale GC pressure. See rung 19 above.
- **Global tessellation/template cache keyed by tessellation level**, including the dedup mapping, so
  thousands of builds share one immutable tessellation instead of rebuilding it each construction.
  Helps small proteins and batch runs most.
- ~~**`float32` opt-in variant.**~~ **Built** as `FloatNumericalSurface` (float verdict, double
  positions/areas; tolerance-tested). The bandwidth/lane argument held only partially: it is ~1.05–1.14x
  on GraalVM and neutral-to-slower on HotSpot C2, because the per-atom `double→float` narrowing is not
  free and the neighbor build is a large unchanged share. See lesson 5.
