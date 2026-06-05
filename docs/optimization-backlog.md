# Optimization Backlog & Research Directions

Forward-looking companion to [`performance-lessons.md`](performance-lessons.md) (which records what
has been *done*). This file records overlooked optimization opportunities **not yet tried**, the
measurement work that must precede them, and a ready-to-use external deep-research prompt. It is the
output of a multi-agent review (7 independent reviewers: algorithmic, JVM/SIMD, two out-of-the-box,
two internet surveys of libraries + literature, one benchmarking-methodology). Ideas flagged
"converged" were independently surfaced by multiple reviewers — a stronger signal.

Nothing here is committed work; each item needs the measurement gate in §2 before it can be trusted.

## Framing: where we stand

We appear to be **at or beyond the public state of the art for bit-exact numerical Shrake-Rupley.**
The fastest open competitor surveyed (RustSASA) independently converged on our exact architecture
(SoA + cell list + symmetric stencil + SIMD) and is *missing* our per-pair occlusion prune
(lesson 7). So the remaining wins fall into three buckets: (a) a few unexploited *classic* tricks;
(b) one structural reframe of the inner loop; (c) jumping algorithm families (analytic). The kernel's
current bottleneck at p2rank's tess 2 is the **neighbor build** (~47% of CPU at 16 threads), not the
occlusion scan.

A recurring geometric insight drives most Tier-A ideas: **a neighbor occludes a spherical cap of
directions** (axis `diff/|diff|`, cosine half-angle `thresh/|diff|`). The buried test is therefore a
point-vs-cap / cap-coverage problem on a *fixed, small* direction set (42 at tess 2), not inherently a
stream of dot-products.

---

## 1. Optimizations to try

Effort/confidence are the reviewers' estimates; treat as hypotheses to measure, not promises.

### Tier A — bit-exact, ship on the default path if they win

| ID | Idea | Bit-exact | Effort | Notes / first experiment |
|----|------|-----------|--------|--------------------------|
| **A1** | **Neighbor-major "bitmask" occlusion scan** *(converged ×4)*. **Instrumented (tess 2): naive A1 would LOSE** — the direction-major early-exit already resolves a direction in ~7.9 of 42 neighbor-tests (50% via the hint in 1 test), saving 81% of tests; A1 discards that. It is viable **only** with a cheap precomputed cap→direction-mask LUT **plus** a 'stop when all 42 directions buried' atom-level early-exit — which is the only way to capture A3's headroom (see A3). High effort, needs the sphere-index LUT (deep-research Q2). | Yes | High | Don't build naive; needs the cap-mask LUT + mask-full early-exit. |
| ~~**A2**~~ | ~~**Tighter cell-grid cell size**~~ — **TRIED, negative result** (`DevSurfaceV20TightGrid`, ~12–14% slower; bit-exact). The premise was wrong: `CellGrid`'s `2*radius` already *equals* the cutoff, so the grid was optimally sized; a finer grid only adds stencil + counting-sort overhead. See `performance-lessons.md` rung 20. | Yes | — | Do not re-try. |
| **A3** | **Whole-atom trivial reject** *(converged ×3)*. **Instrumented (tess 2): the biggest remaining prize — 61% of atoms are fully buried and consume 57% of the scan's neighbor-tests** (~19% of tess-2 CPU). But the direction-major scan can't atom-early-exit; the unlock is a *sound cheap* coverage certificate, OR the bitmask scan's 'stop when all directions buried' (see A1) — both hard. The cheap sound single-neighbor 'engulfment' test (`d < Rj−Ri`) is low-yield. | Yes, if the skip is sound | High | Same machinery as A1 (cap-mask LUT + mask-full exit), or the cap-coverage algorithm (deep-research Q1). |
| ~~**A4**~~ | ~~**Coverage-based neighbor prune.**~~ **Instrumented → DE-PRIORITIZED: only 6.1% of neighbors bury zero directions at tess 2**, so even a free, perfect prune removes ~6% of scan input — and the detection isn't free. Not worth it (would be another A2/A7-style near-wash). | Yes | — | Headroom too small. |
| **A5** | **DCLM second lattice over the dots** (Eisenhaber 1995) *(converged ×2 — the one classic published trick we lack)*. Spatially bin an atom's surface dots so each neighbor tests only the dots in its cap. Composes with our dedup; **pays at tess 3–4**, less at tess 2. | Yes | Medium | Study the GROMACS `sasa.cpp` / `nsc.c` implementation (LGPL). Related to A1/A4. |
| **A6** | **Vectorize the neighbor-build distance pass.** — **IMPLEMENTED & CONFIRMED** as `DevSurfaceV21SimdBuild` (bit-exact vs V2): **+4.8% at tess 2 single-thread, +3.8% at 16 threads**, neutral at tess 4 (clean JMH on idle box). See `performance-lessons.md` rung 21. V16–V18 fixed the build's memory access but `d² < sumR²` is still scalar. Repack cell-sorted candidates to SoA, masked SIMD compare + compress-store of surviving edges. Attacks the measured #1 bottleneck. | Yes (single compare, no reduction) | Med-High | Compress-store of low-density survivors is the tricky part. |
| ~~**A7**~~ | ~~**Padded-tail elimination** in the scan~~ — **TRIED on the double path, negative result** (`DevSurfaceV22PaddedTail` = V21+A7, ~11% slower single-thread than V21; the padding costs more than the cheap 0–3-element tail). Bit-exact. May still help the 8-wide float scan (untested). See `performance-lessons.md` rung 22. | Yes | — | Don't re-try on the double path. |

### Tier B — opt-in variants / bigger bets (need a tolerance scorecard + p2rank validation)

| ID | Idea | Status | Notes |
|----|------|--------|-------|
| **B1** | **Exact analytic per-atom SASA** *(converged ×4)* — Richmond 1984 Gauss-Bonnet on the neighbor-circle arrangement, or Edelsbrunner power-diagram / alpha-complex (cf. 2024 dSASA GPU paper). `O(k²)`/`O(N log N)`, exact, **wins at tess 3–4 and for area-only consumers**. | Contract-changing (true area ≠ sampled area) | Cheap gate: compute analytic area, compare to discretized tess-2/3/4 on the corpus; pursue only if the gap is within p2rank's feature noise. Needs robust circle-circle / weighted-Delaunay geometry (the hard dependency in Java). |
| **B2** | **Neighbor-list reuse** (Verlet-skin) across solvent-radius sweeps or perturbed/docked poses of one structure. Could make the dominant build ~free. | Bit-exact when applicable | **Gated by a grep of p2rank's call pattern** — if strictly one-shot per structure, dead. Do this 10-min check first. |
| **B3** | **FMA + tree reduction inside the existing float variant.** Lesson 11 rejects FMA for the bit-exact default, but the float path has no exactness to protect (FMA is *more* accurate there). | Within the existing float tolerance | ~1.05–1.15× on GraalVM; re-validate the ≤1.4e-5 bound. Confirm `FloatVector.fma` intrinsifies on Graal. |
| **B4** | **Golden-spiral / Fibonacci sphere points** for a fast approximate mode (all fast tools use it; no dedup needed, free choice of N). | Non-bit-exact | Opt-in only; helps only if p2rank accepts an approximate surface. |

### Investigated and deprioritized (recorded so we don't relearn)
- **GPU offload** — poor fit: 300–5000 atoms is tiny; PCIe/launch latency dominates p2rank's
  per-protein-per-core pattern. Only viable as a separate whole-dataset batch mode.
- **ML/statistical surrogate** for buried-fraction — highest variance; a wrong "confident" prediction
  silently corrupts area. Cheap offline gate: fit a 5-feature model to corpus areas; pursue only if
  residual tails are tiny.
- **Morton/Z-order atom relabeling** — bit-exact but modest, and the output inverse-permutation/point
  re-sort may eat the locality gain.
- **Project Valhalla value types** — the SoA flattening (V12/V19) already captured this; no per-element
  objects remain on the hot path.
- **Lee–Richards slice integration, LCPO/POPS pairwise approximation, ML accessibility predictors** —
  wrong surface/contract or too coarse; produce area only, no point cloud.

---

## 1b. UNFINISHED TRACK — the fully-buried-atom prize (LUT-bitmask scan)

**Status: open, scoped, the single biggest remaining bit-exact opportunity. Resume here.**

The feasibility instrumentation (`ScanInstrumentationTest`, tess 2) established:
- **61% of atoms are fully buried** (emit zero surface points) yet consume **57% of the scan's
  neighbor-tests** (~19% of tess-2 CPU at 16 threads).
- The direction-major scan **cannot** exploit this — it must test all 42 directions to learn an atom is
  fully buried. Its per-direction early-exit + last-occluder hint is already very effective (a direction
  resolves in ~7.9 of 42 neighbor-tests, 50% via the hint in 1 test), so **naive A1 (bitmask) and A4
  (neighbor prune) both lose** — measured, do not re-try them in those forms.

**The idea that could capture it:** a **neighbor-major bitmask occlusion scan** that, for each neighbor,
ORs a precomputed *cap→direction bitmask* into a per-atom accumulator, and **stops as soon as all
directions are buried** (mask full) — an *atom-level* early-exit that fires for exactly the 61% fully-buried
atoms. At tess 2 the accumulator is one `long` (42 bits). Bit-exact if each set bit is confirmed by the
exact `diff·p > thresh` test (the cap mask only *narrows* candidates).

**The hard dependency (why it's not a quick scaffold):** building each neighbor's cap mask must cost
**far less than ~8 dot-tests**, or it loses to the existing early-exit. That needs a cheap, *sound* map
from a neighbor's cap `(axis = diff/|diff|, cos halfangle = thresh/|diff|)` to "which of the fixed 42
directions it buries" — a sphere index / LUT over the fixed tessellation directions (HEALPix / octant grid;
deep-research §3 Q2). A conservative over-inclusive mask + exact confirmation keeps it bit-exact.

**Resume checklist:**
1. Prototype the cap→direction mask: for the fixed tess-2 42-direction set, a structure that, given a cap,
   returns a 64-bit superset mask of buried directions cheaply (target ≪ 8 ops/neighbor).
2. `BitmaskOcclusionScan` (new `OcclusionScan`): OR neighbor masks, confirm set bits exactly, **break when
   all 42 bits set**, `emitWeighted` the survivors. Validate bit-exact vs V3 over the full corpus×config.
3. Benchmark vs V3 at tess 2 (single + 16 threads) on an idle box. The win, if any, comes from the
   fully-buried atoms; expect little at tess 4 (fewer fully-buried, bigger direction set per `long`).
4. If it wins, promote to a production surface; else document as a negative result (like A2/A7).

Related: A5 (DCLM dot-lattice) is the tess-3/4 analogue; B1 (analytic SASA) sidesteps sampling entirely.
The bit-exact *incremental* well is otherwise tapped — A6 (now shipped as `DistinctPackedNumericalSurfaceV3`)
was the catch.

---

## 1c. Additional ideas (post-A6, surfaced this session)

Beyond the original A/B lists. The first three are being implemented + measured now as Dev rungs:

- **C1 — Float-precision neighbor build** (`DevSurfaceV24FloatBuild`, tolerance). **WIN: ~1.6% single / 2.8% 16t at tess 2** (bandwidth signature confirmed). Candidate opt-in fast variant. A6 found the 16-thread build is
  *bandwidth-bound*; doing the build's `d²<sumR²` test in single precision (float candidate coords, 8-wide)
  halves the coordinate traffic — the one lever expected to scale *better* at 16 threads. Non-bit-exact (a
  pair at the cutoff boundary can flip), so tolerance-tested vs V3; the misclassified pairs have ~empty caps,
  so the area error should be tiny.
- **C2 — A7 on the float scan** (`DevSurfaceV25FloatPaddedScan`). **NEGATIVE: ~neutral at tess 2/3, JIT-unstable ~5× slower at tess 4** — the padded loop defeats Vector intrinsics. A7 now dead on double AND float. Padded-tail elimination lost on the 4-wide double scan
  (rung 22) but the scalar tail is proportionally larger on the 8-wide `FloatNumericalSurface` scan.
  Sentinels (`thresh=+∞`) never bury, so it's **bit-for-bit identical to `FloatNumericalSurface`**.
- **C3 — Region-binned last-occluder hint** (`DevSurfaceV23RegionHint`, bit-exact). **NEGATIVE: ~+2–3% slower** — single hint was already enough. The hint already resolves ~50% of
  directions in 1 test; caching a dominant occluder *per sphere octant* (not just the single last one) could
  raise that. Bit-exact (only reorders which neighbor is tested first).

Deferred / deployment-shaped:
- **C4 — Whole-dataset GPU/batch throughput mode.** Single-protein GPU is dead (launch latency), but a
  separate batch mode over thousands of proteins could win at dataset scale. A distinct offering.
- **C5 — Intra-protein parallelism for single-large-protein *latency*** (adaptive gating so it doesn't fight
  the per-protein-per-core throughput path lesson 14 chose).
- **C6 — Solvent-radius-sweep reuse** (subset of B2): reuse the grid + pruned pairs across a radius sweep on
  fixed coordinates. Bit-exact; gated by whether p2rank actually sweeps.
- **C9 — Diagnose the float-scan tess-4 multi-thread collapse.** Re-measurement found the float
  occlusion scan goes ~32× slower at tess 4 / 16 threads (38→1230 ms; double scan flat). Suspected
  denormal-float microcode assists or AVX downclock. Confirm via `perf stat` (fp-assist counters) or a
  flush-to-zero (FTZ/DAZ) experiment; if denormals, set FTZ for the scan thread or clamp tiny floats.
  Until fixed, the float surfaces are tess-2-only. (perf-lessons lesson 5 correction.)
- **C7 — Boundary-adaptive tessellation** / **C8 — Monte-Carlo importance-sampled directions** — approximate
  fast modes; non-bit-exact, need the tolerance scorecard.

---

## 2. Measurement gate — do this BEFORE the next optimization round

This is a prerequisite, not a suggestion. The last three rungs (V16 = 1.03×, V17 = 1.057×,
V18 = 1.031×) are **inside the noise floor of the harness that measured them**, and several Tier-A
ideas are similarly sub-1.1× refinements — they cannot be evaluated with the current setup. The
*correctness* gate (golden + equivalence harness) is production-grade; the *measurement* gate is not.

> **Status: largely built.** The JMH harness (`SurfaceBench` + `bench.sh` + `bench-table.py`), the shared
> `SurfaceCatalog` registry, the fidelity-tiered `scorecard` (accuracy + oracle-free quality metrics:
> duplicate ratio, too-close ratio, point evenness), and the machine-pinning/env-stamp script now exist
> (items 1, 3, 5, 6 below; profiling via JMH's built-in `profilers`, item 2). The legacy median-of-3
> harness is kept for historical reproducibility. A first JMH-measured production table (item 4) is now
> recorded in `performance-lessons.md` (GraalVM 25, governor-pinned). **Remaining nice-to-have:** a fully
> turbo-disabled, longer-iteration run for the official numbers (this box can't disable turbo via its
> cpufreq driver; and re-measure the one `FLOAT/POINTS/tess4` outlier cell).

**Must-do-first:**
1. **JMH with `@Fork ≥ 3` + variance/CI reporting.** The current median-of-3, single-fork,
   `volatile`-blackhole harness can't separate a 1.03× win from turbo jitter. `@Fork` (separate JVMs)
   is the biggest missing thing. Refuse to quote any ratio whose CI overlaps 1.0. (The docs already
   record one warmup artifact: the "1.45×" float win that was really ~1.12×.)
2. **Scripted profiling** — async-profiler + `perf stat -e cycles,instructions,cache-misses,LLC-load-misses`.
   Every microarch claim in the ladder ("gather-bound", "47% in the build", "AVX-512 downclock",
   "bandwidth-bound") is asserted from prose, not instrumented. Can't pick the next lever without it.
3. **Pin the machine** — governor `performance`, disable turbo, `taskset` to fixed cores; stamp every
   result with CPU/microcode/JVM/flags/git-SHA. On a turbo/downclock-bound box, unpinned frequency
   swamps a 1.03× effect.
4. **One canonical ladder table, one harness.** The current table mixes the median-wall and
   throughput harnesses (the doc admits it), so the baseline for the next rung is ambiguous.
5. **Reusable accuracy+speed scorecard for non-bit-exact variants** (generalize
   `FloatNumericalSurfaceTest`): area error, per-atom error distribution, point-set symmetric
   difference, AND JMH speed. This is the gate every Tier-B idea must pass.
6. **Benchmark p2rank's real pattern**, not area-only: drain `surfacePointsXYZ()` /
   `getAllSurfacePoints()` through a Blackhole (the point-materialization path V12/V19 exist for is
   never timed today); treat the 16-thread throughput run as first-class.

**Nice-to-have:** a perf-regression tripwire in CI (nightly, pinned self-hosted runner only — hosted
runners are too noisy for tight thresholds); add one >10k-atom structure to the corpus (all-globular,
≤4779 atoms today) to exercise the build's scaling; report cold-first-build vs warm separately.

---

## 3. External deep-research prompt

Paste verbatim into Claude / ChatGPT / Gemini deep research. Self-contained.

```
I'm optimizing a high-performance solvent-accessible-surface-area (SASA) kernel — numerical
Shrake-Rupley: for each atom, tessellate a unit sphere into points (icosahedral, 42/162/642
DISTINCT directions at levels 2/3/4); a direction is "buried" if any neighboring atom's expanded
sphere occludes it (dot(diff, dir) > threshold); surviving directions form the surface and their
count gives the per-atom area. It's in Java (Vector API / SIMD-256, GraalVM), already heavily
optimized: structure-of-arrays, cell-list neighbor search with a per-pair overlap prune
(d >= R_i + R_j), icosahedral direction DEDUPLICATION (~5.7-6x), a last-occluder hint, and a flat
zero-copy point store. It's ~13x the reference CDK implementation at tessellation level 2 (my real
operating point), and the NEIGHBOR-BUILD now dominates (~47% of CPU at 16 threads). I need bit-exact
reproduction of the numerical result by default; approximate/analytic methods are acceptable only as
separate opt-in variants. Consumer: a protein-pocket-prediction tool that runs this over large
datasets, proteins of 300-5000 atoms, one protein per CPU core.

Research these open questions and give me concrete, cited, implementable answers (algorithms,
pseudocode, complexity, accuracy, and existing open-source code in ANY language I can study or port):

1. CHEAP, SOUND "fully-buried atom" test: a guaranteed-correct criterion to decide that the union of
   a set of spherical caps covers an entire sphere (so the atom contributes zero surface), cheaper
   than sampling. Keywords: union of spherical caps covering a sphere, spherical cap coverage decision
   problem, Eisenhaber buried-atom test, Le Grand-Merz exposed-point culling.

2. BITMASK / loop-inverted occlusion: replacing per-direction dot-products with OR-ing precomputed
   "which of N fixed directions fall in this neighbor's cap" bitmasks (popcount for area). Fast
   point-in-cap queries on a FIXED small direction set (HEALPix / octant / cube-face indexing).
   Software-occlusion-culling and hierarchical-Z analogues (masked occlusion culling) for accumulating
   coverage of a fixed sample set.

3. The DOUBLE CUBIC LATTICE METHOD second-lattice-over-dots (Eisenhaber et al. 1995, JCC 16:273):
   the EXACT mechanism/pseudocode for binning a tessellated sphere's surface dots so each neighbor
   only tests the dots within its occlusion cap, and how to map a cap (axis + angular radius) to the
   overlapping dot-buckets. Read the GROMACS implementation (sasa.cpp / nsc.c) since it's LGPL.

4. EXACT ANALYTIC per-atom SASA from the arrangement of neighbor circles (Gauss-Bonnet): Richmond 1984;
   Fraczkiewicz-Braun 1998 (robust + gradients); Cazals-Loriot 2009 arrangement of circles on a sphere.
   Give the explicit per-atom formula, robust circle-circle intersection ordering, and degeneracy
   handling. What is the practical crossover (neighbor count, tessellation level) where analytic beats
   numerical sampling? How close is true analytic area to icosahedral level-2/3/4 sampled area?

5. POWER DIAGRAM / weighted-alpha-complex exact SASA in O(N log N) (Edelsbrunner; the 2024 dSASA GPU
   paper, arXiv:2401.10462) — and crucially, is there a robust 3D WEIGHTED/REGULAR DELAUNAY
   triangulation usable from Java (pure-Java or via JNI), since that's the hard dependency?

6. NEIGHBOR-LIST REUSE: molecular-dynamics Verlet-list-with-skin and game-physics sort-and-sweep
   broad-phase — applicability to recomputing SASA across perturbed/docked poses of one structure or
   across solvent-radius sweeps (coordinates fixed). Cell-list vs sweep-and-prune vs kd-tree for ~5000
   irregularly-distributed spheres.

7. Survey the FASTEST open-source SASA / molecular-surface implementations in any language (RustSASA,
   FreeSASA, mdtraj, biotite, GROMACS gmx sasa, MSMS, NanoShaper, EDTSurf, dr_sasa, PowerSASA, POPS,
   LCPO) — for each: algorithm, what makes it fast, SIMD/GPU use, and any acceleration (whole-atom
   culling, neighbor reuse, Hilbert/Morton reordering, analytic exposure) that an already-heavily-
   optimized dedup'd SIMD Shrake-Rupley kernel like mine would NOT already have. License-flag anything
   whose code (not just idea) I could reuse.

For each finding: state complexity, expected speedup regime (especially at my tessellation level 2),
whether it preserves bit-exactness or needs a tolerance variant, and a concrete first experiment to
de-risk it. Cite primary sources (papers + code URLs); flag anything you cannot verify rather than
guessing; do not invent benchmark numbers.
```

---

## Survey pointers (verified by the internet reviewers)

- **DCLM second-lattice** — Eisenhaber et al. 1995, *JCC* 16:273 (DOI 10.1002/jcc.540160303); LGPL
  implementation in GROMACS `src/gromacs/trajectoryanalysis/modules/sasa.cpp`.
- **Exact analytic** — Richmond 1984 (*JMB* 178:63); Fraczkiewicz-Braun 1998 (robust + gradients);
  Cazals-Loriot 2009 (arrangement of circles on a sphere); Edelsbrunner alpha-complex; dSASA GPU
  arXiv:2401.10462.
- **Fast implementations** — RustSASA (MIT, Shrake-Rupley + `pulp` SIMD + rayon; *lacks* our per-pair
  prune); FreeSASA (C); mdtraj, biotite; GROMACS DCLM. Approximate/area-only (not applicable to our
  point-cloud contract): LCPO, POPS, pwSASA.

> Citation caveat: some primary PDFs (Eisenhaber full text, Cazals-Loriot) were access-blocked during
> the survey; the DCLM dot-binning mechanism and the circle-arrangement complexity constants should be
> confirmed against the primary papers before implementation.
