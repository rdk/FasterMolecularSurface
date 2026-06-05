# Autoresearch run — final report

**Scope:** maximize throughput of the molecular SASA kernel at **tessellation 2 & 3**, single-thread and
16-thread (p2rank's one-protein-per-core pattern), producing the **surface point set** (not just area),
bit-exact by default or within a documented tolerance.
**Date:** 2026-06-05. **Branch:** `autoresearch`. **Result:** 8 phases, **0 promotable speedups** — the
tess-2/3 single-protein surface is confirmed at a **strong local optimum**. No production code changed
(only test-scope instrumentation, docs, and the journal). Concluded by user decision ("accept the optimum").

Full detail: `autoresearch/LOG.md` (per-phase journal). Raw benchmarks: `autoresearch/results/`.

---

## Phases

| # | Lead | Method | Verdict |
|---|---|---|---|
| 1 | LUT-bitmask "fully-buried atom" prize (§1b — the long-standing "biggest open bit-exact prize") | Count-based feasibility + a sound cap→direction LUT (`BitmaskFeasibilityTest`/`CapDirectionLut`) | **Closed-negative.** Perfect-mask floor is only 11.9%/18% of current scalar tests, but the baseline is already 4-wide SIMD; a sound cheap LUT is 85–92% false positives → gather-bound, loses to SIMD. |
| 2 | Float-scan scaling at tess 3 / many threads | JMH `-prof gc`, `FLOAT_V2` vs `V3`, -t 1/4/8/16; an escape-analysis probe | **C9 CONFIRMED + extended to tess 3.** Float collapses 7×/12.5×/23× at 4/8/16t; `alloc.rate.norm` jumps 13.6 MB/op → 1.084 GB/op at ≥2 threads. Mechanism: Vector-API float (8-lane) intrinsics deopt to boxing under concurrency (EA-off probe ruled out generic escape analysis). Not fixable in source → **float is tess-2-only.** |
| 3 | A5 — DCLM second dot-lattice | Survivor-angle instrumentation (`DclmFeasibilityTest`) | **Closed.** Survivors are 25% of tess-3 scan tests and a spatial bin could cut them ~40×, but that saving is already inside Phase 1's neighbor-major total, which loses to the SIMD baseline. |
| 4 | B1 — analytic per-atom SASA | Area-vs-converged-tessellation gate (`AnalyticGateTest`) | **Closed.** Analytic ≠ sampled (0.3% total / 21% per-atom gap at tess 2); area-only (no point set); O(neighbors²) trig, no SIMD → slower. |
| 5 | Re-profile the build | Phase-decomposed instrumented clone (`BuildProfileTest`) | **Build tapped.** Distance pass is 90% of build, already SIMD (A6) and locality-optimized (V18 cell-sorted coords); 6.27 candidate tests/survivor (16% survive) is inherent uniform-grid geometry; the tighter-grid fix is the A2 negative. |
| 6 | C1 — float build (`V24` = float build + double scan) at tess 3 | Tolerance gate (`FloatBuildToleranceGateTest`) + JMH | **Negative.** Bit-identical to V3 on the corpus, but only ~1.7% faster at 1t and neutral/noisy at 16t — build-only wins dilute at tess 3 where the double scan dominates. Stays a tess-2 variant. |
| 7 | Source novel approaches (generative) | `deep-research` workflow: 5 angles, 19 sources, 24/25 claims adversarially verified | **No new lever.** Power-diagram / weighted-alpha-complex exact SASA is area-only (re-derives B1 via a different method); robust regular/power Delaunay is native-C++-only (CGAL GPLv3 / Geogram BSD-3); fast-impl survey found no accelerator this kernel lacks; Morton/Hilbert reorder closed by convergent evidence (V18 + backlog line 60 + "gains likely small, dilute at tess 3"). |
| 8 | Real p2rank pattern: `consume=POINTS` at 16t / tess 3 | JMH `-prof gc`, AREA vs POINTS (the one untested measurement-gate cell) | **No headroom — last gap closed.** POINTS ≈ AREA in both wall-clock and allocation at every tess/thread cell; the zero-copy `surfacePointsXYZ()`/`PackedSurfaceAccess` path delivers the point set for free. 16t scales near-linearly. |

---

## The throughline (why there's no tess-2/3 win left)

Triangulated by in-tree experiment **and** external literature:

1. **The occlusion scan is a dead well.** The direction-major scan (last-occluder hint + per-direction
   early-exit) is simultaneously near-minimal in scalar dot-tests, 4-wide vectorized, and sequential in
   memory. Every spatial-narrowing alternative (bitmask, DCLM, region-hint) trades that for gather + a
   per-neighbor cap→direction mapping and loses in wall-clock even when it wins on test count.
   **Decisive lens: a scalar-test-count win must beat a 4-wide-SIMD sequential baseline (~4 scalar tests
   ≈ 1 unit of wall-clock; gather is strictly worse than sequential).**
2. **The build is tapped.** Distance pass dominates (90%) and is already SIMD + locality-optimized; the
   candidate waste is inherent uniform-grid geometry; the tighter-grid lever is a documented negative (A2).
3. **Float can't widen the tess-3 SIMD.** The 8-lane float scan's Vector-API intrinsics deopt to boxing
   under concurrency at tess ≥ 3 (C9) — not fixable in source. Float build (C1) helps only at tess 2.
4. **Different-algorithm exact SASA is area-only.** Analytic Gauss-Bonnet and power-diagram both emit no
   point set and differ ~21% per-atom from the sampled surface — they cannot serve p2rank's point consumer.

**Meta-lesson worth carrying forward:** *counting headroom ≠ wall-clock headroom when the baseline is
vectorized*, and *a Vector-API speedup that holds single-threaded can invert under concurrency if the
intrinsics deopt* — always validate vectorized variants at the deployment thread count and watch alloc-rate.

---

## Recommendation (unchanged)

- **Bit-exact default:** `DistinctPackedNumericalSurfaceV3`.
- **Tolerance / float (tess 2 only):** `FloatNumericalSurfaceV2`. Do **not** use float surfaces at tess ≥ 3
  with > 1 thread (C9 collapse). Use double `V3` at tess ≥ 3 — it scales near-linearly to 16 threads.

## What's left (out of the tess-2/3 single-protein point-set scope — needs a scope decision)

- **C4** whole-dataset batch / GPU throughput (a distinct offering).
- **C5** intra-protein parallelism for single-large-protein latency.
- A native **power-diagram** O(N log N) path for an **area-only** consumer (not p2rank).
- A changed **operating point / output contract** (new tessellation, area-only mode) reopens specific leads.

## Artifacts produced

- **Instrumentation (test scope, reusable):** `BitmaskFeasibilityTest` + `BitmaskFeasibilityScan` +
  `CapDirectionLut`, `DclmFeasibilityTest` + `DclmFeasibilityScan`, `AnalyticGateTest`, `BuildProfileTest`,
  `FloatBuildToleranceGateTest`.
- **Benchmark logs:** `autoresearch/results/phase2-float-tess3.txt`, `phase6-floatbuild-tess3.txt`,
  `phase8-points-path.txt`, `phase7-deep-research.json`.
- **Docs updated:** `docs/performance-lessons.md` (lesson 5 / C9 confirmed + mechanism; "not worth doing"
  gains the LUT-bitmask entry), `docs/optimization-backlog.md` (§1b closed, C9 confirmed), `CLAUDE.md`
  (float tess ≥ 3 guidance), this report + `autoresearch/LOG.md` + `autoresearch/KICKOFF.md`.
- **Memory (cross-session):** `lut-bitmask-scan-track` (closed-negative), `scan-is-a-dead-well`.

## To resume

Re-arm `/loop execute the autoresearch kickoff` only after a scope/contract change. Read this report and
the `scan-is-a-dead-well` memory first; do **not** re-try any closed lead (scan bitmask/DCLM/region-hint,
float scan at tess ≥ 3, analytic/power-diagram for points, tighter grid, Morton relabeling, float build at
tess 3).
