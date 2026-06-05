# Autoresearch log — FasterMolecularSurface (tess 2 & 3)

Append-only phase journal. Each phase: hypothesis · cheap kill-experiment · what was done · measurement
(with load-before) · verdict · meta-lesson · commit. Newest at the bottom.

See `KICKOFF.md` for the mission and operating rules (it is rewritten each phase to reflect the latest
state and lead ranking).

---

## Phase 0 — kickoff (setup)

- Created the `autoresearch` branch and this folder. Authored `KICKOFF.md` (self-updating mission +
  operating rules + prioritized leads, scoped to tess 2 & 3).
- Inherited state: bit-exact default `DistinctPackedNumericalSurfaceV3`; tess-2 float
  `FloatNumericalSurfaceV2`. Wins so far: A6 (SIMD build), C1 (float build). Negatives: A2, A7 (double +
  float), C3, naive A1, A4.
- Biggest open lead: the LUT-bitmask / fully-buried-atom prize (61% of atoms / 57% of scan time at tess 2).
- Next: Lead #1 — prototype the cheap cap→direction bitmask + instrument its cost before building a variant.

---

## Phase 1 — Lead #1 LUT-bitmask: the floor / budget kill-experiment (count-based, load-immune)

Machine was loaded at kickoff (1-min load ~25) — benchmarking blocked, but this phase is pure
counting, which is load-immune.

**Hypothesis.** A neighbor-major scan that ORs a precomputed cap→direction bitmask per neighbor and
stops once all directions are confirmed buried can beat the current direction-major early-exit scan,
because it gets a true *atom-level* early-exit on the 61% fully-buried atoms (57% of tess-2 scan tests).

**Why a floor measurement first (cheap kill-experiment).** Any real cap mask only *adds* cost over a
free, perfect mask. So measure the idealized lower bound before designing the mask:
- **Perfect-mask neighbor-major exact-dot-tests == `buriedDirections`.** Each buried direction is
  confirmed exactly once by its first occluding neighbor; survivors are never tested; confirmed
  directions are never re-tested. So the floor on exact dot-tests for the whole scan is simply the
  number of buried directions — already counted by the existing instrumentation.
- The realistic win therefore lives in the gap `currentNeighborTests − buriedDirections`, and the
  cap-mask *build* cost (per neighbor processed, until the atom is confirmed fully buried) must fit
  inside that gap. So also measure **neighbors-processed-until-fully-buried** (the mask-build budget)
  and the **no-mask** neighbor-major cost (the upper bound if the mask is useless).

**Decision rule.** If `buriedDirections` is not materially below `currentNeighborTests` (i.e. the
floor barely beats today), the idea is dead regardless of mask cleverness — record as a documented
negative like A2/A7 and move to Lead #2. If the gap is large, design a concrete cheap mask and check
it fits the per-neighbor budget.

Built `BitmaskFeasibilityTest` (+ `BitmaskFeasibilityScan`, `CapDirectionLut`) — new test classes,
existing instrumentation untouched so its numbers stay reproducible. Measures at tess 2 AND tess 3,
sweeping LUT resolutions (cube-face grid × cosθ bins), and includes a per-(neighbor,direction)
soundness sweep (every truly-buried pair MUST be in the LUT mask, else not bit-exact).

**Run:** `./gradlew scorecard --tests '*BitmaskFeasibilityTest'`.

**Results (corpus, 14555 atoms).**

| | tess 2 | tess 3 |
|---|---|---|
| fully-buried atoms | 61.0% | 52.8% |
| their share of current scan tests | 57.0% | 44.1% |
| current (direction-major SIMD, scalar-equiv tests) | 4,856,041 | 12,417,807 |
| perfect-mask floor (= buriedDirections) | 11.9% | 18.0% |
| no-mask neighbor-major | 161.8% | 244.1% |
| real LUT, grid16/bins32, exact tests | **18.6%** | **28.0%** |
| ↳ false-positive rate | 85.0% | 85.9% |
| ↳ candidates examined / neighbor | 8.1 | 31.1 |
| soundness violations (all configs) | **0** | **0** |

The LUT is provably sound (0 violations over the full corpus after fixing a face-parametrization sign
bug), so a confirm-each-candidate bitmask scan *would* be bit-exact.

**Verdict: NEGATIVE — the LUT-bitmask / fully-buried-atom prize is infeasible. Lead #1 closed.**

Why, decisively:
1. **The ceiling is too low and unreachable.** The perfect-mask floor is only 11.9% of current
   *scalar* tests, but the current scan is already 4-wide SIMD (wall-clock ≈ 25% of its scalar count
   with sequential access). So even a *free, perfect* mask is ~2× on an op that is ~34% of tess-2 CPU
   at 16t → ≤~17% total — and a free perfect mask does not exist.
2. **A cheap mask is far too loose.** 85–92% false positives even at a 0.4–1.2 MB LUT; the rate barely
   moves from grid6→grid16 (91→85%) while the table grows 16×. Real buried caps are small and finely
   shaped; quantized (cell × cosθ-bin) supersets can't resolve them without per-direction precision,
   which is just the exact test again. So the realistic 18.6%/28% scalar-test count is done with
   scattered gathers (random `ddx[d]` by popcount) + a per-neighbor sqrt + 3 divides mask-build + a
   large LUT thrashing shared L3 at 16–32 threads — slower than the SIMD baseline in wall-clock.
3. **The approximate pure-OR-accumulator form is dead too.** At 85% false positives, OR-ing supersets
   marks nearly everything buried → drastic surface-point under-count, far outside the float tolerance
   envelope. No tolerance variant rescues it.

Root cause / meta-lesson: the existing per-direction early-exit + last-occluder hint is already a
near-optimal *scalar-test minimizer* AND it vectorizes; the fully-buried atoms are cheap precisely
because their directions bury fast. The "57% of tests on 61% of atoms" headroom is real in test-count
but is dominated by SIMD width + memory locality, which a gather-heavy neighbor-major scan with a loose
LUT cannot beat. **Counting headroom ≠ wall-clock headroom when the baseline is vectorized.** This is
the same trap A1/A4 fell into; now quantified end-to-end so it is closed, not merely suspected.

Kept the three instrumentation classes as the durable evidence (like `ScanInstrumentationTest`).
No new production variant scaffolded (correctly — the kill-experiment killed it before that cost).

**Next: Lead #2 — float-path scaling at tess 3 / 16 threads (now the top lead). Needs an idle box for
the timing half; the count/alloc-probe half can start anytime.**

---

## Phase 2 — Lead #1 float tess-3/16t scaling (IN PROGRESS — gated on idle box)

Load-immune de-risk done first (box was at 1-min load ~1.5–3, just over the 1.5 gate, so the timing
half is deferred to the next idle window).

**Code read (`Float256WeightedDedupOcclusionScan`).** Its `FloatVector`s (`vx/vy/vz/vth/dot`) are created
and discarded inside the inner per-direction loop — i.e. already method-local, the SAME shape as the
double `Vectorized256WeightedDedupOcclusionScan` that does NOT collapse. Implication: the C9 "restructure
so the vectors stay method-local and scalar-replace" lever (lesson 3) is likely a near no-op — they're
already local. So IF `-prof gc` confirms a float alloc blowup at tess 3, the root cause is the JIT
failing to scalar-replace `FloatVector.SPECIES_256` (8-lane) where it succeeds for `DoubleVector.SPECIES_256`
(4-lane), which is NOT fixable by a trivial source restructure — would need a different formulation
(e.g. avoid `VectorMask`/`firstTrue` boxing, or a reduce-based verdict).

**Next (idle box):** `-prof gc` alloc.rate.norm for `FloatNumericalSurfaceV2` vs `V3` at tess 3,
`-t 1/4/8/16`, `-p consume=AREA`. Confirm/refute: (a) does float wall-clock collapse begin at tess 3
(not just tess 4)? (b) does alloc.rate.norm scale with thread count (the boxing signature)? Then verdict.

### Phase 2 — RESULTS (benchmark ran in the idle window, load-before 1.13–1.83; alloc numbers load-immune)

`FLOAT_V2` vs `DISTINCT_PACKED_V3` @ tess 3, consume=AREA, `-prof gc` (full log:
`autoresearch/results/phase2-float-tess3.txt`):

| threads | FLOAT_V2 ms/op | V3 ms/op | float÷double | FLOAT_V2 alloc B/op | V3 alloc B/op |
|---|---|---|---|---|---|
| 1  | 21.12  | 21.60 | 0.98× (float faster) | 13.6 MB | 13.7 MB |
| 4  | 163.50 | 22.85 | **7.2×**  | **1.084 GB** | 13.7 MB |
| 8  | 302.04 | 24.26 | **12.5×** | **1.084 GB** | 13.7 MB |
| 16 | 611.15 | 26.74 | **22.9×** | **1.084 GB** | 13.7 MB |

**Verdict: the float-scan collapse begins at TESS 3 (not just tess 4), and the C9 boxing hypothesis is
CONFIRMED.** At 1 thread the `FloatVector`s cost nothing extra (13.6 MB/op, marginally faster than
double); at ≥2 threads per-op allocation jumps **80× to 1.084 GB/op**, alloc rate ~27 GB/s (matching
C9's figure), and wall-clock collapses 7–23×. Double V3 is flat (21.6→26.7 ms, 1→16t — near-linear
throughput scaling). The `alloc.rate.norm` signal is load-immune and the interleaved float/double ratio
controls for contamination, so the result is robust.

**Mechanism nailed (cheap EA probe, single thread):** `-XX:-DoEscapeAnalysis` does NOT change the float
allocation at 1 thread (13.56 MB/op both ways). So the low single-thread allocation is NOT generic
escape-analysis scalar-replacement — it's **Vector-API intrinsics** keeping the float vectors
register-resident. The multithread collapse is therefore the **Vector-API float (8-lane) intrinsics
failing to fire under concurrency**, falling back to the boxing Java implementation (~1.08 GB/op of
`FloatVector` objects). This is NOT fixable by EA, and NOT by keeping the vectors method-local (the
Phase-2 code read confirmed they already are — same shape as the non-collapsing 4-lane double scan).

**Fixability: NONE in the Vector-API float path.** No source restructure makes the JVM re-intrinsify;
a scalar-float scan would drop the only reason to go float (8-wide SIMD). So `FloatNumericalSurfaceV2`
cannot become a clean tess-3 win. Lead resolved as a negative-diagnosis.

**Action / takeaways:**
- Confirms + sharpens the existing guidance: float surfaces are **tess-2-only** — tess 3 already
  collapses, starting at just 2–4 threads. Updated CLAUDE.md ("tess ≥ 4" → "tess ≥ 3"), perf-lessons
  lesson 5 / C9 (hypothesis → CONFIRMED, with the EA-off mechanism), backlog C9.
- Strong positive corollary, now data-backed: **double `DistinctPackedNumericalSurfaceV3` is the correct
  tess-3 recommendation** — it scales near-linearly to 16 threads where float is catastrophic.
- Meta-lesson: a Vector-API speedup that holds single-threaded can *invert* under concurrency if the
  intrinsics stop firing — always validate vectorized variants at the deployment thread count, and
  diagnose alloc-rate (boxing) not just wall-clock. `-prof gc` alloc.rate.norm is the load-immune probe.

**Next: Lead #1 is now A5 (DCLM second dot-lattice, bit-exact, tess-3 oriented) or B1 (analytic SASA,
tolerance). Both bit-exact-scan wells look tapped; remaining upside is a different algorithm (A5) or the
build. See KICKOFF re-rank.**

---

## Phase 3 — Lead #1 A5 (DCLM second dot-lattice): survivor-angle kill-experiment (count-based, load-immune)

Box loaded (1-min 5.9) — counting is load-immune.

**Framing.** A5/DCLM (bin the fixed tessellation dots so each neighbor tests only the dots in its cap)
is the SAME neighbor-major spatial-narrowing family that Phase 1 closed at tess 3 (real-LUT neighbor scan
= 28% of current scalar tests but gather-heavy + 86% false positives → loses to the SIMD direction-major
baseline). The DCLM lattice is just a different spatial index than Phase 1's cube-face LUT; it doesn't
escape the gather/false-positive penalty. So most of A5 is already decided.

The ONE angle Phase 1 didn't isolate, and A5's strongest case: **surviving directions.** The current
direction-major scan resolves *buried* directions fast (hint + early-exit, ~1 test via hint 50% of the
time), but a *surviving* direction has no occluder, so it must be tested against ALL `numNeighbors`
neighbors to confirm survival. DCLM's spatial bin could cut that to only the neighbors whose cap is near
the direction. (C3/region-hint already tried — and lost — on the *buried* side; the survivor side is
distinct.)

**Decision rule.** Measure (a) the share of total scan tests spent on surviving directions, and (b) per
survivor, the DCLM candidate count (neighbors whose cap-LUT includes the direction) vs `numNeighbors`.
If survivor cost is a small share, OR survivors already have ~all neighbors as candidates, A5 has no
opening and is closed by Phase 1 + this. If survivor cost is large AND candidates are few, there may be a
narrow win worth a real DCLM variant — but it must still clear the Phase-1 SIMD-baseline bar.

**Results (`DclmFeasibilityTest`, finest LUT grid16/bins32, corpus).**

| | tess 2 | tess 3 |
|---|---|---|
| surviving directions | 5.4% | 5.4% |
| survivor tests as share of all scan tests | 16.7% | 25.2% |
| current tests per survivor (all neighbors) | 24.5 | 24.5 |
| DCLM candidates per survivor | 0.6 | 0.6 |
| survivor saving (DCLM) as share of ALL scan tests | 16.3% | 24.6% |

**Verdict: A5/DCLM CLOSED (negative).** The survivor angle is real and striking — survivors are 25% of
tess-3 scan tests (each scans all ~24 neighbors, no early-exit), and a sound spatial bin would cut that
~40× (24.5 → 0.6 candidates). BUT this is NOT headroom beyond Phase 1: Phase 1's neighbor-major real-LUT
total was already 28% of current scalar tests at tess 3, a figure that ALREADY excludes most survivor
tests (survivors are rarely candidates). So the survivor saving is *part of* Phase 1's 28%, not additional.
Whichever way DCLM flows (neighbor-major OR direction-major-with-per-atom-neighbor-bin), it still pays the
per-neighbor cap→cells mapping (sqrt + divides + marking) and ends in gather-based candidate tests — the
exact components Phase 1 showed lose to the 4-wide-SIMD sequential baseline. A cheap coarse pre-filter
(octant histogram) can't isolate the 0.6-candidate survivors; that selectivity required the fine LUT.
DCLM's classic win is vs a NAIVE O(D×N) scan, a premise our SIMD + hint + early-exit baseline already
defeats. Kept the instrumentation (`DclmFeasibilityScan`/`DclmFeasibilityTest`) as evidence.

Meta-lesson (reinforces Phase 1): every attempt to beat the *scan* loses because the scan is already
SIMD + sequential + early-exit + hint; any spatial-narrowing alternative trades that for gather + a
per-neighbor cap mapping. Survivors are the scan's worst case but can't be skipped without first paying
to identify them.

---

## State of the search (after 3 consecutive closes: §1b bitmask, float tess-3/C9, A5 DCLM)

Three leads closed with no win. The unifying finding: **the bit-exact SCAN is a dead well.** The
direction-major SIMD scan with last-occluder hint + per-direction early-exit is at once (a) near-minimal
in scalar dot-tests, (b) vectorized 4-wide, and (c) sequential in memory. Every alternative explored —
neighbor-major bitmask (§1b), DCLM spatial dot-lattice (A5), region-binned hints (C3, earlier) — replaces
that with gather-based access + a per-neighbor cap→direction mapping, and loses in wall-clock even when it
wins in scalar-test count. **The decisive lens is always: a scalar-test-count win must beat a 4-wide-SIMD
sequential baseline, so ~4 scalar tests ≈ 1 unit of baseline wall-clock, and gather is strictly worse than
sequential.** Float precision (8-wide) is the one thing that could widen SIMD, but its Vector-API
intrinsics deopt to boxing under concurrency at tess ≥ 3 (Phase 2) — so float is tess-2-only.

**Pivot. Stop attacking the scan.** Remaining productive directions, re-ranked:
1. **B1 — analytic per-atom SASA (tolerance).** Sidesteps sampling entirely — a fundamentally different
   algorithm, not a scan refinement, so it escapes the dead well. Cheap load-immune gate FIRST: implement
   only the analytic *area* (per-atom inclusion-exclusion of cap overlaps / Gauss-Bonnet arcs), compare its
   per-atom + total area to sampled tess-2/3 over the corpus — is the gap within p2rank's tolerance
   (total ≤ 1e-4 rel, per-atom ≤ 2%)? If yes, it's a candidate fast opt-in surface; if the gap is too big,
   closed cheaply. (No fast impl needed for the gate, just a correct one.)
2. **The BUILD, not the scan.** A6 (SIMD build) was the last real win and the build is the larger share at
   tess 2 / 16t (bandwidth-bound). Re-profile where build time actually goes now (grid construction vs
   distance pass vs neighbor materialization) before assuming it's tapped.
3. **Generate new ideas** outside the scan/build dichotomy if 1–2 stall.

---

## Phase 4 — Lead #1 B1 (analytic per-atom SASA): cheap gate (load-immune, ran on idle box)

**Hypothesis / lead.** Analytic SASA sidesteps sampling — a different algorithm that escapes the dead
scan well; could be a fast opt-in area path. Gate first: is the gap between analytic and sampled tess-2/3
within p2rank tolerance?

**Method.** Analytic SASA is the converged limit of the dot-sampled area as tessellation refines, so
tess-5 (2562 directions) is a faithful "analytic truth" proxy. `AnalyticGateTest` compares sampled
tess-2/3/4 area to tess-5 over the corpus (area-only — no point keying — so it scales). 8.9 s, load 1.0.

**Results vs tess-5 (≈analytic):**

| tess | total area rel err | per-atom rel err mean | per-atom max |
|---|---|---|---|
| 2 | 3.04e-3 | 21.1% | 455% |
| 3 | 1.63e-3 | 8.0% | 141% |
| 4 | 5.95e-4 | 2.8% | 61% |

p2rank tolerance is total ≤ 1e-4, per-atom ≤ 2%. Sampled tess-2/3 is **16–30× over on total and 4–10×
over per-atom** vs analytic.

**Verdict: B1 CLOSED (negative for this project's goals).** Three converging reasons:
1. **Not interchangeable with the sampled surface.** Analytic and sampled differ by 0.3% total / 21%
   per-atom at tess 2 — genuinely different numbers, well outside any drop-in tolerance. So analytic
   cannot satisfy the bit-exact OR the tolerance-vs-sampled-oracle contract; it's a different output.
2. **Area-only, no point cloud.** p2rank consumes surface POINTS (feature extraction); analytic SASA
   yields per-atom scalars only. It can't serve the primary consumer.
3. **Almost certainly slower.** Exact per-atom SASA is O(neighbors²) pairwise-arc trig (atan2/acos),
   no SIMD; for ~24 neighbors that's ~hundreds of transcendental ops/atom vs the dot scan's vectorized
   tests — and the dot scan is already fast. No speed win even for the area-only niche.

Interesting (non-actionable) corollary: dot-sampled per-atom areas are quite inaccurate vs true SASA
(21% mean at tess 2) — inherent to Shrake-Rupley sampling, presumably fine for p2rank (point-based,
tuned on the sampled surface). The total area converges much faster (errors cancel) than per-atom.

Kept `AnalyticGateTest` as evidence. No variant built (gate killed it cheaply). 4th consecutive close.

**Next: Lead #2 — re-profile the BUILD (not the scan).** The scan is a dead well (Phases 1–3) and the
two alt-output ideas (float tess-3 = C9; analytic = B1) are closed. The build is the remaining lever:
A6 (SIMD build) was the last real win, the build is the larger share at tess 2 / 16t (bandwidth-bound),
and it has NOT been re-profiled since. Start by measuring where build time actually goes now (grid
construction vs the distance pass vs neighbor materialization), load-immune where possible.

---

## Phase 5 — Lead #1 re-profile the BUILD (load-immune phase decomposition)

The scan is a dead well (Phases 1–3) and both alt-output ideas are closed (float tess-3 = C9, analytic
= B1). The build (neighbor-list construction) is the only remaining lever — A6/SIMD build was the last
win, the build is the larger share at tess 2 / 16t (bandwidth-bound, lesson 12), and it hasn't been
re-profiled since A6 shipped.

**Hypothesis.** After A6 vectorized the distance pass (was ~47% of CPU), the remaining build cost has
shifted. Suspect the `edges[]` intermediate in `SimdDistanceCellGridNeighborList`: every surviving pair
is written to a growing `int[]` and then re-scattered into the CSR `adj[]` — double the adjacency write
traffic, which matters in a bandwidth-bound build.

**Cheap kill-experiment.** A test-scope instrumented CLONE of the V3 build (frozen class untouched) with
`nanoTime` around each phase — (A) CellGrid construction, (B) expandedR, (C) SoA repack, (D) distance
pass, (E) edges→CSR — run over the corpus (warmup + measured), reporting the % breakdown. Single-thread
shows the compute distribution; the bandwidth bottleneck is a 16t effect, but this localizes the target.

**Results (V3 SIMD build clone, corpus, single-thread, 5 measured passes).**

| phase | share | ms/pass |
|---|---|---|
| grid construction | 0.9% | 0.11 |
| expandedR | 1.6% | 0.19 |
| SoA repack | 0.8% | 0.10 |
| **distance pass** | **90.2%** | **10.7** |
| edges→CSR | 6.5% | 0.78 |

distance-pass candidate tests = 1.93M for 308k surviving pairs → **6.27 candidates/survivor (16% survive)**.

**Verdict (bit-exact build): tapped within the current grid approach.** The distance pass dominates (90%);
the `edges→CSR` intermediate I suspected is only 6.5% (not worth attacking). The pass already does the
exact pairwise prune and is SIMD-vectorized (A6). The 6.27× candidate waste (only 16% survive) is the
inherent geometry of a uniform grid (cell = max-pair-cutoff = 2·(maxVdw+solvent); cubic stencil vs
spherical cutoff + radius variation), and the one lever to cut it — a tighter grid — is the A2 documented
negative. No obvious NEW bit-exact build win single-thread.

**But — a promising untested follow-up surfaced:** C1 = `DevSurfaceV24FloatBuild` (float build + DOUBLE
scan) cuts the bandwidth-bound distance pass's coordinate traffic, won +1.6–2.8% at tess 2, and keeps the
double scan → it has NO C9 float-scan collapse at tess 3. It has never been benchmarked at tess 3. If it
wins there (it should: build traffic cut + no C9), it's a promotable tess-3 tolerance variant — distinct
from FloatNumericalSurfaceV2 (which collapses at tess 3). **Next: benchmark V24 vs V3 at tess 3, -t 1/16
(idle box), + a tess-3 tolerance gate before any speed claim.** Grabbing the current idle window for the
benchmark.

---

## Phase 6 — C1 float-build (V24) at tess 3: tolerance gate + benchmark (IN PROGRESS)

Testing whether `DevSurfaceV24FloatBuild` (float-precision SIMD neighbor build + DOUBLE scan, idea C1) is
a promotable tess-3 win: it cuts the bandwidth-bound distance pass's coordinate traffic (the 90%-of-build
hot phase, Phase 5) AND keeps the double scan so it has no C9 collapse at tess 3.

**Tolerance gate (load-immune, `FloatBuildToleranceGateTest`, vs V3 over the corpus): PASS — and stronger
than required.** Total area rel err, per-atom rel err, and point-set symdiff are ALL exactly 0.000e+00 at
tess 2 AND tess 3. The float build produces a BIT-IDENTICAL neighbor set to the double build on every
corpus structure — the float vs double `d²<sumR²` test never disagrees at these atom spacings. So V24 is
empirically exact on the corpus (still nominally "tolerance" since a pathological boundary pair could flip,
but zero error observed). Accuracy is a non-issue; the question reduces to speed.

**Benchmark (idle box, running): V24 vs V3 @ tess 3, -t 1/16, consume=AREA, -f 3.** Results pending.

### Phase 6 — RESULTS (idle box, load-before 1.16/1.95)

V24 (float build + double scan) vs V3 @ tess 3, consume=AREA, -f 3 (raw:
`autoresearch/results/phase6-floatbuild-tess3.txt`):

| threads | V24 ms/op | V3 ms/op | Δ |
|---|---|---|---|
| 1  | 21.178 ± 0.264 | 21.537 ± 0.134 | V24 ~1.7% faster (CIs barely touch) |
| 16 | 26.923 ± 1.808 | 26.540 ± 0.147 | V24 ~1.4% slower, V24 noisy (bimodal) |

**Verdict: NEGATIVE for tess-3 promotion.** A marginal single-thread gain (~1.7%) that vanishes — turns
slightly negative and noisy — at 16t, the deployment regime. The mechanism is the key insight: **at tess 3
the build is a SMALLER fraction of total** (the double scan over 162 directions dominates), so a build-only
optimization (C1/float build) dilutes to nothing. C1 wins at tess 2 *precisely because* the build is the
larger share there (~47%); that premise doesn't hold at tess 3. V24 stays a tess-2 variant; no tess-3
promotion. (Tolerance gate had passed bit-identically, so accuracy was never the issue — purely a
speed-share dilution.) Kept `FloatBuildToleranceGateTest` as evidence.

---

## State of the search — round 2 (after 6 phases, 0 promotable wins)

The tess-2/3 surface is **well-optimized**; this round found no shippable win, but tightened the map:
- **Scan: dead well** (Phases 1–3) — SIMD + hint + early-exit is near-optimal; every spatial-narrowing
  alternative loses to it (gather vs sequential-SIMD).
- **Build: tapped** (Phase 5) — distance pass is 90% and already SIMD (A6); the candidate waste (6.27×,
  16% survival) is inherent uniform-grid geometry and the tighter-grid fix is the A2 negative.
- **Float precision** is the only SIMD-widening lever, but: float SCAN collapses at tess ≥ 3 (C9, Phase 2),
  and float BUILD (C1) is a tess-2-only win — it dilutes at tess 3 where the scan dominates (Phase 6).
- **Alt algorithms:** analytic SASA (B1) closed — different output, area-only, slower.

**The structural reason there's no tess-3 win left:** the double scan is the dominant, already-optimal
tess-3 cost, and float can't widen it (C9). The remaining ideas are out-of-scope (batch/GPU C4–C5,
deployment-shaped) or need genuinely new external input.

**Next: the only productive moves left are (a) the backlog §3 deep-research prompt for novel approaches,
or (b) accept the surface is at a strong local optimum and stop.** No more cheap in-tree kill-experiments
remain for tess 2/3 single-protein throughput.

---

## Phase 7 — generative: deep-research for novel approaches (backlog §3, scoped to the still-open questions)

All cheap in-tree leads are exhausted (Phases 1–6, 0 promotable wins; surface at a strong local optimum).
Per the kickoff guardrail, the productive move is generative: source approaches NOT yet tried. The §3
prompt's Q1–Q4 (buried-atom test, bitmask, DCLM, analytic) are now empirically CLOSED this session, so
this run is scoped to the open external questions:
- Q5: power-diagram / weighted-alpha-complex exact SASA O(N log N), and the hard dependency — is a robust
  3D weighted/regular Delaunay usable from Java (pure-Java or JNI)?
- Q6: neighbor-list reuse (Verlet skin / sort-and-sweep) — applicability given coords differ per protein.
- Q7: what do the FASTEST open-source SASA impls (FreeSASA, RustSASA, GROMACS, MSMS, NanoShaper, dr_sasa,
  PowerSASA, …) do that an already-dedup'd SIMD cell-list Shrake-Rupley kernel does NOT have — especially
  spatial reordering (Morton/Hilbert) for the bandwidth-bound build, whole-atom culling, or batch reuse.

Hypothesis: most likely outcomes are (a) confirmation the kernel is near-SOTA for sampled SASA, or (b) one
new accelerator — the prime candidate being atom reordering along a space-filling curve for the
bandwidth-bound distance pass (note: V18 SortedCoords already tried coordinate sorting — check its result
before treating Morton as novel). Each returned idea gets a cheap in-tree kill-experiment, held to the
Phase-1 lens (beat 4-wide-SIMD sequential) and the tess-3 share insight (build-only wins dilute; C9 blocks
float). Running the `deep-research` skill now.

### Phase 7 — RESULTS (deep-research, 5 angles, 19 sources, 25 claims adversarially verified 24/1)

Full report: `autoresearch/results/phase7-deep-research.json`. Triage against the baseline (must produce a
POINT SET; must beat 4-wide-SIMD sequential; build-only wins dilute at tess 3):

- **Power-diagram / weighted-alpha-complex exact SASA (dSASA 2024) — DEAD, independently confirmed.**
  Verified (3-0): it produces ONLY scalar area + analytical derivatives, emits NO surface point set —
  categorically fails p2rank's point-set need. This re-derives the B1 closure (Phase 4) via a totally
  different method. (Also: robust 3D regular/power Delaunay exists only in native C++ (CGAL GPLv3 / Geogram
  BSD-3) — JNI-only from the JVM — but moot since the output is area-only.)
- **Fast-SASA impl survey (Q3) + neighbor reuse (Q4) — no distinguishing accelerator found.** No claim
  survived verification identifying a trick (whole-atom culling, reuse, reordering, adaptive sampling) that
  this dedup'd SIMD cell-list Shrake-Rupley kernel lacks.
- **Spatial reordering / Morton-Hilbert (Q1) — the one "maybe", CLOSED by convergent evidence.** The
  mechanism is real (LAMMPS, GROMACS reorder atoms) but the verified caveat is decisive: the gains
  (LAMMPS 0–3.4×, SFC 25–75% cache misses / up to 50%) come from large DRAM-resident many-timestep
  workloads; my build is one-shot, cache-resident (tens of KB), ALREADY cell-bucket-ordered — "the
  working-set mechanism that drives those gains is largely absent ... expect small gains," and it dilutes
  at tess 3. Convergent in-tree evidence closes it without building a variant: (i) **V18 SortedCoords
  already** makes the candidate-coordinate reads sequential (cell-sorted SoA — it eliminated exactly the
  random gather Morton would target; it was the compute champion); (ii) **backlog line 60** already records
  "Morton/Z-order atom relabeling — bit-exact but modest, and the output inverse-permutation/point re-sort
  may eat the locality gain." Three independent lines agree → not worth a variant. (Also perf counters are
  blocked on this box, so the cache-miss delta the research wants couldn't even be measured directly.)

**Verdict: deep-research surfaced no new viable lever; it independently confirms the search is complete.**

---

## SESSION CONCLUSION — the tess-2/3 single-protein surface is at a strong local optimum

Seven phases, zero promotable wins, but a sharply tightened map and several long-open questions closed
with hard data. The structural picture is now complete and triangulated (in-tree experiments + external
deep-research agree):

1. **The occlusion SCAN is a dead well.** Direction-major SIMD + last-occluder hint + per-direction
   early-exit is near-minimal in scalar tests AND vectorizes AND is sequential. Every spatial-narrowing
   alternative (bitmask §1b/Phase 1, DCLM A5/Phase 3, region-hint C3) trades that for gather + a
   per-neighbor cap→direction mapping and loses in wall-clock. Decisive lens: a scalar-test-count win must
   beat a 4-wide-SIMD sequential baseline (~4 scalar ≈ 1 unit wall-clock; gather is strictly worse).
2. **The BUILD is tapped.** Distance pass is 90% (Phase 5), already SIMD (A6) and locality-optimized
   (V18 cell-sorted coords); the 6.27× candidate waste is inherent uniform-grid geometry (tighter grid =
   A2 negative); Morton relabeling is modest (Phase 7). Float build (C1) helps only at tess 2 — it dilutes
   at tess 3 where the scan dominates (Phase 6).
3. **Float can't widen the tess-3 SIMD.** The float (8-lane) scan's Vector-API intrinsics deopt to boxing
   under concurrency at tess ≥ 3 (C9 CONFIRMED, Phase 2: 13.6 MB/op → 1.084 GB/op, 7–23× slower); not
   fixable in source. Float surfaces are tess-2-only.
4. **Different-algorithm exact SASA is area-only.** Analytic Gauss-Bonnet (B1/Phase 4) and power-diagram
   (Phase 7) both emit no point set and differ ~21% per-atom from the sampled surface — they cannot serve
   p2rank's point consumer.

**The recommended surfaces stand:** bit-exact default `DistinctPackedNumericalSurfaceV3`; tess-2 float
`FloatNumericalSurfaceV2`. No production change this session.

**Only directions left are out of the tess-2/3 single-protein scope:** whole-dataset batch / GPU
throughput (C4), intra-protein parallelism for single-large-protein latency (C5), or a native power-diagram
path for an AREA-ONLY consumer (not p2rank). These need a scope/requirements decision from the human, not
another in-tree kill-experiment.

**Stopping the autoresearch loop here** (kickoff stop condition met: leads exhausted + several new ideas
generated and tested, all documented). Re-arm `/loop execute the autoresearch kickoff` if a new direction
or a scope change (batch mode, area-only consumer, a different operating point) opens up.

---

## Phase 8 — the real p2rank pattern: consume=POINTS at 16 threads (the one untested measurement-gate cell)

Loop re-armed by the user after the Phase-7 conclusion → pursue the one genuinely UNMEASURED, in-scope gap
(§2 measurement-gate item 6): all this session's benchmarks used consume=AREA, but p2rank consumes the
POINT SET. The production table already shows POINTS ≈ AREA at tess 2 single-thread (V2: 13.05 vs 13.08 ms)
— zero-copy materialization is free per-call. UNTESTED: POINTS at **16 threads** (the GC-pressure regime the
flat store `DistinctFlatSurfacePointStoreV2` / `PackedSurfaceAccess` exists for) and at **tess 3**.

**Hypothesis.** If the point-materialization path has headroom, it shows as POINTS ≫ AREA in wall-clock OR a
large `gc.alloc.rate.norm` jump at 16t (GC pressure across concurrent builds). If POINTS ≈ AREA with
controlled allocation, the flat store already handles it → no headroom (the last gap closes).

**Benchmark (idle box, running): V3 AREA vs POINTS @ tess 2 & 3, -t 1/16, -prof gc.** Results pending.

### Phase 8 — RESULTS (idle box, load-before 1.01/2.02): POINTS path is FREE — no headroom, last gap closed

V3 AREA vs POINTS, -prof gc (raw: `autoresearch/results/phase8-points-path.txt`):

| cell | AREA ms/op | POINTS ms/op | AREA alloc B/op | POINTS alloc B/op |
|---|---|---|---|---|
| tess2 -t1  | 12.914 | 12.756 | 10,607,665 | 10,607,706 |
| tess3 -t1  | 21.536 | 21.581 | 13,660,961 | 13,671,408 |
| tess2 -t16 | 16.274 | 16.171 | 10,605,905 | 10,605,802 |
| tess3 -t16 | 26.835 | 26.904 | 13,658,114 | 13,658,111 |

**Verdict: the point-materialization path has NO headroom.** POINTS ≈ AREA in BOTH wall-clock and
allocation-per-op at every tess/thread cell — consuming the full surface point set (the real p2rank
pattern, via the zero-copy `surfacePointsXYZ()`/`PackedSurfaceAccess` path) costs the same as computing
just the scalar area, and allocates nothing extra. The flat-store design (V19/`DistinctFlatSurfacePointStoreV2`)
achieves its goal: zero marginal cost for point delivery. 16t scales near-linearly (tess2 1t→16t:
12.9→16.3 ms = 1.26× latency for 16× threads; GC is ~213 ms of the run, not a bottleneck). This closes
the last untested measurement-gate cell (§2 item 6) and confirms the surface is at a local optimum even
under the real consumer pattern — not just the AREA proxy used in Phases 1–7.

---

## SESSION CONCLUSION — round 2 (after Phase 8: the real p2rank POINTS pattern also confirms the optimum)

Phase 8 (run after the user re-armed the loop) tested the one in-scope cell never measured — the actual
p2rank `consume=POINTS` pattern at 16 threads / tess 3 — and found it free (POINTS ≈ AREA, same allocation).
This removes the last "but we never measured the real pattern" caveat: the tess-2/3 single-protein surface
is confirmed at a strong local optimum under the consumer's real workload, not just the AREA proxy.

**There are no remaining in-scope (tess-2/3, single-protein, point-set) optimization leads.** Every
in-tree lever (scan, build, float, store/emit, sampling) and the external literature (deep-research) have
been exhausted and documented. Continuing to launch benchmarks would re-confirm the optimum — i.e.
manufacture churn, which the kickoff explicitly forbids.

**The genuine blocker is now a SCOPE decision (human's to make), not an experiment.** The productive
directions all require widening the scope: whole-dataset batch/GPU throughput (C4); intra-protein
parallelism for single-large-protein latency (C5); a native power-diagram path for an AREA-ONLY consumer;
or a changed operating point / output contract. Surfacing that decision to the user and stopping the loop.

---

## Autoresearch CONCLUDED (user decision, 2026-06-05)

Asked the scope question after Phase 8; user chose **"Stop — accept the optimum."** No scope change (no
batch/GPU C4, no latency parallelism C5, no area-only path). The effort is concluded.

**Final state.** The tess-2/3 single-protein point-set surface is at a confirmed strong local optimum
(8 phases, 0 promotable wins, all leads + the real POINTS consumer pattern measured and documented).
Recommended surfaces stand, unchanged: bit-exact default `DistinctPackedNumericalSurfaceV3`; tess-2 float
`FloatNumericalSurfaceV2`. No production code changed this entire effort — only test-scope instrumentation,
docs, and the autoresearch journal.

To reopen: re-arm `/loop execute the autoresearch kickoff` only after a scope/contract change (batch mode,
single-protein latency, an area-only consumer, or a new operating point). Read this SESSION CONCLUSION and
the `scan-is-a-dead-well` memory first; do not re-try any closed lead.
