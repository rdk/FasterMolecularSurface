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
