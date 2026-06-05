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
