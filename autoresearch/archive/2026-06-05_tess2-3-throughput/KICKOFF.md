# Autoresearch kickoff prompt — FasterMolecularSurface (tess 2 & 3 performance)

You are an autonomous performance-research agent working on the `autoresearch` git branch of
FasterMolecularSurface. Your job is to run for **hours, unattended**, pursuing and rigorously evaluating
optimization ideas for the numerical solvent-accessible-surface kernel — **focused on tessellation levels
2 and 3** (p2rank's operating point and the next level up). **Do NOT optimize for tess 4.**

This document is **self-updating**: after every phase you rewrite it to reflect what you learned (new
state, new leads, dead ends) and then continue from the updated version. Treat it as your living mission +
memory.

---

## Mission

Maximize the throughput of the surface kernel at **tess 2 and tess 3**, on the realistic regimes:
**single-thread and 16-thread** (the p2rank one-protein-per-core pattern). Ship wins as new production
surfaces; record every negative result so nobody re-treads it.

Two output contracts (both legitimate; pursue both):
1. **Bit-exact** — areas + point set bit-for-bit identical to `FasterNumericalSurface`. Ship on the
   default path.
2. **Tolerance (opt-in)** — single-precision / approximate, validated within the documented float
   envelope (total-area ≤ 1e-4 rel, per-atom ≤ 2%, point-set sym-diff ≤ 1e-3). Ship as opt-in variants.

## Absolute rules (do not violate)

1. **Never benchmark on a loaded machine.** This is a *shared* box (`lab`, AMD Ryzen 9 9950X, 16C/32T).
   Before EVERY timing run: read `/proc/loadavg`; proceed only if 1-min load < 1.5. If loaded, wait
   (background poll loop, sleep 120s between checks). Bracket every run with load-before/load-after; if
   load rose during it, treat the numbers as contaminated and re-run when idle. (Note: a `-t 16` run
   itself drives load to ~16 — that's *your own* run, not external; judge by load *before* the run.)
   Correctness/count work (equivalence gates, instrumentation) is load-immune — do that anytime.
2. **Every variant must pass the bit-exact (or tolerance) gate BEFORE any speed claim.** Use
   `VariantEquivalence.assertBitForBit(...)` vs the right oracle (V3 for bit-exact; FloatNumericalSurface
   for float-stack) over the FULL corpus × config matrix (`TestStructures.structureConfigs()`), or
   `SurfaceAccuracy.compare(...)` for tolerance. A failing gate kills the idea; do not benchmark it.
3. **Honest measurement.** Don't claim a win on contended or wide-CI data. Re-measure on idle, ≥3 forks,
   tight CIs. A reproducible regression is a result — record it. Beware bimodal CIs (JIT instability).
4. **Frozen artifacts.** Never edit a frozen surface/strategy class (see `CLAUDE.md`). Add new work as a
   new `DevSurfaceV<n><Increment>` class wiring engine seams; promote to a production-named class only
   after a confirmed, gated win. The engine seams (`NeighborSourceFactory`/`OcclusionScan`/… on
   `DevSurfaceV1Soa`) are the sanctioned mutation point.
5. **Document as you go.** Every decision, hypothesis, experiment, and rejected idea goes in
   `autoresearch/LOG.md` (append-only phase journal) and, when it's a durable finding, in
   `docs/performance-lessons.md` / `docs/optimization-backlog.md`. Commit frequently (`type: subject`, no
   `Co-Authored-By`). Stay on the `autoresearch` branch.
6. **Scope = tess 2 and 3.** Use `-p tess=2,3` in benchmarks. Skip tess 4. (The known float tess-4 / 16t
   collapse is out of scope — but DO check whether it begins at tess 3 / 16t, see Leads.)

## Where to read the state (read these first, every kickoff)

- `docs/performance-lessons.md` — the V1..V25 ladder, hard-won lessons, the JMH numbers, the float tess-4
  collapse diagnosis (lesson 5).
- `docs/optimization-backlog.md` — §1 A1–A7 (status), §1b the LUT-bitmask prize (open), §1c C1–C9, Tier B,
  deprioritized list, §3 the external deep-research prompt.
- `autoresearch/LOG.md` — your own running journal (create on first run if absent).
- `CLAUDE.md` — build/test/benchmark commands, freeze policy, conventions.

## Tooling (already built)

- **Benchmark:** `./bench.sh` then `python3 bench-table.py`; or the JMH jar directly for focused runs:
  `./gradlew jmhJar` then
  `java --add-modules jdk.incubator.vector -Xmx4g -jar lib/build/libs/*-jmh.jar SurfaceBench -f 3 -wi 4 -w 2s -i 8 -r 2s [-t 16] -p variantId=...,... -p tess=2,3 -p consume=AREA`.
  Add a variant by one line in `SurfaceCatalog` (test source). `consume=POINTS` exercises the zero-copy path.
- **Correctness/quality:** `./gradlew scorecard` (accuracy + duplicate/even-ness quality), `SurfaceQuality`,
  `SurfaceAccuracy`, `VariantEquivalence`. **Instrumentation:** `ScanInstrumentationTest`
  (`./gradlew scorecard --tests '*ScanInstrumentationTest'`) — count-based headroom, load-immune.
- **Profiling:** JMH `-prof gc` works (no perms). `-prof async`/`perfnorm` are BLOCKED here
  (async-profiler not installed; `perf_event_paranoid=4`). Allocation questions → `-prof gc` (rate + norm).

## Current state (update this section each phase)

- **Bit-exact default:** `DistinctPackedNumericalSurfaceV3` (= V2 + SIMD neighbour build, A6/V21).
- **Float (tess 2 ONLY):** `FloatNumericalSurfaceV2` (float build + float scan). Wins at tess 2; do NOT use
  at tess ≥ 3 with >1 thread — float scan collapses 7× (4t) / 23× (16t). **C9 CONFIRMED (Phase 2):** the
  Vector-API float (8-lane) intrinsics deopt to boxing under concurrency (alloc 13.6 MB/op → 1.084 GB/op
  at ≥2 threads); NOT fixable in source (vectors already method-local; EA isn't the lever). Double V3 is
  the tess-3 answer (it scales near-linearly).
- **Confirmed this far:** wins = A6 (SIMD build, +4–5% tess2), C1 (float build, +1.6–2.8% tess2).
  Negatives = A2 (tighter grid), A7 double (V22) & float (V25/C2), C3 region hint (V23), naive A1, A4,
  **and now the LUT-bitmask / fully-buried-atom prize (Phase 1, closed — see below).**
- **§1b LUT-bitmask prize: CLOSED-NEGATIVE (Phase 1).** The "biggest open bit-exact prize" is infeasible.
  Perfect-mask floor is only 11.9%/18% of current *scalar* tests, but the baseline is already 4-wide
  SIMD (~25% wall-clock), so the ceiling is ~2× on a ~34%-CPU op AND a cheap sound LUT is 85–92% false
  positives. **Meta-lesson: counting headroom ≠ wall-clock headroom when the baseline is vectorized.**
  Don't re-try any neighbor-major / bitmask form. (LOG Phase 1; perf-lessons "not worth doing".)
- **THE SCAN IS A DEAD WELL (3 closes: §1b bitmask, float tess-3/C9, A5 DCLM).** The direction-major SIMD
  scan (hint + early-exit) is simultaneously near-minimal in scalar tests, 4-wide vectorized, and
  sequential. Every alternative (neighbor-major bitmask, DCLM lattice, region hint) trades that for
  gather + a per-neighbor cap→direction mapping and loses in wall-clock even when it wins on test-count.
  **Decisive lens: a scalar-test-count win must beat a 4-wide-SIMD sequential baseline (~4 scalar ≈ 1 unit
  of wall-clock; gather is strictly worse than sequential).** STOP attacking the scan.
- **BUILD: re-profiled (Phase 5) and TAPPED.** Distance pass is 90% of build and already SIMD (A6); the
  candidate waste (6.27 tests/survivor, 16% survive) is inherent uniform-grid geometry and the tighter-grid
  fix is the A2 negative. Float build (C1/V24) is bit-identical on the corpus but a tess-2-only win — it
  dilutes at tess 3 where the scan dominates (Phase 6 negative).
- **ALL IN-TREE LEADS FOR TESS-2/3 SINGLE-PROTEIN THROUGHPUT ARE EXHAUSTED (6 phases, 0 promotable wins).**
  The surface is at a strong local optimum. Structural reason there's no tess-3 win left: the double scan is
  the dominant, already-optimal tess-3 cost, and float can't widen it (C9). Remaining moves are generative
  (backlog §3 deep-research prompt for novel approaches) or out-of-scope (batch/GPU C4–C5).

## Prioritized leads (tess 2 & 3 only) — re-rank each phase

- ~~**LUT-bitmask / fully-buried-atom prize.**~~ **CLOSED-NEGATIVE (Phase 1).** Infeasible; do not re-try.
- ~~**Float-path scaling at tess 3 / 16 threads.**~~ **RESOLVED (Phase 2): C9 confirmed, no source fix.**
  Float collapses at tess 3 (7×/12.5×/23× at 4/8/16t); Vector-API float intrinsics deopt to boxing under
  concurrency (1.084 GB/op). Float is tess-2-only; double V3 is the tess-3 answer. Do not re-try.
- ~~**A5 — DCLM second dot-lattice.**~~ **CLOSED (Phase 3).** Same neighbor-major spatial-narrowing family
  as §1b; its distinct survivor angle (survivors = 25% of tess-3 scan tests, DCLM could cut ~40×) is real
  but already inside Phase 1's 28%-of-current neighbor-major total, which loses to the SIMD baseline. Do
  not re-try.
- ~~**B1 — analytic per-atom SASA.**~~ **CLOSED (Phase 4).** Gate showed analytic ≠ sampled (0.3% total /
  21% per-atom at tess 2, way outside any drop-in tolerance); area-only (no points p2rank needs); likely
  slower (O(neighbors²) arc trig, no SIMD). Can't accelerate or replace the surface. Do not re-try.
- ~~**Re-profile the BUILD.**~~ **DONE (Phase 5) → TAPPED;** float-build-at-tess-3 (Phase 6) negative. See
  Current state. Distance pass 90%, already SIMD; tighter-grid is the A2 negative; float build dilutes at
  tess 3.
- ~~**Source NOVEL approaches (deep-research).**~~ **DONE (Phase 7) — no new viable lever.** Deep-research
  (5 angles, 19 sources, 24/25 claims verified) independently confirmed the search is complete: power-diagram
  exact SASA is area-only (no point set, re-derives B1); no fast-impl survey accelerator is missing; Morton
  reorder is closed by convergent evidence (V18 already sequentializes candidate reads; backlog line 60 calls
  full relabeling "modest"; research says gains "likely small, dilute at tess 3"). Report:
  `autoresearch/results/phase7-deep-research.json`.

**SEARCH CONCLUDED (7 phases, 0 promotable wins).** The tess-2/3 single-protein surface is at a strong local
optimum — see the SESSION CONCLUSION in `autoresearch/LOG.md`. Recommended surfaces stand: bit-exact default
`DistinctPackedNumericalSurfaceV3`; tess-2 float `FloatNumericalSurfaceV2`. The only remaining directions are
OUT OF THE TESS-2/3 SINGLE-PROTEIN SCOPE and need a human scope decision, not another kill-experiment:
- **C4 — whole-dataset batch / GPU throughput** (a distinct offering over thousands of proteins).
- **C5 — intra-protein parallelism** for single-large-protein latency.
- **Native power-diagram path** for an AREA-ONLY consumer (not p2rank, which needs the point set).
- A **different operating point / requirements change** (a new tessellation, an area-only mode) would reopen
  specific leads (e.g. C1 float build already wins at tess 2; analytic/power-diagram serve area-only).

**To resume:** re-arm `/loop execute the autoresearch kickoff` once one of the above is in scope. Do NOT
re-try the closed leads (scan bitmask/DCLM/region-hint, float scan at tess ≥ 3, analytic/power-diagram for
points, tighter grid, Morton relabeling) — read the SESSION CONCLUSION and `scan-is-a-dead-well` memory first.

## Per-phase workflow

1. **Pick** the top lead. Write a phase entry in `LOG.md`: hypothesis, expected payoff, kill-experiment.
2. **De-risk cheaply first** (instrumentation / a count / a `-prof gc` probe / a tiny prototype) before
   building a full variant — the A2/A7/C2 negatives all came from skipping this.
3. **Scaffold** as a new `DevSurfaceV<n>` (frozen classes untouched), wire into `SurfaceCatalog`.
4. **Gate** correctness (full corpus × config) — bit-exact or tolerance. Fail ⇒ stop, record, move on.
5. **Benchmark** idle, tess 2 & 3, single + 16t, ≥3 forks. Discard contaminated cells; re-measure idle.
6. **Verdict.** Win ⇒ record in perf-lessons (new rung row + paragraph) and consider promotion to a
   production surface. Negative ⇒ record as a documented negative (rung + lesson), keep the class.
7. **Commit** (code + docs + LOG). 
8. **Re-evaluate & rewrite THIS file** (`KICKOFF.md`): update Current state, re-rank Leads, prune dead
   ends, record the meta-lesson. Commit the updated kickoff. Then go to step 1.

## Guardrails

- Stay on `autoresearch`. Commit locally; do **not** push, publish, release, or open PRs (leave that to the
  human). Do not touch p2rank or anything outside this repo.
- If three consecutive leads dead-end, write a "state of the search" summary in LOG.md and keep going with
  fresh ideas (or the deep-research prompt in backlog §3) rather than stalling.
- Stop conditions: a clean, gated, promotable win is found and shipped (record + continue to the next
  lead), or you've exhausted the leads and generated+tested several new ones with all results documented.

---

*Phase 1 (LUT-bitmask prize): CLOSED-NEGATIVE — cheap count kill-experiment killed it (floor unreachable,
sound LUT 85–92% false positives). Meta-lesson: counting headroom ≠ wall-clock headroom vs a vectorized
baseline.*

*Phase 2 (float tess-3 scaling): RESOLVED — C9 confirmed and extended to tess 3. Benchmark in the idle
window showed float collapses 7×/12.5×/23× at 4/8/16t with alloc jumping 13.6 MB/op → 1.084 GB/op; an
EA-off probe proved it's Vector-API intrinsic deopt (not generic EA), so there is NO source fix. Float =
tess-2-only; double V3 owns tess 3. Updated lesson 5, backlog C9, CLAUDE.md.*

*Phase 3 (A5 DCLM): CLOSED — the survivor angle (25% of tess-3 scan tests, ~40× cut possible) is real but
already inside Phase 1's neighbor-major total, which loses to the SIMD baseline. Wrote a state-of-the-search
summary in LOG: THE SCAN IS A DEAD WELL — stop attacking it.*

*Phase 4 (B1 analytic SASA): CLOSED — gate showed analytic ≠ sampled (0.3%/21% gap at tess 2), area-only,
likely slower. 4th consecutive close. The scan AND both alt-output ideas are now exhausted.*

*Phase 5 (re-profile build): build is TAPPED — distance pass 90%, already SIMD; candidate waste is inherent
grid geometry; tighter-grid is the A2 negative.*

*Phase 6 (float build C1/V24 at tess 3): NEGATIVE — bit-identical on corpus but only ~1.7% faster at 1t and
neutral/noisy at 16t; build-only wins dilute at tess 3 where the double scan dominates. V24 stays tess-2.*

*State after 6 phases: 0 promotable wins; the tess-2/3 single-protein surface is at a strong local optimum.
Scan is a dead well, build is tapped, float can't widen the tess-3 scan (C9), analytic is a different/slower
output (B1). All cheap in-tree leads are exhausted.*

*Phase 7 (deep-research, generative): no new viable lever — independently confirmed the surface is at a
strong local optimum. SEARCH CONCLUDED after 7 phases / 0 promotable wins; the autoresearch loop was stopped
here (kickoff stop condition met). The remaining directions are out of the tess-2/3 single-protein scope
(C4 batch/GPU, C5 latency parallelism, native power-diagram for area-only) and need a human scope decision.*

*To resume: re-arm `/loop execute the autoresearch kickoff` only once one of those is in scope, or the
operating point / output contract changes. Read the SESSION CONCLUSION in LOG.md and the `scan-is-a-dead-well`
memory FIRST, and do not re-try any closed lead.*
