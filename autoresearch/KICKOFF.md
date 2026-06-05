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
- **Open headroom now: the BUILD.** The scan is dead (Phases 1–3) and both alt-output ideas are closed
  (float tess-3 = C9 confirmed Phase 2; analytic SASA = B1 closed Phase 4 — not interchangeable with the
  sampled surface, area-only, likely slower). The build (neighbor-list construction) is the remaining
  lever: A6/SIMD build was the last real win, it's the larger share at tess 2 / 16t (bandwidth-bound),
  and it has NOT been re-profiled since. That's the next target.

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
1. **Re-profile the BUILD — TOP LEAD.** The only remaining lever. A6 (SIMD build, V21) was the last win
   and the build is the larger share at tess 2 / 16t (bandwidth-bound, lesson 12). **De-risk first
   (mostly load-immune):** read the build path (`*CellGrid*NeighborList`, `SimdDistanceCellGridNeighborList`,
   `DevSurfaceV21SimdBuild` wiring) and break build time into phases — grid construction, the d²<sumR²
   distance pass, neighbor materialization/sort. Then one idle-box `-prof gc` + timing run isolating build
   vs scan (e.g. compare a tess where scan≈0 work, or instrument). Find the dominant build sub-cost BEFORE
   proposing a variant. Candidate levers once profiled: cheaper grid binning, fewer passes, better memory
   layout for the distance pass (it's bandwidth-bound, so traffic reduction is the lever A6/C1 exploited).
2. **Generate your own** if the build re-profile shows no headroom (e.g. backlog §3 deep-research prompt,
   or the deployment-shaped C4/C5 batch ideas — though those are out of the tess-2/3 single-protein scope).

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

*Next agent: start at Lead #1 = re-profile the BUILD (the only remaining lever; the scan is a dead well).
De-risk load-immune first: read the build path and decompose build time (grid / distance pass / neighbor
materialization); then ONE idle-box `-prof gc` + timing run to isolate the dominant build sub-cost before
proposing any variant. Heed lesson 12 (build is bandwidth-bound at 16t → traffic reduction is the lever,
as A6/C1 found) and the Phase-1 meta-lesson (count wins must beat the SIMD+sequential baseline). Timing
needs load < 1.5 — this shared box swings 1.0–25, so poll /proc/loadavg first; counting/reading is
load-immune. If the build also proves tapped, switch to generating fresh ideas (backlog §3).*
