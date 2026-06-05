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
- **Float (tess 2):** `FloatNumericalSurfaceV2` (float build + float scan). Wins at tess 2; do not use at
  tess ≥ 4 / many threads (float scan collapse — boxing/bandwidth hypothesis, C9).
- **Confirmed this far:** wins = A6 (SIMD build, +4–5% tess2), C1 (float build, +1.6–2.8% tess2).
  Negatives = A2 (tighter grid), A7 double (V22) & float (V25/C2), C3 region hint (V23), naive A1, A4,
  **and now the LUT-bitmask / fully-buried-atom prize (Phase 1, closed — see below).**
- **§1b LUT-bitmask prize: CLOSED-NEGATIVE (Phase 1).** The "biggest open bit-exact prize" is infeasible.
  Perfect-mask floor is only 11.9%/18% of current *scalar* tests, but the baseline is already 4-wide
  SIMD (~25% wall-clock), so the ceiling is ~2× on a ~34%-CPU op AND a cheap sound LUT is 85–92% false
  positives. **Meta-lesson: counting headroom ≠ wall-clock headroom when the baseline is vectorized.**
  Don't re-try any neighbor-major / bitmask form. (LOG Phase 1; perf-lessons "not worth doing".)
- **Open headroom now:** float-path tess-3/16t scaling (Lead #1 below), A5 DCLM dot-lattice (tess 3),
  B1 analytic SASA (tolerance). The bit-exact incremental scan well looks tapped — remaining bit-exact
  upside is most likely in the *build* or in a genuinely different algorithm (A5), not the scan.

## Prioritized leads (tess 2 & 3 only) — re-rank each phase

0. ~~**LUT-bitmask / fully-buried-atom prize.**~~ **CLOSED-NEGATIVE (Phase 1).** Infeasible; do not
   re-try. See Current state above and LOG Phase 1.
1. **Float-path scaling at tess 3 / 16 threads.** Does the FloatVector boxing/bandwidth collapse (C9) begin
   at tess 3? Measure FloatNumericalSurfaceV2 vs V3 at tess 3, -t 1/4/8/16 (idle). If tess-3 float also
   collapses at 16t, that's in scope: confirm the boxing hypothesis (`-prof gc` alloc.rate.norm) and try
   restructuring the float scan so its `FloatVector`s stay method-local and scalar-replace (lesson 3). A
   fix would make FloatNumericalSurfaceV2 a clean tess-3 win too.
3. **A5 — DCLM second dot-lattice (backlog A5).** Bin an atom's tessellation dots so each neighbor tests
   only the dots in its cap. Pays most at tess 3 (162 directions). Bit-exact. Study GROMACS `sasa.cpp`.
4. **B1 — analytic per-atom SASA (backlog B1).** Opt-in; could win at tess 3. Cheap gate first: compute
   analytic area, compare to sampled tess-2/3 on the corpus (is the gap within p2rank's tolerance?).
5. **Generate your own.** You are encouraged to invent and test new ideas; record each in the LOG with a
   hypothesis and a cheap kill-experiment before investing.

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

*Phase 1 done: the LUT-bitmask / fully-buried-atom prize (old Lead #1) is CLOSED-NEGATIVE — a cheap
count-based kill-experiment (`BitmaskFeasibilityTest`) showed the floor is unreachable and a sound cheap
LUT is 85–92% false positives, slower than the SIMD baseline. No variant scaffolded (the cheap probe
saved that cost). Meta-lesson recorded: counting headroom ≠ wall-clock headroom vs a vectorized baseline.*

*Next agent: start at Lead #1 (float-path tess-3/16t scaling). The count/alloc-probe half (`-prof gc`
alloc.rate.norm of FloatNumericalSurfaceV2 vs V3 at tess 3, and ScanInstrumentation-style counts) is
load-immune and can start now; the timing half needs an idle box (load < 1.5 — it was ~25 during
Phase 1, so the bench is still pending). De-risk cheaply first (the A2/A7/C2/§1b pattern).*
