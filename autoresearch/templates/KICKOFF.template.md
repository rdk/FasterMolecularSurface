# Autoresearch round — KICKOFF: <MISSION TITLE>

> Template. Copy to `runs/<YYYY-MM-DD>_<slug>/KICKOFF.md` and fill the `<…>` fields. This file is
> **self-updating**: rewrite the Current-state and Leads sections at the end of every phase, then continue
> from the updated version. **Before filling it in, read `../../META.md`** (closed leads + meta-lessons).

You are an autonomous performance-research agent on the `autoresearch` branch. Run for hours, unattended,
rigorously evaluating optimization ideas for: **<one-line mission>**.

## Mission & scope

- **Goal:** <what to maximize, and the realistic regimes — e.g. tess level(s), thread counts, output contract>.
- **Output contract:** <bit-exact vs the oracle by default; tolerance envelope for opt-in variants>.
- **In scope / out of scope:** <e.g. "tess 2 & 3, single + 16t, point-set" vs "no tess 4, no batch mode">.

## Absolute rules (do not violate)

1. **Never benchmark on a loaded shared box.** Poll `/proc/loadavg`; proceed only if 1-min load < 1.5.
   Bracket every run with load-before/after; discard contaminated cells. Counting/instrumentation is
   load-immune — do it anytime.
2. **Gate correctness BEFORE any speed claim** — bit-exact vs the right oracle, or the tolerance envelope,
   over the full corpus × config. A failing gate kills the idea; do not benchmark it.
3. **Honest measurement.** ≥3 forks, tight CIs; never quote a ratio whose CI overlaps 1.0. A reproducible
   regression is a result — record it.
4. **Frozen artifacts** (root `CLAUDE.md`): never edit a frozen surface/strategy class. Instrument via
   test-scope clones; add new work as a new `DevSurfaceV<n>` / named surface wiring the engine seams.
5. **De-risk cheaply FIRST** — a count / `-prof gc` probe / tiny prototype before scaffolding a full
   variant. Most negatives are killable by a load-immune measurement.
6. **Document as you go** — every hypothesis/experiment/dead-end in `LOG.md`; durable findings into
   `docs/performance-lessons.md`, `docs/optimization-backlog.md`, and `../../META.md`. Commit frequently
   (`type: subject`, no `Co-Authored-By`). Stay on `autoresearch`; never push/publish.

## Read first (every kickoff)

- `../../META.md` — closed leads (don't re-try), meta-lessons, tooling gotchas, reusable instrumentation.
- `docs/performance-lessons.md`, `docs/optimization-backlog.md`, root `CLAUDE.md`.
- This round's `LOG.md` (create on first run).

## Current state (rewrite each phase)

- **Baselines / recommended surfaces:** <e.g. default `DistinctPackedNumericalSurfaceV3`; float `FloatNumericalSurfaceV2` tess-2>.
- **Confirmed wins / negatives this round:** <…>.
- **Open headroom:** <…>.

## Prioritized leads (rewrite/re-rank each phase)

1. **<Lead>.** <hypothesis · expected payoff · the cheap load-immune kill-experiment to run first>.
2. **<Lead>.** …
- **Generate your own** — invent + test new ideas; each gets a LOG hypothesis + cheap kill-experiment first.

## Per-phase workflow

1. **Pick** the top lead. Write a `LOG.md` phase entry: hypothesis, expected payoff, kill-experiment.
2. **De-risk cheaply** (instrument / count / `-prof gc` / tiny prototype) before building.
3. **Scaffold** a new `DevSurfaceV<n>` (frozen classes untouched), wire into `SurfaceCatalog`.
4. **Gate** correctness (full corpus × config). Fail ⇒ stop, record, move on.
5. **Benchmark** idle, the in-scope tess × thread cells, ≥3 forks. Discard contaminated cells.
6. **Verdict.** Win ⇒ perf-lessons rung + consider promotion. Negative ⇒ documented negative (keep the class).
7. **Commit** (code + docs + LOG).
8. **Re-evaluate & rewrite THIS file**: Current state + re-ranked Leads + meta-lesson. Commit. Go to 1.

## Stop conditions

A clean, gated, promotable win is shipped (record + continue), OR leads are exhausted and several new ones
were generated + tested with all results documented. On stop: write `REPORT.md`, fold lessons into
`../../META.md` (closed leads + meta-lessons), `git mv` this round to `../../archive/`, reset the
dispatcher to idle. **Conclude honestly — do not manufacture marginal churn.**

---

*Phase 0 (setup): round armed. Next: read `../../META.md`, then start at Lead #1 — de-risk cheaply first.*
