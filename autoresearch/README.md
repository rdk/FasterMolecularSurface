# autoresearch/ — the performance-research harness

This directory drives **autonomous, multi-hour performance-research rounds** on the SASA kernel, run via
`/loop execute the autoresearch kickoff` on the `autoresearch` branch. It is built to be *reused* across
rounds: each round is self-contained and archived; the durable lessons accumulate at the top level so a
new round starts from hard-won state instead of a blank page.

## Read order for a new agent

1. **`META.md`** — accumulated cross-round lessons, permanently-closed leads (do **not** re-try these),
   machine/tooling gotchas, and the reusable instrumentation catalog. **Read this first, every time.**
2. **`KICKOFF.md`** — the live dispatcher. Tells you whether a round is active and, if not, how to arm one.
3. The active round's `runs/<id>/KICKOFF.md` (if one is active) — its mission, current state, and ranked leads.

## Layout

```
autoresearch/
  README.md                  # this file — conventions
  META.md                    # durable lessons + closed leads + tooling notes + archive index (READ FIRST)
  KICKOFF.md                 # live dispatcher: active round pointer, or how to arm a new one
  templates/
    KICKOFF.template.md      # skeleton for a new round's mission/state/leads/workflow
  runs/<id>/                 # an ACTIVE round (KICKOFF.md + LOG.md + results/), promoted to archive/ when done
  archive/<date>_<slug>/     # COMPLETED rounds, frozen for the record
    KICKOFF.md  LOG.md  REPORT.md  results/
```

A round id is `<YYYY-MM-DD>_<short-slug>` (e.g. `2026-06-05_tess2-3-throughput`).

## Lifecycle of a round

1. **Arm.** Run the **`autoresearch-new-round`** skill (`/autoresearch-new-round <scope>`) — it reads
   `META.md`, confirms the scope, copies `templates/KICKOFF.template.md` to `runs/<id>/KICKOFF.md`, fills in
   the mission/state/leads (seeded with the closed leads), starts `runs/<id>/LOG.md`, points the dispatcher
   `KICKOFF.md` at the active round, and commits. (Or do those steps by hand.)
2. **Run** (per-phase loop, see the template's workflow): pick the top lead → write a LOG phase entry
   (hypothesis · cheap kill-experiment · expected payoff) → **de-risk cheaply first** (load-immune
   counting/instrumentation, a `-prof gc` probe, a tiny prototype) → scaffold a new variant only if it
   survives → **gate correctness before any speed claim** → benchmark on an **idle** box → verdict →
   commit → rewrite the round's KICKOFF (state + re-ranked leads). Repeat.
3. **Conclude.** Write a `REPORT.md` (consolidated summary). Fold durable findings into
   `docs/performance-lessons.md` / `docs/optimization-backlog.md` and into **`META.md`** (closed leads +
   any new meta-lesson). `git mv runs/<id>` → `archive/<id>`. Reset `KICKOFF.md` to the idle dispatcher.

## Hard rules (also in `META.md`, repeated because they matter)

- **Never benchmark on a loaded shared box.** Poll `/proc/loadavg`; proceed only if 1-min load < 1.5.
  Counting/instrumentation/correctness work is load-immune — do it anytime.
- **Gate correctness before any speed claim** (bit-exact vs the oracle, or the tolerance envelope, over
  the full corpus × config).
- **Never edit a frozen surface/strategy class** (see root `CLAUDE.md`). Instrument via test-scope clones;
  add new work as a new `DevSurfaceV<n>` / named surface.
- **Commit frequently**, `type: subject`, no `Co-Authored-By`. Stay on `autoresearch`; never push/publish.
- **De-risk cheaply before building.** Most negatives this project found were killable by a count, not a
  variant. Don't scaffold + benchmark what a load-immune measurement can refute.
- **Conclude honestly.** If leads are exhausted, say so and stop — don't manufacture marginal churn.
