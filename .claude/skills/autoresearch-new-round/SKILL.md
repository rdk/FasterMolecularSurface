---
name: autoresearch-new-round
description: >-
  Arm a new autoresearch performance-research round for this repo's SASA kernel: pick/confirm a scope,
  scaffold autoresearch/runs/<date>_<slug>/ from the template (seeded with the already-closed leads from
  META.md so they're not re-tried), point the dispatcher KICKOFF at it, start the LOG, commit, and ready
  the /loop. Use when the user wants to start, arm, or set up a new autoresearch round / a new optimization
  investigation round on the molecular-surface kernel. NOT for executing an already-armed round (that is
  `/loop execute the autoresearch kickoff`).
---

# Arm a new autoresearch round

Sets up a fresh, self-contained research round under `autoresearch/runs/<id>/` and points the dispatcher
at it, so `/loop execute the autoresearch kickoff` runs the new mission. See `autoresearch/README.md` for
the directory conventions and `autoresearch/META.md` for the accumulated lessons.

`$ARGUMENTS` (optional) is the scope/mission for the round (free text, e.g. "whole-dataset batch throughput
mode" or "tess-4 single-protein latency"). If empty or ambiguous, ask the user for it before scaffolding.

## Procedure

1. **Read the accumulated state FIRST** (non-negotiable — this is the whole point of the harness):
   - `autoresearch/META.md` — the **permanently-closed leads** (do NOT seed a round around them), the
     meta-lessons, and the machine/tooling gotchas.
   - The most recent `autoresearch/archive/*/REPORT.md` — what the last round concluded.
   - Confirm there is **no active round** already: check `autoresearch/KICKOFF.md` (the dispatcher). If it
     already points to an active `runs/<id>/`, stop and tell the user a round is in progress (offer to
     resume it via `/loop` or to conclude+archive it first).

2. **Determine the scope.** If `$ARGUMENTS` gives a clear scope, use it. Otherwise ask the user with
   `AskUserQuestion` — the productive out-of-prior-scope directions are listed in the dispatcher
   (batch/GPU throughput C4; single-protein latency C5; an area-only fast path; a new operating point).
   **Reject a scope that only re-opens a closed lead** (per META.md §1) unless the user explicitly wants to
   re-validate one — name the closure and confirm before proceeding.

3. **Pick the round id:** `<YYYY-MM-DD>_<short-kebab-slug>` (date = today from the environment's current
   date; slug = 2–4 words of the scope, e.g. `2026-07-01_batch-throughput`).

4. **Scaffold the round** (these are mostly mechanical; do them, then verify):
   - `mkdir -p autoresearch/runs/<id>`
   - `cp autoresearch/templates/KICKOFF.template.md autoresearch/runs/<id>/KICKOFF.md`
   - Fill in the template's `<…>` fields in `runs/<id>/KICKOFF.md`: **Mission & scope** (from step 2),
     **Output contract** (bit-exact vs the oracle by default; the tolerance envelope for opt-in variants),
     in/out of scope, **Current state** (seed the baselines + the relevant closed leads from META.md so the
     round starts informed), and a first **ranked Leads** list (3–5, each with a cheap load-immune
     kill-experiment). Keep the rules / per-phase workflow / stop conditions sections as-is.
   - Create `autoresearch/runs/<id>/LOG.md` with a header and a **Phase 0 — setup** entry (mission, inherited
     state, the top lead, and "next: de-risk cheaply first").
   - Create `autoresearch/runs/<id>/results/` (empty; benchmark logs land here).

5. **Point the dispatcher at the active round.** Edit `autoresearch/KICKOFF.md` so its top clearly states
   the active round and instructs the loop to **execute `runs/<id>/KICKOFF.md`** (replace the "no active
   round / idle" status with an `ACTIVE ROUND: runs/<id>/KICKOFF.md` pointer + one-line scope). Keep the
   operating-rules footer.

6. **Commit** (`type: subject`, no `Co-Authored-By`; stay on the `autoresearch` branch; do not push):
   `git add -A && git commit -m "autoresearch: arm round <id> — <scope>"`.

7. **Hand off to the loop.** Tell the user the round is armed and that running
   **`/loop execute the autoresearch kickoff`** will start it (self-paced: load-immune work anytime,
   timing gated on an idle box). Do **not** start the multi-hour loop yourself unless the user asks —
   arming and running are separate steps.

## Guardrails

- Read `autoresearch/META.md` closed-leads before seeding leads; never scaffold a round whose top leads are
  already closed.
- Mechanical file ops (mkdir/cp/edit/commit) only — **never** edit a frozen surface/strategy class (root
  `CLAUDE.md`); the round itself instruments via test-scope clones.
- One active round at a time. Conclude + `git mv runs/<id> archive/<id>` + reset the dispatcher to idle
  before arming the next (the round's own stop-condition workflow covers this).
