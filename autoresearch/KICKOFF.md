# Autoresearch kickoff — DISPATCHER (no active round)

This is the entry point for `/loop execute the autoresearch kickoff`. It is currently in the **idle**
state: **there is no active research round.**

## If you are an agent executing this with no active round

1. **Read `META.md` first** — closed leads (do not re-try them), meta-lessons, tooling gotchas.
2. **Do NOT invent a round or re-run a concluded mission.** The last round
   (`archive/2026-06-05_tess2-3-throughput/`) concluded that the tess-2/3 single-protein point-set surface
   is at a **strong local optimum** — 8 phases, 0 promotable wins, all in-tree leads + the real
   `consume=POINTS` consumer pattern exhausted, plus an external deep-research sweep. Re-running that scope
   would only re-confirm it (= manufacture churn, which the rules forbid).
3. **A new round needs a SCOPE the human chooses** (it changes the contract/operating point). The
   productive directions left, all out of the prior scope:
   - **C4** whole-dataset batch / GPU throughput (a distinct offering over thousands of proteins).
   - **C5** intra-protein parallelism for single-large-protein latency.
   - A native **power-diagram** O(N log N) path for an **area-only** consumer (not p2rank — needs points).
   - A different **operating point / output contract** (new tessellation, area-only mode).
4. **Action:** if the user has named a scope, **arm a round** (below). If not, surface the scope choice and
   **stop** (do not schedule another wake-up). Re-confirming a known optimum is not work.

## Arming a new round

1. `cp templates/KICKOFF.template.md runs/<YYYY-MM-DD>_<slug>/KICKOFF.md`; fill in Mission, Scope, Current
   state, and the ranked Leads (seed "Current state" from `META.md` §1 closed leads + the prior REPORT).
2. Start `runs/<id>/LOG.md` (append-only journal; first entry = Phase 0 setup).
3. Replace this dispatcher's "no active round" status with a one-line pointer to `runs/<id>/KICKOFF.md`,
   and from then on the loop executes **that** file (self-updating per phase). On completion, write
   `runs/<id>/REPORT.md`, fold lessons into `META.md`, `git mv runs/<id> archive/<id>`, and reset this
   dispatcher to idle.

## Operating rules (full set in `README.md` / `META.md`)

Never benchmark on a loaded box (poll `/proc/loadavg`, 1-min < 1.5). Gate correctness before any speed
claim. Never edit frozen classes — instrument via test-scope clones. De-risk cheaply before building.
Commit often (`type: subject`, no `Co-Authored-By`), stay on `autoresearch`, never push/publish. Conclude
honestly; don't manufacture churn.
