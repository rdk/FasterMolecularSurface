# CLAUDE.md

Guidance for Claude Code when working in this repository. This is a single-person scientific
project: optimize for correct, reproducible R&D, not team/process ceremony.

## What this is

A Java library with optimized implementations of CDK's `NumericalSurface`, behind the
`MolecularSurface` interface. Coordinates: `group=cz.cuni.cusbg`, `artifactId=faster-molecular-surface`,
`version=1.5`. Source in `lib/src/main/java/cz/cuni/cusbg/surface/`, tests alongside in
`lib/src/test/java/cz/cuni/cusbg/surface/`. Baseline JVM is **Java 17**; the deployment/benchmark JVM is
GraalVM 25. Recommended default surface is `DistinctPackedNumericalSurfaceV3` (see `README.md`); the full
optimization history is `docs/performance-lessons.md`.

## Build / test / benchmark

- `./gradlew build` — assemble + check (compile + the default test scope).
- `./gradlew test` — **default** test scope. It deliberately **excludes** the `DevSurfaceV*` per-rung
  contract tests and the cross-variant equivalence harness (`GridSoaEquivalenceTest`, `SoaEquivalenceTest`)
  — the bulk of the runtime. So a green `./gradlew test` does **not** mean the full ladder was verified.
- `./unit-test.sh` — the curated default run (single fork). `./unit-test-all.sh` — the **full** suite
  (`-PallTests -PtestForks=auto`): runs the DevSurface ladder + equivalence harness. Run the full suite
  **locally before committing** a shared-engine / rung / surface change — it's the fast feedback loop for
  the equivalence guarantee. CI also runs it (the `full-suite` job), so it is the backstop too, but don't
  rely on CI to first tell you an engine change broke bit-exactness.
- Flags: `-PallTests` (include the ladder + equivalence tests), `-PnoVector` (drop
  `jdk.incubator.vector` to exercise the scalar fallback), `-PtestForks=N|auto`,
  `-PjavaToolchain=<ver>` (toolchain defaults to 17; CI overrides per matrix entry).
- `./bench.sh` then `python3 bench-table.py` — the **primary** benchmark: JMH (`src/jmh/java/SurfaceBench`,
  forked + CI) over the shared `SurfaceCatalog` registry; `bench.sh` pins the machine + stamps env, the
  python script renders the ×CDK ladder from the CSV. `./gradlew benchmark` is the **legacy** median-of-3
  harness, kept only to reproduce the historical numbers (don't trust its sub-1.1× deltas).
- `./gradlew scorecard` — opt-in accuracy + oracle-free quality report (duplicate ratio, too-close ratio,
  point evenness) per variant, dispatched by fidelity tier. Future non-bit-exact sampling surfaces plug in
  here without needing a point-set oracle. `SurfaceQualityTest` asserts the invariants in the default suite.
- Regenerate the golden baseline (only after an **intentional** algorithm change):
  `./gradlew test --tests '*GoldenValuesGenerator' -Dgolden.regenerate=true`.

**Vector API quirk (important).** Compilation *always* adds `--add-modules jdk.incubator.vector`; the
test/run JVM adds it unless `-PnoVector`. When the module is absent the surfaces' `probeVector()` returns
false and they take the scalar scan path — this is **expected**, not a regression. Do not "fix" a
`-PnoVector` run that reports the scalar path.

**CI / platform note.** The matrix spans Java 17/21/25/26, GraalVM, macOS, Windows, a scalar-fallback
job, and a `full-suite` (`-PallTests`) job that runs the ladder + equivalence harness. The float-verdict path (`FloatNumericalSurface`) is intentionally JIT/platform-dependent and **not**
bit-exact — small per-platform float differences there are by design, not bugs. The CDK-exact /
area-exact surfaces are deterministic across the matrix.

## Surface implementations are frozen artifacts

Each concrete surface implementation is a **benchmarked, documented data point** whose measured
performance/accuracy is recorded in `docs/performance-lessons.md` and the golden baseline. Keeping them
side by side preserves the comparison; rewriting one invalidates its recorded benchmark.

**The frozen set (do not edit their behavior or shape):**
- The whole `DevSurfaceV1..V19` ladder (`DevSurfaceV*.java`), including the V15 negative-result branch.
- The named production surfaces: `FasterNumericalSurface`, `PackedNumericalSurface`,
  `DistinctPackedNumericalSurface`, `DistinctPackedNumericalSurfaceV2`, `DistinctPackedNumericalSurfaceV3`
  (the recommended default), `DistinctFasterNumericalSurface`, `FloatNumericalSurface`,
  `FloatNumericalSurfaceV2` (recommended approximate/float variant **at tess 2**: float build + float scan;
  the float scan collapses with many threads starting at **tess ≥ 3** (7× at 4t, 23× at 16t — Vector-API
  float intrinsics deopt to boxing under concurrency, confirmed autoresearch Phase 2) — use
  double-precision V3 at tess ≥ 3).
- The concrete strategy implementations they wire together once benchmarked (see the closure rule below).

**Classification rule — use this, NOT "is it in the perf doc".** A concrete class is frozen iff it is
constructed (`new`) or method-referenced (`::`) from a frozen surface's constructor — i.e. the transitive
closure out of the `DevSurfaceV*` rungs and the named surfaces above. Concretely that covers every
concrete `*OcclusionScan` (e.g. `Vectorized256WeightedDedupOcclusionScan`, `WeightedDedupOcclusionScan`,
`DedupVectorized256OcclusionScan`, `GlobalDedup*`, `Float256WeightedDedupOcclusionScan`), every
`*SurfacePointStore`/`*Store` (`FlatSurfacePointStore`, `DistinctFlatSurfacePointStore[V2]`,
`ListSurfacePointStore`), every `*CellGridNeighborList` / `CellGrid` / `DirectNeighborSource`, and the
shared helpers (`VdwRadiusCache`, `DirectionMapping`, `EngineScratch`, `Tessellation`). **Do not** rely on
"has a row in `docs/performance-lessons.md`" — most of these strategy classes have no doc row yet are wired
into the frozen default surface, so editing them silently invalidates its benchmark.

**What is NOT frozen — the seam *interfaces*** are the extension points: `MolecularSurface`, `OcclusionScan`,
`NeighborSource` / `NeighborSourceFactory`, `NeighborOrdering`, `SurfacePointStore` /
`SurfacePointStoreFactory`, `TessellationProvider`, `PackedSurfaceAccess`, `SurfacePointSink`. New work adds
a *new* concrete implementation of a seam (plus a new surface that wires it) — never an in-place edit to an
existing frozen concrete class.

**The one nuance:** `DevSurfaceV1Soa` is *both* a frozen rung *and* the shared engine. Its **strategy
seams** (the `NeighborSourceFactory` / `NeighborOrdering` / `OcclusionScan` / `TessellationProvider` /
`SurfacePointStoreFactory` interfaces and the `useArena`/`vdwRadius`/`directNeighbors` flags it accepts)
are the **sanctioned mutation point** — adding a *new* seam implementation or a new constructor wiring is
how the ladder grows. What is frozen is each rung's *behavior*: do not change how an existing rung
composes those seams, and keep engine changes additive (new flag/strategy defaulting to the old path) so
existing rungs stay byte-identical.

Therefore:

- **Treat the frozen set as immutable.** No clean-ups, renames, dedup-against-another-variant, or
  micro-optimizations. Add new work as a **new** `DevSurfaceVN`/named surface that inherits the shared
  contract test — never by editing the existing best in place. This is how a new default supersedes the
  old one (add the implementation, re-point the docs); the current default `DistinctPackedNumericalSurfaceV3`
  is superseded the same way, not rewritten.
- **Code reviews of frozen surfaces:** suppress *style / structure / micro-optimization* suggestions —
  they're out of scope. **Always still surface** correctness, security, resource-leak, or
  benchmark-invalidating defects — but propose the *fix as a new variant*, not an in-place edit to the
  frozen class. Reviews of the editable engine seams and any non-frozen code remain fully in scope.

## Adding a new variant

Subclass `AbstractMolecularSurfaceContractTest` (return the variant from its two factory methods) to
inherit the full behavioural battery, validate against `FasterNumericalSurface` (the bit-exact CDK
oracle) and the pinned golden baseline, then run `./unit-test-all.sh`. Record the measured speedup in
`docs/performance-lessons.md` (a new ladder row + paragraph). See `README.md` (Testing) for the
contract-test details and `docs/performance-lessons.md` for the measurement methodology.

**Before starting an optimization round**, read `docs/optimization-backlog.md`: it lists the untried
ideas (with bit-exactness and effort notes) and — importantly — the **measurement gate** that must be
rebuilt first. The last few rungs (V16–V18, all <1.1×) are inside the current benchmark's noise floor,
so adopting JMH (`@Fork`), variance reporting, profiling, and machine pinning is a prerequisite for
trusting the next round, not an afterthought.

## Conventions

- Commit messages follow the existing `type: subject` style in the history (`docs:`, `ci:`,
  `ci+test:`, or `Add <Class>: <result>` for a new surface). Do **not** add `Co-Authored-By` trailers.
- Do not run `./gradlew publish` or create releases without being asked — those are manual (and hard-blocked
  in `.claude/settings.json`'s `deny` list). `git push` is permitted (no longer denied), but still push only
  when asked or when the change is clearly ready to share — don't push speculative/WIP work unprompted.
