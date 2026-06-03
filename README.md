
# Introduction

FasterMolecularSurface is a Java library with optimized implementation of NumericalSurface from CDK.

The public surface API is captured by the `MolecularSurface` interface.

## Recommended implementation

**`DistinctPackedNumericalSurfaceV2` is the current best implementation — use it by default.**

It runs the fastest compute pipeline (cell-sorted pruned neighbour build, process-cached tessellation
and van der Waals radii, copy-free CSR neighbour access) and emits **one point per distinct surviving
direction** instead of the ~5.7–6× exact-coincident duplicates the icosahedral tessellation otherwise
produces (the factor rises with tessellation level) — so downstream consumers need no sparsification step. The per-atom and total **surface areas
are bit-for-bit identical to `FasterNumericalSurface`** (the dropped duplicates are multiplicity-weighted
into the area). The occlusion scan is SIMD-vectorised (256-bit) when the JVM provides the Vector API,
with a scalar fallback otherwise.

```java
MolecularSurface surface = new DistinctPackedNumericalSurfaceV2(atomContainer, 1.4 /* solventRadius */, 4 /* tessLevel */);
double area = surface.getTotalSurfaceArea();
```

Note: because it omits the coincident duplicates, it is deliberately **not point-set-identical to CDK's
`NumericalSurface`** (the surviving point set equals what sparsification at the exact-coincidence
distance produces; the areas are exact).

### Zero-copy bulk coordinate access

A consumer that reads raw coordinates in bulk (e.g. p2rank, which turns each point into its own
primitive triple) should use the `PackedSurfaceAccess` capability instead of `getAllSurfacePoints()`:
it exposes the points as one flat `double[]` with no per-point `Point3d` allocation. The
recommended/packed surfaces implement it; detect support with `instanceof`:

```java
if (surface instanceof PackedSurfaceAccess packed) {
    double[] xyz = packed.surfacePointsXYZ();   // x0,y0,z0, x1,y1,z1, ...  (buffer may be over-allocated)
    int n = packed.surfacePointCount();         // valid data is exactly xyz[0 .. 3*n)
    for (int i = 0; i < n; i++) {
        double x = xyz[3*i], y = xyz[3*i + 1], z = xyz[3*i + 2];
        // ...
    }
}
```
The array is returned by reference for zero copy — treat it as read-only and do not retain it beyond
the surface's lifetime.

### Runtime note: the SIMD path needs an incubator module

The occlusion scan is SIMD-vectorised through the incubator Vector API, and the published artifact
cannot carry the required JVM flag. To get the vectorised path — and the speedups recorded in
[`docs/performance-lessons.md`](docs/performance-lessons.md) — the **consuming** application must launch
its JVM with `--add-modules jdk.incubator.vector`. Without it the surfaces still produce identical
results via the scalar fallback, only slower.

## Other implementations

- **`FasterNumericalSurface`** — the bit-exact reference: reproduces CDK's `NumericalSurface` point set
  and areas exactly. Use it when byte-for-byte CDK equivalence (including the full point multiplicity)
  is required; it is also the oracle the test suite checks every variant against.
- **`PackedNumericalSurface`** — bit-exact full-multiplicity production surface with a zero-copy
  `surfacePointsXYZ()` delivery path; use when CDK-identical output *and* raw-coordinate bulk access are
  both wanted.
- **`FloatNumericalSurface`** — single-precision-verdict variant of the recommended surface; ~1.05–1.14×
  faster on GraalVM (neutral on HotSpot) at a small, tolerance-bounded accuracy cost. Not bit-exact.
- **`DistinctFasterNumericalSurface`** — the `FasterNumericalSurface`-based counterpart of the distinct
  surfaces: the same CDK-exact pipeline as `FasterNumericalSurface`, but emitting one point per distinct
  direction (areas bit-for-bit identical, point set deduplicated). Primarily a reference/cross-check for
  the faster `DistinctPacked*` engine path.

The full optimization history (the `DevSurfaceV1..V19` ladder) and the measured rationale behind these
choices are documented in [`docs/performance-lessons.md`](docs/performance-lessons.md). A proposed
auto-selection factory and density-sampling API — **planned, not yet implemented** — is sketched in
[`docs/surface-api-evolution-plan.md`](docs/surface-api-evolution-plan.md). Untried optimization ideas,
the measurement work that should precede them, and an external deep-research prompt are collected in
[`docs/optimization-backlog.md`](docs/optimization-backlog.md).

# Testing

The test suite is built to be reused by future, differently-optimized variants of the surface
algorithm. It runs over a corpus of ten representative PDB structures (327–4779 atoms) defined in
`TestStructures`.

- **`AbstractMolecularSurfaceContractTest`** — implementation-independent behavioural contract
  (internal consistency, exposed-atom map, per-sphere geometry, translation invariance,
  determinism, and a pinned golden baseline). A new variant inherits the whole battery by
  subclassing it and returning the variant from the two factory methods (see
  `FasterNumericalSurfaceContractTest`).
- **`CdkEquivalenceTest`** — the oracle: proves `FasterNumericalSurface` reproduces CDK's
  reference `NumericalSurface` within `1e-6` (3CI3 excluded — CDK throws an NPE on its cobalt atom,
  whose van der Waals radius is null in CDK; the gap the VdW fallback addresses).
- **`GoldenValues` / `surface-golden.csv`** — pinned per-structure baseline. Regenerate after an
  intentional algorithm change:
  ```
  ./gradlew test --tests '*GoldenValuesGenerator' -Dgolden.regenerate=true
  ```
- **`NeighborListTest`** — validates the spatial index against a brute-force reference.
- **`ApiContractTest`** — index bounds, constructor parameters, collection immutability, input
  validation.

Run the suite:
```
./gradlew test          # default scope (see the caveat below), or ./unit-test.sh
./unit-test-all.sh      # full suite: the DevSurfaceV* ladder + cross-variant equivalence harness
```

**The default `./gradlew test` is a reduced scope.** To keep everyday runs fast it *excludes* the
per-rung `DevSurfaceV*ContractTest`s and the cross-variant equivalence harness (`GridSoaEquivalenceTest`,
`SoaEquivalenceTest`) — the bulk of the runtime. The production surfaces and `FasterNumericalSurface`
keep full contract coverage, and the distinct/float surfaces keep their equivalence-vs-Faster checks, so
the default run is a real gate; but a green `./gradlew test` does **not** mean the whole ladder was
re-verified. Run the full suite (via `./unit-test-all.sh`, i.e. `-PallTests -PtestForks=auto`) before
claiming a rung result or after touching the shared engine or any surface. Useful flags:

- `-PallTests` — include the ladder contract tests + equivalence harness.
- `-PtestForks=N|auto` — parallel test forks (`auto` = CPU cores / 2); default 1.
- `-PnoVector` — drop `jdk.incubator.vector` to exercise the scalar fallback scan.
- `-PjavaToolchain=<ver>` — build/test on a specific JDK (default 17; CI overrides per matrix entry).

# Benchmarking

The primary harness is **JMH** (`src/jmh/java/SurfaceBench`) — forked, warmed up, with run-to-run
confidence intervals, so sub-1.1× differences can be told apart from turbo jitter:

```
./bench.sh                 # pins CPU governor / disables turbo (needs root), stamps env, runs JMH -> CSV
python3 bench-table.py     # renders the speedup ladder (×CDK) from the CSV
```

`SurfaceBench` is one parameterized benchmark over `variantId` × `tess` × `consume` (`AREA` vs `POINTS`,
the latter draining the zero-copy `surfacePointsXYZ()` path the way p2rank does); thread scaling is JMH's
`-t` flag. Every variant comes from the shared `SurfaceCatalog` registry (one entry per surface), so the
default run covers the production surfaces + champion and the full `DevSurfaceV1..V19` ladder is available
via `-p variantId=V1,V9,V18,…`. Profiling is JMH's built-in `profilers` (async-profiler / `perfnorm`),
commented in the `jmh { }` block of `lib/build.gradle`.

**Accuracy + intrinsic quality scorecard** (`./gradlew scorecard`) reports, per registry variant and by
**fidelity tier** (`REFERENCE` / `BIT_EXACT` / `TOLERANCE` / `SAMPLING`): area agreement vs the exact
oracle, plus oracle-free quality metrics — **duplicate-point ratio, too-close ratio, point evenness
(NN-distance CV, min/mean)**. These need no reference, so they also judge non-bit-exact families: the
planned density samplers (see [`docs/surface-api-evolution-plan.md`](docs/surface-api-evolution-plan.md))
are scored on area convergence + quality, never required to be bit-exact. `SurfaceQualityTest` asserts the
family-relative invariants in the default suite.

**Legacy harness, kept for reproducibility:** `./gradlew benchmark` still runs the original
`@Tag("benchmark")` median-of-3 classes (`SurfaceBenchmark`, `DistinctSurfaceBenchmark`,
`DistinctPackedV2Benchmark`, `DistinctPackedV2ThreadedBenchmark`) — the method that produced the
historical ladder numbers. Its two measurement regimes (per-structure median-wall vs steady-state
aggregate-throughput) are **not interchangeable**. The historical numbers, methodology, and
GraalVM-vs-HotSpot notes live in [`docs/performance-lessons.md`](docs/performance-lessons.md); untried
ideas and the measurement-gate rationale are in [`docs/optimization-backlog.md`](docs/optimization-backlog.md).

# Molecular volume (not implemented)

CDK's `NumericalSurface` carried a `volumes[]` field that was computed but never exposed through a
getter; `FasterNumericalSurface` inherited it. The computation was removed because it was dead
(no accessor), had no oracle to validate against, and its dominant term was incorrectly normalized
(missing the `vconst` factor — values were off by a factor of roughly `3·pointDensity/(4π)`).

If a molecular-volume API is wanted later, implement it as follows. The volume enclosed by the
solvent-accessible surface is obtained from the divergence theorem, `V = (1/3) ∮ (r · n) dA`,
split into a per-atom contribution. For atom `i` with expanded radius `R = vdwRadius + solventRadius`,
each accessible tessellation point represents a surface patch of area `dA = 4πR²/pointDensity` with
outward normal `n = p` (the unit tessellation direction) at position `r = R·p + atomCenter − cp`,
where `cp` is the molecule's geometric centre (the mean of all atom coordinates). This gives

```
vconst = (4/3)·π / pointDensity
dotp1  = (atomCenter − cp) · Σp           // Σp = sum of the unit tessellation directions of accessible points
V_i    = vconst·R³·nPoints + vconst·R²·dotp1
```

and the total volume is `Σ V_i`. Note **both** terms carry `vconst` (this is the bug to avoid).
Implementation steps:

1. Restore the geometric-centre pass in `init()` (mean of atom coordinates → `cp`).
2. Keep the unit tessellation direction alongside each accepted point in `collectPoints` (the
   removed code stored it as the second element of a `Point3d[2]`); accumulate `Σp` per atom.
3. Compute `V_i` per atom with the formula above; expose `getTotalVolume()` /
   `getAllSurfaceVolumes()` (and add them to the `MolecularSurface` interface so variants share them).
4. Validate with an isolated atom: `cp == atomCenter` ⇒ `dotp1 = 0`, so `V` reduces to
   `vconst·R³·nPoints ≈ (4/3)·π·R³`, the analytic sphere volume (within tessellation tolerance).
   Add the result to the golden baseline once it checks out.
