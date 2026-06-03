
# Introduction

FasterMolecularSurface is a Java library with optimized implementation of NumericalSurface from CDK.

The public surface API is captured by the `MolecularSurface` interface.

## Recommended implementation

**`DistinctPackedNumericalSurfaceV2` is the current best implementation — use it by default.**

It runs the fastest compute pipeline (cell-sorted pruned neighbour build, process-cached tessellation
and van der Waals radii, copy-free CSR neighbour access) and emits **one point per distinct surviving
direction** instead of the ~5.7× exact-coincident duplicates the icosahedral tessellation otherwise
produces — so downstream consumers need no sparsification step. The per-atom and total **surface areas
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

## Other implementations

- **`FasterNumericalSurface`** — the bit-exact reference: reproduces CDK's `NumericalSurface` point set
  and areas exactly. Use it when byte-for-byte CDK equivalence (including the full point multiplicity)
  is required; it is also the oracle the test suite checks every variant against.
- **`PackedNumericalSurface`** — bit-exact full-multiplicity production surface with a zero-copy
  `surfacePointsXYZ()` delivery path; use when CDK-identical output *and* raw-coordinate bulk access are
  both wanted.
- **`FloatNumericalSurface`** — single-precision-verdict variant of the recommended surface; ~1.05–1.14×
  faster on GraalVM (neutral on HotSpot) at a small, tolerance-bounded accuracy cost. Not bit-exact.

The full optimization history (the `DevSurfaceV1..V19` ladder) and the measured rationale behind these
choices are documented in [`docs/performance-lessons.md`](docs/performance-lessons.md).

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
  reference `NumericalSurface` within `1e-6` (3CI3 excluded — CDK throws on its cobalt atom, the
  gap the VdW fallback addresses).
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
./gradlew test
```

# Benchmarking

A coarse, opt-in timing harness (`SurfaceBenchmark`, excluded from the normal test run) compares a
variant against CDK's reference over the same corpus:
```
./gradlew benchmark
```

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
