
# Introduction

FasterMolecularSurface is a Java library with optimized implementation of NumericalSurface from CDK.

The public surface API is captured by the `MolecularSurface` interface; `FasterNumericalSurface`
is its current, optimized implementation.

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
