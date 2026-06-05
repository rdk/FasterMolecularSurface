# Surface API Evolution Plan

**Direction:** evolve this from "a faster CDK `NumericalSurface`" into a **standalone, high-performance
molecular-surface library**. Bit-exact CDK reproduction was the bootstrap and stays as one fidelity mode
(validation oracle + compat), but it is no longer the destination. The forward-looking design covers:

1. A future-proof factory that auto-selects a surface implementation (express *what*, not *how*).
2. Generalizing the engine so a new family of **custom point-sampling** surfaces (density-per-Å² instead
   of tessellation) can be added without breaking the public API.
3. **Decoupling the public API from CDK / `javax.vecmath` types** — accept plain coordinates + radii and
   return primitive arrays, so the library is usable without a CDK dependency (CDK becomes an optional
   adapter, not a hard dependency of the core API).
4. **Richer surface output** beyond points + per-atom areas — e.g. per-point normals, atom-of-origin, and
   exposed/buried classification — for downstream consumers (p2rank and others).

This is a plan, not shipped code. The guiding principle: **make the public API additive-only forever** by
generalizing the tessellation-shaped seams *before* the new implementations land. Items 1–2 are detailed
below; items 3–4 are tracked here as committed aims and will get their own design sections as they are
scoped (the new factory's `SurfaceSpec`/`Access` and the per-point output seam are the natural anchors).

---

## 1. Current landscape (as of v1.5)

- `MolecularSurface` — public interface: `getAllSurfacePoints`, `getAtomSurfaceMap`,
  `getAllSurfaceAreas`, `getTotalSurfaceArea`.
- `PackedSurfaceAccess` — optional zero-copy capability (`surfacePointsXYZ`, `surfacePointCount`),
  detected via `instanceof`. `DevSurfaceV1Soa implements MolecularSurface, PackedSurfaceAccess`, so
  **every production surface that extends the engine already supports zero-copy**.
- Production implementations, all sharing the `(IAtomContainer, double solventRadius, int tessLevel)`
  constructor:
  - `FasterNumericalSurface` — bit-exact CDK reference / test oracle (no zero-copy).
  - `PackedNumericalSurface` — bit-exact full-multiplicity + zero-copy.
  - `DistinctPackedNumericalSurfaceV2` — **recommended**: distinct points (~5.7× fewer), areas
    bit-exact vs Faster, SIMD weighted-dedup scan.
  - `FloatNumericalSurface` — single-precision-verdict variant; ~1.05–1.14× on GraalVM, neutral on
    HotSpot, not bit-exact.
- Engine strategy seams (on `DevSurfaceV1Soa`): `NeighborSourceFactory`, `NeighborOrdering`,
  `OcclusionScan`, `TessellationProvider`, `SurfacePointStoreFactory`, plus `useArena`, `vdwRadius`,
  `directNeighbors` flags.
- The `DevSurfaceV1..V19` ladder are experimental baselines — **out of scope for the factory**.

---

## 2. The auto-selection factory (additive)

Express **what** the caller needs, not **how**; let the factory resolve it against runtime
capabilities.

```java
/** WHAT the caller needs. Immutable; builder with defaults = forward-compatible. */
public final class SurfaceSpec {
    final double     solventRadius;     // default 1.4
    final Resolution resolution;        // default TessLevel(4)  -- see §3
    final Fidelity   fidelity;          // default AREA_EXACT
    final Access     access;            // default POINT3D
    final boolean    allowApproximate;  // default false -- explicit gate for the float/approx tiers
    final long       seed;              // for quasi-random samplers; ignored otherwise
    // builder / wither methods (withTessLevel, withDensity, withFidelity, allowingApproximate, ...)
    public static SurfaceSpec defaults() { ... }
}

/** Output requirement, independent of speed. Semantics are family-relative (see §3, §5). */
public enum Fidelity {
    CDK_EXACT,    // byte-identical to CDK NumericalSurface (full multiplicity)
    AREA_EXACT,   // distinct points; per-atom + total areas bit-exact vs Faster (recommended default)
    APPROXIMATE   // tolerance-bounded; fastest where it actually pays
}

/** Delivery requirement (orthogonal to fidelity). */
public enum Access { POINT3D, RAW_XYZ /* wants zero-copy PackedSurfaceAccess */ }

/** Runtime capabilities, probed once and cached. Consolidates the scattered isVectorized() probes. */
final class SurfaceCapabilities {
    static final boolean VECTOR_API; // try { new ...256Scan(4); true } catch (Throwable) { false }
    static final boolean GRAAL_JIT;  // best-effort: java.vm.name / jvmci property
}

public final class MolecularSurfaces {
    public static MolecularSurface create(IAtomContainer mol) { return create(mol, SurfaceSpec.defaults()); }

    public static MolecularSurface create(IAtomContainer mol, SurfaceSpec s) {
        // branch on family (resolution) FIRST, then fidelity/access/capabilities within it.
        return Selector.resolve(s).build(mol, s);
    }

    /** Decide without building, for reproducibility logging: {implName, reason}. */
    public static Selection plan(SurfaceSpec s) { ... }
}
```

Selection sketch (tessellation family):

```java
switch (s.fidelity) {
    case CDK_EXACT  -> s.access == Access.RAW_XYZ ? new PackedNumericalSurface(mol, r, lvl)
                                                  : new FasterNumericalSurface(mol, r, lvl);
    case AREA_EXACT -> new DistinctPackedNumericalSurfaceV2(mol, r, lvl);   // vectorized w/ scalar fallback
    case APPROXIMATE -> (s.allowApproximate && CAP.VECTOR_API && CAP.GRAAL_JIT && lvl >= 3)
                            ? new FloatNumericalSurface(mol, r, lvl)        // measured net win only here
                            : new DistinctPackedNumericalSurfaceV2(mol, r, lvl);
}
```

**Future-proofing:** stable `MolecularSurface` return type; new impls attach to a `Fidelity`/family
value with no signature change; defaulted builder makes new axes additive; the selection **policy**
should be overridable (`create(mol, spec, SurfacePolicy)` SPI) so heuristics can evolve or be pinned.

---

## 3. The new family: custom density-based point sampling

Planned: surfaces that sample sphere points from a **selectable custom function** at an **approximate
expected density per 1 Å²**, instead of the icosahedral tessellation.

### Why it is structurally different

| | Tessellation (today) | Custom sampler (planned) |
|---|---|---|
| Parameter | discrete **level** | continuous **density** (pts/Å²) + sampler fn |
| Points per atom | same for every atom (`numTess`) | scales with r² (≈ `density·4πr²`) → varies per atom |
| Direction set | one shared set | per-atom (radius-dependent), possibly quasi-random/seeded |
| Area | `4πr²·count / pointDensity` (per sphere) | `count / density` (per-point area = `1/density`) |
| Duplicates | ~5.7× coincident → dedup is the big win | points already distinct → no dedup, no weighting |

Two consequences:

1. **The dedup/distinct family is tessellation-only.** `DirectionMapping`, the weighted scans, and the
   distinct stores collapse the icosahedral 5.7× duplication. Sampling points are already distinct, so
   the sampling family's champion is a **plain** SIMD scan + flat store — *not* the distinct machinery.
   The factory must not assume "distinct is universally best."
2. **`OcclusionScan.collect(tx, ty, tz, numTess, …)` already takes directions per call**, so per-atom
   variable direction sets need **no scan-interface change** — the plain scans iterate whatever they
   are handed. The only tessellation assumption is that today's engine computes one direction set
   before the atom loop.

### The two seams to generalize

**(a) Per-point area weight** — replace the per-sphere normalization in `DevSurfaceV1Soa` (currently
`areas[i] = 4π·r²·count / pointDensity`) with a sum of per-point area weights over survivors. This one
generalization subsumes all three families:
- tessellation: `w = 4πr² / pointDensity`
- sampling at density d: `w = 1/d`
- weighted-dedup: `w × multiplicity`

**(b) A sphere point-set provider** generalizing `TessellationProvider`, called **per atom** (vs once):

```java
/** Supplies an atom's sphere sample. Tessellation is one impl; density samplers are others. */
public interface SpherePointSet {
    int directions(double radius, DirectionsOut out); // unit dirs into engine scratch; count may depend on r
    double pointArea(double radius);                   // tess: 4πr²/density-per-sphere; sampler: 1/densityPerÅ²
    boolean isDeterministic();                          // false for unseeded quasi-random
}
```

Tessellation wraps the existing cached directions (count independent of radius). A
`FibonacciSampler(densityPerÅ², seed)` returns `round(d·4πr²)` spiral points. Prefer deterministic
samplers (Fibonacci/spiral) as defaults to keep output reproducible.

### Resolution parameter generalization (do this NOW)

Replace the bare `int tessLevel` in `SurfaceSpec` with a family-tagged resolution so adding sampling
later is additive rather than a breaking change:

```java
public sealed interface Resolution {
    record TessLevel(int level)          implements Resolution {} // tessellation family
    record Density(double pointsPerSqA)  implements Resolution {} // sampling family
    record Sampler(SpherePointSet sampler) implements Resolution {} // fully custom
}
```

The factory branches on `resolution` first (which family), then fidelity/access/capabilities within it.

---

## 4. Refactor scope: do we need to change all APIs?

**No — not "all APIs," but more than the pure-additive factory.**

- **Additive / unchanged:** `MolecularSurface`, `SurfacePointStore`/`Sink`, the plain `OcclusionScan`
  contract (already per-call directions), `PackedSurfaceAccess`, all existing constructors, the entire
  dedup family.
- **One contained engine-seam refactor:** (a) per-point **area weights** instead of the per-sphere
  `pointDensity` formula, and (b) generalize `TessellationProvider` → `SpherePointSet`, called **per
  atom** in the engine loop. Touches `DevSurfaceV1Soa` + one new interface — shared infra, but localized.
- **Optional polish:** add `PackedSurfaceAccess` to `FasterNumericalSurface` so zero-copy is uniform;
  consolidate the three per-class `isVectorized()` probes into `SurfaceCapabilities`.

**Recommended order:** generalize the **`Resolution` spec param** and the **per-point area-weight seam**
first (insurance: locks the public API), then add `SpherePointSet` + the first sampler, then expand the
factory's selection to the sampling family.

---

## 5. Issues / risks

1. **Reproducibility across environments.** Auto-selection can return different output for the same
   `(mol, spec)`. Confined to `APPROXIMATE`/sampling tiers — `CDK_EXACT` and `AREA_EXACT` (tessellation)
   are deterministic everywhere (vectorized scan is bit-exact to scalar; only float flips boundary
   points). Mitigations: default to exact/deterministic tiers; gate float behind `allowApproximate`;
   make `seed` first-class for samplers; expose the chosen impl + seed via `plan()`; offer a "pinned"
   mode that ignores capabilities for archival runs.
2. **Heuristic rot.** "Float on GraalVM at tess≥3" is one-box, benchmark-derived (we already corrected a
   1.45×→1.12× artifact). Keep the policy conservative, overridable, re-validated against the benchmark.
3. **JIT-vendor detection is not a stable contract.** GraalVM-vs-HotSpot via `java.vm.name`/jvmci is
   best-effort; a wrong guess costs a little perf, never correctness. Allow override.
4. **Construction does the work** (eager `init()` in the constructor). The factory must decide from
   spec + static capabilities before building — fine for heuristics, rules out adaptive
   build-measure-pick without a laziness refactor. Future option, not now.
5. **Spec needs two orthogonal axes**, not one. `CDK_EXACT` alone is ambiguous (Faster vs Packed); the
   differentiator is delivery (`Access`). Keep fidelity and access separate.
6. **No CDK oracle for the sampling family.** The bit-exact equivalence harness does not apply. Needs a
   new correctness story: convergence to analytic sphere area for an isolated atom, and statistical
   agreement (within tolerance) against a high-density reference.
7. **`Fidelity` semantics diverge per family.** Sampling has no `CDK_EXACT`; "fidelity" becomes the
   target density itself. Interpret fidelity relative to the family, or split into family-specific knobs.
8. **Density→count edge cases.** Rounding, a floor for tiny radii, and whether "approximate density" is
   an expected or exact guarantee must be pinned in the `SpherePointSet` contract.
9. **Load imbalance.** Per-atom point count ∝ r², so large atoms dominate the scan. Throughput-fine;
   the neighbor build is sampling-independent.
10. **Scope the factory's universe.** Exclude the `DevSurfaceV1..V19` ladder (experimental baselines).

---

## 6. Concrete next steps (suggested)

1. Generalize the engine area math to per-point weights; verify bit-exact against all current variants
   (tessellation weight = `4πr²/pointDensity` reproduces today's numbers exactly).
2. Introduce `SpherePointSet`; reimplement `TessellationProvider` as one impl on top of it; engine calls
   it per atom. Re-run the equivalence + golden suites (must stay bit-exact).
3. Add `SurfaceSpec` with `Resolution`, `Fidelity`, `Access`, `seed`; add `MolecularSurfaces.create` +
   `plan`; consolidate capability probes. Cover the tessellation family only at first.
4. Implement the first density sampler (deterministic Fibonacci/spiral) as a new surface; add its
   statistical correctness tests; wire it into the factory under `Resolution.Density`.
5. Document the chosen selection policy and the reproducibility contract in the README + this doc.

See also [`performance-lessons.md`](performance-lessons.md) for the measured rationale behind the
existing variant choices (esp. lesson 5 on the JVM-dependent float result).
