package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.List;

/**
 * Single source of truth for the surface variants under test, shared by BOTH the JMH benchmark
 * (src/jmh) and the test-side accuracy/quality scorecard (src/test). Adding a variant is one entry
 * here, not edits across parallel arrays in several harness classes (which is what the legacy
 * {@code SurfaceBenchmark} required).
 *
 * <p>Each variant carries a {@link Fidelity} tier telling the scorecard HOW to judge it, so the
 * framework is future-proof for non-bit-exact families (e.g. the planned density samplers from
 * {@code docs/surface-api-evolution-plan.md}): a {@code SAMPLING} variant is never required to match a
 * point-set or be bit-exact — it is judged by area agreement against a reference plus the intrinsic
 * {@link SurfaceQuality} metrics, which need no oracle at all.
 */
final class SurfaceCatalog {

    private SurfaceCatalog() {}

    /** How a variant's correctness is judged (generalizes the evolution plan's {@code Fidelity}). */
    enum Fidelity {
        /** CDK {@code NumericalSurface} — the baseline everything is timed against. */
        REFERENCE,
        /** Per-atom + total areas bit-for-bit identical to {@code FasterNumericalSurface}. */
        BIT_EXACT,
        /** Close to exact within a bounded tolerance (e.g. the single-precision float verdict). */
        TOLERANCE,
        /** Alternative unit-sphere sampling: NO point-set / bit-exact oracle. Judged by area agreement
         *  against a reference + intrinsic quality only. (No such surface exists yet; the tier is ready.) */
        SAMPLING
    }

    /** Builds a surface. {@code resolution} is the tessellation level today; a density for samplers later. */
    @FunctionalInterface
    interface Factory {
        MolecularSurface build(IAtomContainer mol, double solventRadius, int resolution);
    }

    record Variant(String id, String label, Fidelity fidelity, Factory factory) {}

    /** The production + champion subset the default JMH run and scorecard cover (full ladder via {@code -p}). */
    static final List<String> DEFAULT_IDS = List.of(
            "CDK", "FASTER", "PACKED", "DISTINCT_PACKED_V2", "DISTINCT_PACKED_V3", "FLOAT", "V18", "V19");

    /** Every registered variant. Add a row to benchmark / score a new surface (incl. future samplers). */
    static final List<Variant> ALL = List.of(
            new Variant("CDK",                "CDK NumericalSurface",            Fidelity.REFERENCE, CdkSurfaceAdapter::new),
            new Variant("FASTER",             "FasterNumericalSurface",          Fidelity.BIT_EXACT, FasterNumericalSurface::new),
            new Variant("PACKED",             "PackedNumericalSurface",          Fidelity.BIT_EXACT, PackedNumericalSurface::new),
            new Variant("DISTINCT_PACKED",    "DistinctPackedNumericalSurface",  Fidelity.BIT_EXACT, DistinctPackedNumericalSurface::new),
            new Variant("DISTINCT_PACKED_V2", "DistinctPackedNumericalSurfaceV2", Fidelity.BIT_EXACT, DistinctPackedNumericalSurfaceV2::new),
            new Variant("DISTINCT_PACKED_V3", "DistinctPackedNumericalSurfaceV3", Fidelity.BIT_EXACT, DistinctPackedNumericalSurfaceV3::new),
            new Variant("DISTINCT_FASTER",    "DistinctFasterNumericalSurface",  Fidelity.BIT_EXACT, DistinctFasterNumericalSurface::new),
            new Variant("FLOAT",              "FloatNumericalSurface",           Fidelity.TOLERANCE, FloatNumericalSurface::new),
            // The DevSurfaceV1..V19 optimization ladder (available on demand: -p variantId=V1,V18,...).
            new Variant("V1",  "DevSurfaceV1Soa",          Fidelity.BIT_EXACT, DevSurfaceV1Soa::new),
            new Variant("V2",  "DevSurfaceV2Grid",         Fidelity.BIT_EXACT, DevSurfaceV2Grid::new),
            new Variant("V3",  "DevSurfaceV3Sorted",       Fidelity.BIT_EXACT, DevSurfaceV3Sorted::new),
            new Variant("V4",  "DevSurfaceV4Hinted",       Fidelity.BIT_EXACT, DevSurfaceV4Hinted::new),
            new Variant("V5",  "DevSurfaceV5Symmetric",    Fidelity.BIT_EXACT, DevSurfaceV5Symmetric::new),
            new Variant("V6",  "DevSurfaceV6LowAlloc",     Fidelity.BIT_EXACT, DevSurfaceV6LowAlloc::new),
            new Variant("V7",  "DevSurfaceV7Simd",         Fidelity.BIT_EXACT, DevSurfaceV7Simd::new),
            new Variant("V8",  "DevSurfaceV8Pruned",       Fidelity.BIT_EXACT, DevSurfaceV8Pruned::new),
            new Variant("V9",  "DevSurfaceV9Dedup",        Fidelity.BIT_EXACT, DevSurfaceV9Dedup::new),
            new Variant("V10", "DevSurfaceV10CachedMap",   Fidelity.BIT_EXACT, DevSurfaceV10CachedMap::new),
            new Variant("V11", "DevSurfaceV11CachedTess",  Fidelity.BIT_EXACT, DevSurfaceV11CachedTess::new),
            new Variant("V12", "DevSurfaceV12Flat",        Fidelity.BIT_EXACT, DevSurfaceV12Flat::new),
            new Variant("V13", "DevSurfaceV13Arena",       Fidelity.BIT_EXACT, DevSurfaceV13Arena::new),
            new Variant("V14", "DevSurfaceV14ArenaFlat",   Fidelity.BIT_EXACT, DevSurfaceV14ArenaFlat::new),
            new Variant("V15", "DevSurfaceV15LeanNbr",     Fidelity.BIT_EXACT, DevSurfaceV15LeanNbr::new),
            new Variant("V16", "DevSurfaceV16DirectNbr",   Fidelity.BIT_EXACT, DevSurfaceV16DirectNbr::new),
            new Variant("V17", "DevSurfaceV17PackedNbr",   Fidelity.BIT_EXACT, DevSurfaceV17PackedNbr::new),
            new Variant("V18", "DevSurfaceV18SortedCoords", Fidelity.BIT_EXACT, DevSurfaceV18SortedCoords::new),
            new Variant("V19", "DevSurfaceV19FlatStore",   Fidelity.BIT_EXACT, DevSurfaceV19FlatStore::new),
            new Variant("V20", "DevSurfaceV20TightGrid",   Fidelity.BIT_EXACT, DevSurfaceV20TightGrid::new),
            new Variant("V21", "DevSurfaceV21SimdBuild",   Fidelity.BIT_EXACT, DevSurfaceV21SimdBuild::new),
            new Variant("V22", "DevSurfaceV22PaddedTail",  Fidelity.BIT_EXACT, DevSurfaceV22PaddedTail::new),
            new Variant("V23", "DevSurfaceV23RegionHint",  Fidelity.BIT_EXACT, DevSurfaceV23RegionHint::new),
            new Variant("V24", "DevSurfaceV24FloatBuild",  Fidelity.TOLERANCE, DevSurfaceV24FloatBuild::new),
            new Variant("V25", "DevSurfaceV25FloatPaddedScan", Fidelity.TOLERANCE, DevSurfaceV25FloatPaddedScan::new)
    );

    static Variant byId(String id) {
        for (Variant v : ALL) if (v.id().equals(id)) return v;
        throw new IllegalArgumentException("Unknown surface variant id: " + id);
    }
}
