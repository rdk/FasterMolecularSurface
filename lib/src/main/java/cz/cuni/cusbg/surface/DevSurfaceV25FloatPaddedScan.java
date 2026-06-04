package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * {@link FloatNumericalSurface} with the 8-wide float occlusion scan's scalar tail eliminated by padding
 * (backlog idea C2, {@link PaddedTailFloat256WeightedDedupOcclusionScan}). Measures whether A7 — which lost
 * on the 4-wide double scan (rung 22) — pays on the 8-wide float scan, where the tail is proportionally
 * larger. Bit-for-bit identical to {@link FloatNumericalSurface} (sentinels never bury). Scalar fallback
 * (no Vector API) uses the exact scalar weighted-dedup scan, as {@link FloatNumericalSurface} does.
 */
public class DevSurfaceV25FloatPaddedScan extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try { new PaddedTailFloat256WeightedDedupOcclusionScan(4); return true; }
        catch (Throwable t) { return false; }
    }

    public static boolean isVectorized() { return VECTOR_AVAILABLE; }

    public DevSurfaceV25FloatPaddedScan(IAtomContainer atomContainer) { this(atomContainer, 1.4, 4); }

    public DevSurfaceV25FloatPaddedScan(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) ->
                        new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new PaddedTailFloat256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);
    }
}
