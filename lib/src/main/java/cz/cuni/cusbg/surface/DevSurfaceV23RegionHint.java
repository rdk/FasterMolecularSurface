package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * {@link DistinctPackedNumericalSurfaceV3} with the occlusion scan's single last-occluder hint replaced
 * by a per-octant hint (backlog idea C3, {@link RegionHintedVectorized256WeightedDedupOcclusionScan}).
 * Bit-for-bit identical to V3; measures whether a region-binned hint beats the single hint.
 * Scalar fallback (no Vector API) uses the exact scalar build + scalar weighted-dedup scan.
 */
public class DevSurfaceV23RegionHint extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try { new RegionHintedVectorized256WeightedDedupOcclusionScan(4); return true; }
        catch (Throwable t) { return false; }
    }

    public static boolean isVectorized() { return VECTOR_AVAILABLE; }

    public DevSurfaceV23RegionHint(IAtomContainer atomContainer) { this(atomContainer, 1.4, 4); }

    public DevSurfaceV23RegionHint(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                VECTOR_AVAILABLE
                        ? (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new SimdDistanceCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get)
                        : (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new RegionHintedVectorized256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);
    }
}
