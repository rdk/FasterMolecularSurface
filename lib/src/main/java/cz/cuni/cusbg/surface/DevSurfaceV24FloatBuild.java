package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * {@link DistinctPackedNumericalSurfaceV3} with the SIMD neighbor build done in single precision (backlog
 * idea C1, {@link FloatSimdDistanceCellGridNeighborList}) — half the candidate-coordinate traffic to attack
 * the bandwidth ceiling A6 hit at 16 threads. The occlusion scan and everything else stay V3's.
 *
 * <p><b>Not bit-exact:</b> the float distance test can classify a cutoff-boundary pair differently, so the
 * neighbor set (and surface) can differ slightly from V3 — tolerance-tested, not via the bit-exact harness.
 * Scalar fallback (no Vector API) uses the exact double build, so it is bit-exact there.
 */
public class DevSurfaceV24FloatBuild extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try { new Vectorized256WeightedDedupOcclusionScan(4); return true; }
        catch (Throwable t) { return false; }
    }

    public static boolean isVectorized() { return VECTOR_AVAILABLE; }

    public DevSurfaceV24FloatBuild(IAtomContainer atomContainer) { this(atomContainer, 1.4, 4); }

    public DevSurfaceV24FloatBuild(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                VECTOR_AVAILABLE
                        ? (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new FloatSimdDistanceCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get)
                        : (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new Vectorized256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);
    }
}
