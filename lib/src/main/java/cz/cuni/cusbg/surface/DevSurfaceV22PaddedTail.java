package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Stacks backlog ideas A6 + A7 on {@link DistinctPackedNumericalSurfaceV2}: the SIMD-vectorized neighbor
 * build of {@link DevSurfaceV21SimdBuild} ({@link SimdDistanceCellGridNeighborList}, A6) <em>and</em> the
 * padded-tail occlusion scan {@link PaddedTailVectorized256WeightedDedupOcclusionScan} (A7), which removes
 * the scalar remainder loop from the per-direction neighbor scan. Everything else — cached
 * tessellation/VdW, distinct flat store, copy-free CSR access — is V2's.
 *
 * <p>Output is <b>bit-for-bit identical</b> to {@link DistinctPackedNumericalSurfaceV2}: A6 yields the
 * same neighbor set and A7's sentinels never bury a direction. On a JVM without {@code jdk.incubator.vector}
 * both SIMD pieces fall back to V2's scalar build + scalar weighted-dedup scan, so the result is unchanged.
 * {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV22PaddedTail extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new PaddedTailVectorized256WeightedDedupOcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the SIMD build + padded SIMD scan are in use; false on a JVM without the Vector API. */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public DevSurfaceV22PaddedTail(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV22PaddedTail(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                VECTOR_AVAILABLE
                        ? (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new SimdDistanceCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get)
                        : (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new PaddedTailVectorized256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access
    }
}
