package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * {@link DistinctPackedNumericalSurfaceV2} with its neighbor build's distance pass SIMD-vectorized
 * (backlog idea A6): the build is ~47% of CPU at 16 threads / tess 2 on V2 (p2rank's regime), and its
 * inner {@code d2 < (Ri+Rj)²} test was scalar. {@link SimdDistanceCellGridNeighborList} keeps the same
 * grid, stencil, and prune (so the neighbor set is identical) but tests a lane of candidates at once.
 * Everything else — the weighted dedup SIMD scan, cached tessellation/VdW, distinct flat store, copy-free
 * CSR access — is V2's.
 *
 * <p>Output is <b>bit-for-bit identical</b> to {@link DistinctPackedNumericalSurfaceV2}. On a JVM without
 * {@code jdk.incubator.vector} the SIMD build cannot load, so the neighbor factory falls back to the
 * scalar {@link CoordSortedPrunedSymmetricCellGridNeighborList} (still V2's build) and the result is
 * unchanged. {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV21SimdBuild extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new Vectorized256WeightedDedupOcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the vectorized scan + vectorized build are in use; false on a JVM without the Vector API. */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public DevSurfaceV21SimdBuild(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV21SimdBuild(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                // SIMD distance pass when the Vector API is available; else V2's scalar build (the SIMD
                // class is referenced only inside the branch taken when VECTOR_AVAILABLE, so it never
                // loads on a JVM lacking the incubator module).
                VECTOR_AVAILABLE
                        ? (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new SimdDistanceCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get)
                        : (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new Vectorized256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access
    }
}
