package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * V3 of {@link DistinctPackedNumericalSurface}: identical output to {@link DistinctPackedNumericalSurfaceV2}
 * (the same distinct-point surface, areas bit-for-bit identical to {@link FasterNumericalSurface}), with one
 * change — the neighbor build's distance pass is SIMD-vectorized (`{@link SimdDistanceCellGridNeighborList}`,
 * backlog idea A6, proven as the experimental rung {@link DevSurfaceV21SimdBuild}).
 *
 * <p>The neighbor build is ~47% of CPU at 16 threads / tess 2 (p2rank's operating point), and its inner
 * {@code d2 < (Ri+Rj)²} test was scalar; testing a 256-bit lane of candidates at once makes it **~4–5%
 * faster than V2 at tess 2** (single-thread and 16 threads; neutral at higher tess where the build is a
 * smaller share), measured by JMH on an idle, governor-pinned box. The grid, ±1 stencil, and per-pair prune
 * are unchanged, so the neighbor set — and hence the entire surface — is **bit-for-bit identical to V2**,
 * verified over the full corpus × (solvent, tess) matrix. On a JVM without {@code jdk.incubator.vector} the
 * SIMD build cannot load and the factory falls back to V2's scalar build, so the result is unchanged.
 *
 * <p>This is the recommended production surface; prefer it over V2. {@link #isVectorized()} reports the
 * active path. Like V2 it is deliberately NOT point-set-identical to CDK {@code NumericalSurface} (it omits
 * the coincident duplicates; areas are exact).
 */
public class DistinctPackedNumericalSurfaceV3 extends DevSurfaceV1Soa {

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

    public DistinctPackedNumericalSurfaceV3(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DistinctPackedNumericalSurfaceV3(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                // SIMD distance pass when the Vector API is available; else V2's scalar build.
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
