package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * V2 of {@link FloatNumericalSurface}: the single-precision distinct surface with <b>both</b> float
 * optimizations stacked — the float-verdict occlusion scan ({@link Float256WeightedDedupOcclusionScan},
 * as in V1) <em>and</em> the float-precision SIMD neighbour build ({@link FloatSimdDistanceCellGridNeighborList},
 * backlog idea C1, proven as the rung {@link DevSurfaceV24FloatBuild}). V1 floated only the scan; this also
 * floats the build's distance test, halving the candidate-coordinate traffic — the bandwidth lever that
 * scales with thread count.
 *
 * <p>So both halves of the kernel now run in single precision: the build's {@code d²<(Ri+Rj)²} test (C1,
 * wins where the build dominates — tess 2, many threads) and the scan's buried verdict (V1, wins where the
 * scan dominates — tess 3–4). Point POSITIONS and area weights stay {@code double}
 * ({@link DirectionMapping#emitWeighted}).
 *
 * <p><b>Not bit-exact</b> (opt-in, tolerance-tested): a pair at the cutoff boundary, or a direction within
 * float epsilon of the occlusion boundary, may be classified differently — both effects are tiny (boundary
 * pairs have near-empty caps; boundary directions are a handful out of hundreds of thousands), within SAS
 * discretization error. On a JVM without {@code jdk.incubator.vector} it falls back to the exact scalar
 * build + scalar scan, so there it is bit-exact. {@link #isVectorized()} reports the active path.
 *
 * <p>Prefer this over {@link FloatNumericalSurface} as the approximate/fast variant. Like the other
 * distinct surfaces it is deliberately not point-set-identical to CDK {@code NumericalSurface}.
 */
public class FloatNumericalSurfaceV2 extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new Float256WeightedDedupOcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the float build + float scan are in use; false on a JVM without the Vector API. */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public FloatNumericalSurfaceV2(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public FloatNumericalSurfaceV2(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                VECTOR_AVAILABLE
                        ? (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new FloatSimdDistanceCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get)
                        : (NeighborSourceFactory) (atoms, ax, ay, az, radius) ->
                                new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new Float256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access
    }
}
