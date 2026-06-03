package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Single-precision distinct-point solvent-accessible surface: the {@code float} analog of
 * {@link DistinctPackedNumericalSurfaceV2}. It runs the identical pipeline (cell-sorted pruned
 * neighbour build, cached tessellation, cached VdW, copy-free CSR access, distinct weighted dedup),
 * but the occlusion VERDICT is computed in {@code float} via {@link Float256WeightedDedupOcclusionScan}
 * - 8 SIMD lanes instead of 4 and half the neighbor-scratch memory traffic - attacking the scan that
 * profiling showed is ~72% of V2's CPU at tessellation level 4 (the throughput-bound hot spot).
 *
 * <p><b>Not bit-exact.</b> Point POSITIONS and per-atom areas are still computed in {@code double}
 * (the scan only narrows the comparison, and the coordinate subtraction stays in {@code double} to
 * avoid cancellation), so every surviving point's coordinates are identical to V2's. But a tessellation
 * point within {@code float} epsilon of the occlusion boundary may flip survival, so the surviving SET
 * - and therefore the area - can differ from {@link FasterNumericalSurface} by a small tolerance, well
 * within SAS tessellation discretization error. This surface ships with tolerance tests rather than the
 * bit-exact equivalence harness the other variants use.
 *
 * <p>Falls back to the exact scalar {@code double} scan ({@link WeightedDedupOcclusionScan}) when the
 * Vector API is unavailable - single-precision scalar offers no benefit, and the exact path trivially
 * satisfies the tolerance. {@link #isVectorized()} reports which path is active.
 *
 * <p>Like the other distinct surfaces this is deliberately NOT point-set-identical to CDK
 * {@code NumericalSurface}: it omits the coincident duplicates and may drop/keep a few boundary points.
 */
public class FloatNumericalSurface extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new Float256WeightedDedupOcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the float SIMD scan is in use; false on a JVM without the Vector API (exact scalar fallback). */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public FloatNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public FloatNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
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
