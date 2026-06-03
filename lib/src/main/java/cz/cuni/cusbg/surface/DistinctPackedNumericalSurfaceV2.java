package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * V2 of {@link DistinctPackedNumericalSurface}: the same distinct-point surface - one point per
 * distinct surviving direction (~5.7x fewer points than the full-multiplicity surfaces, no need for
 * downstream sparsification), with <em>areas bit-for-bit identical</em> to
 * {@link FasterNumericalSurface} - but with three improvements that matter for code built on top of it:
 *
 * <ul>
 *   <li><b>SIMD scan.</b> Uses {@link Vectorized256WeightedDedupOcclusionScan} when the Vector API is
 *       available (matching {@link PackedNumericalSurface}'s vectorized path), falling back to the
 *       scalar {@link WeightedDedupOcclusionScan} otherwise. The original
 *       {@link DistinctPackedNumericalSurface} was scalar-only.</li>
 *   <li><b>Right-sized output buffer.</b> Pairs with {@link DistinctFlatSurfacePointStoreV2}, which
 *       sizes its flat buffer for distinct survival instead of full multiplicity (the original store
 *       over-allocated ~6-8x, defeating the surface's footprint advantage).</li>
 *   <li><b>Shared direction mapping.</b> Both scan paths resolve the distinct-direction collapse from
 *       one cached {@link DirectionMapping} instead of each scan carrying its own build.</li>
 * </ul>
 *
 * <p>The compute pipeline (cell-sorted pruned neighbour build, cached tessellation, cached VdW,
 * copy-free CSR access) mirrors {@link PackedNumericalSurface}; the only axis that differs is the
 * {@code (scan, store)} pair, which is exactly what makes this the distinct surface. This class is the
 * single assembly that wires the weighted scan to the deduplicating store, so the pairing invariant
 * (weighted scan &lt;-&gt; distinct store) cannot be mismatched by a caller.
 *
 * <p>Like {@link DistinctPackedNumericalSurface} this is deliberately NOT point-set-identical to CDK
 * {@code NumericalSurface} - it omits the coincident duplicates; the surviving point set equals what
 * sparsification at the exact-coincidence distance produces, and the areas are exact.
 */
public class DistinctPackedNumericalSurfaceV2 extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new Vectorized256WeightedDedupOcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the vectorized weighted dedup scan is in use; false on a JVM without the Vector API (scalar fallback). */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public DistinctPackedNumericalSurfaceV2(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DistinctPackedNumericalSurfaceV2(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
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
