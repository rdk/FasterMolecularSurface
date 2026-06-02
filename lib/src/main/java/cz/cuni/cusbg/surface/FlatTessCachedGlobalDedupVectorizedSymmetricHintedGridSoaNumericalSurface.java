package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization steps D + B on top of {@link TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface}:
 * <ul>
 *   <li><b>B - flat point output.</b> Surface points are stored in a flat {@code double[]}
 *       ({@link FlatSurfacePointStore}) instead of one {@code Point3d} object per point;
 *       {@code Point3d} are materialized lazily only when an object accessor is called. Removes the
 *       per-point allocation (a top allocation source), which mainly helps GC/bandwidth at scale.</li>
 *   <li><b>D - cached VdW radii.</b> The pruned neighbor source resolves van der Waals radii through
 *       {@link VdwRadiusCache} (memoized per element symbol) instead of a per-atom periodic-table
 *       lookup.</li>
 * </ul>
 *
 * <p>Stacks the full engine: pruned per-pair occlusion cutoff, last-occluder hint, process-cached
 * 256-bit dedup scan, and cached tessellation. Both B and D are bit-for-bit output-preserving (same
 * coordinates/areas; only storage and a memoized lookup change), so the result is identical to
 * {@link FasterNumericalSurface}. Scalar fallback when {@code jdk.incubator.vector} is absent.
 */
public class FlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface extends SoaNumericalSurface {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new GlobalDedupVectorized256OcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the dedup SIMD scan is in use; false if this JVM lacks the Vector API (scalar fallback). */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public FlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public FlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                FlatSurfacePointStore::new);
    }
}
