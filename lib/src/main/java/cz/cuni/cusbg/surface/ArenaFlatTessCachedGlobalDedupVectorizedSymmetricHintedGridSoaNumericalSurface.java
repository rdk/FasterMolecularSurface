package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization A applied to {@link FlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface}
 * (the flat-output, cached-VdW variant): the engine's transient per-build scratch is reused from a
 * per-thread {@link EngineScratch} arena instead of reallocated per surface. The full stack with both
 * the flat point store (B) and the arena (A).
 *
 * <p>Output-identical (bit-for-bit vs {@link FasterNumericalSurface}); only allocation changes. As with
 * the list-store arena variant, the flat point store itself is not pooled (it is retained per surface),
 * so the arena pools the engine-side transient scratch (coordinates/radii, neighbor list, diff/thresh).
 */
public class ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface extends SoaNumericalSurface {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new GlobalDedupVectorized256OcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                FlatSurfacePointStore::new,
                true);
    }
}
