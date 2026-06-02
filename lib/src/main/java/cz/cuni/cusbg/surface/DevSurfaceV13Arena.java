package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization A applied to {@link DevSurfaceV11CachedTess}
 * (the previous best, list-based point store): the engine's transient per-build scratch
 * (coordinate/radius arrays, neighbor list, diff/thresh) is reused from a per-thread {@link EngineScratch}
 * arena instead of reallocated per surface.
 *
 * <p>Output-identical to the non-arena variant (only allocation changes), so bit-for-bit vs
 * {@link FasterNumericalSurface}. The arena's main benefit is at 16-thread scale, where it cuts the
 * per-build allocation that competes for memory bandwidth. The output point store and the neighbor
 * source's internal arrays are not pooled (the store is retained per surface; the source encapsulates
 * its arrays), so this pools the engine-side transient scratch only.
 */
public class DevSurfaceV13Arena extends DevSurfaceV1Soa {

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

    public DevSurfaceV13Arena(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV13Arena(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                true);
    }
}
