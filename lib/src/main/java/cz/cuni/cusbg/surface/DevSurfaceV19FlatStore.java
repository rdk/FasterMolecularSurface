package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * {@link DevSurfaceV18SortedCoords} with the flat point store (optimization B) instead of the
 * {@code Point3d}-per-point list store. The surface computation is identical to V18 (same neighbor
 * build, scan, tessellation - bit-for-bit); only point <em>storage</em> changes:
 * {@link FlatSurfacePointStore} keeps coordinates in one flat {@code double[]} and materializes
 * {@code Point3d} lazily.
 *
 * <p>The point of this variant is the {@link PackedSurfaceAccess} path: because the store is flat, a
 * caller that wants raw coordinates (e.g. p2rank, which immediately turns each point into its own
 * three-double {@code Point}) gets them via {@code surfacePointsXYZ()} <em>with no copy</em> and never
 * allocates a {@code Point3d} or the {@code Point3d[]}. Through the {@code Point3d}-valued
 * {@link MolecularSurface} accessors it behaves like V18 (lazily materializing on demand), so on its own
 * it is throughput-neutral here (the box is bandwidth-bound, not GC-bound; see docs); its value is the
 * zero-copy delivery path and the lower allocation footprint, not single-thread speed.
 *
 * <p>Output is bit-for-bit identical to {@link FasterNumericalSurface}. Scalar fallback when
 * {@code jdk.incubator.vector} is absent; {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV19FlatStore extends DevSurfaceV1Soa {

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

    public DevSurfaceV19FlatStore(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV19FlatStore(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                FlatSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access (inherited from V16)
    }
}
