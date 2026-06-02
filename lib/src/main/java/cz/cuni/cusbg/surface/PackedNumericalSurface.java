package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Production solvent-accessible-surface implementation optimized for callers that consume the surface
 * points as raw coordinates (e.g. building their own point objects) rather than as {@code Point3d}.
 *
 * <p>It is bit-for-bit identical to {@link FasterNumericalSurface} and CDK's {@code NumericalSurface}
 * (same tessellation, neighbour sets, occlusion test and arithmetic, same van der Waals radii including
 * the metal fallback), but it pairs the fastest compute pipeline (structure-of-arrays scratch,
 * cell-sorted pruned neighbour build, process-cached dedup SIMD occlusion scan, cached tessellation)
 * with a flat {@code double[]} point store and the zero-copy {@link PackedSurfaceAccess} delivery path.
 * A caller that reads {@link #surfacePointsXYZ()} / {@link #surfacePointCount()} never materializes a
 * {@code Point3d} per point; the legacy {@link MolecularSurface} accessors still work (materializing
 * {@code Point3d} lazily) for compatibility.
 *
 * <p>This is the stable, production-named form of the optimization ladder's flat-store champion (the
 * {@code DevSurfaceV*} classes are the experimental ladder rungs). Scalar fallback when
 * {@code jdk.incubator.vector} is unavailable; {@link #isVectorized()} reports the active scan path.
 */
public class PackedNumericalSurface extends DevSurfaceV1Soa {

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

    public PackedNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public PackedNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                FlatSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access
    }
}
