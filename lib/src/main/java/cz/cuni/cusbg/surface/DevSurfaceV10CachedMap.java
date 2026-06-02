package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 9: {@link DevSurfaceV9Dedup} with the
 * distinct-direction mapping cached process-wide instead of rebuilt per surface.
 *
 * <p>Same full stack (SoA, per-pair occlusion cutoff, last-occluder hint, 256-bit dedup SIMD scan) but
 * the scan is {@link GlobalDedupVectorized256OcclusionScan}, which resolves the tessellation's
 * distinct-direction mapping from a process-wide cache keyed by tessellation level. The per-build dedup
 * variant rebuilt that mapping (a {@code String}-keyed {@code HashMap}) on every surface; profiling at
 * 16 threads showed that rebuild was the largest allocation source and a scaling tax once the dedup had
 * made the scan cheap. The mapping is identical for every build at a level and immutable, so it is
 * built once and shared read-only.
 *
 * <p>The dedup scan keeps a tiny per-build verdict scratch, so it is created fresh per surface (the
 * shared state is only the immutable cached mapping). When {@code jdk.incubator.vector} is absent the
 * scan is unavailable and this variant falls back to the scalar {@link OcclusionScan#LAST_OCCLUDER_FIRST}
 * (pruning still applies; dedup does not). Output is bit-for-bit identical to
 * {@link FasterNumericalSurface}. {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV10CachedMap extends DevSurfaceV1Soa {

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

    public DevSurfaceV10CachedMap(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV10CachedMap(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST);
    }
}
