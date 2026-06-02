package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step on top of {@link DevSurfaceV11CachedTess}, targeting the neighbor-list cost the
 * V11 profile flagged, but via "less work" rather than "less allocation" (the lesson from the
 * {@link DevSurfaceV15LeanNbr} regression):
 * <ul>
 *   <li><b>#2 Copy-free CSR neighbor access.</b> V11's {@link PrunedSymmetricCellGridNeighborList}
 *       stores every atom's neighbors in a single CSR array, but the {@link NeighborSource} contract
 *       made the engine copy each atom's slice into a reused {@code IntArrayList} per query, element by
 *       element ({@code IntArrayList.add} was ~15% of single-thread CPU, ~2/3 of it this per-query copy).
 *       The source now also implements {@link DirectNeighborSource}, and this variant enables the
 *       engine's {@code directNeighbors} path so the neighbor pre-pass reads the CSR array directly with
 *       no per-query copy. No added arithmetic (unlike V15's two-pass build), so it cannot regress the
 *       build the way trading compute for allocation did.</li>
 *   <li><b>#3 Cached VdW radii (engine + neighbor source).</b> The van der Waals radius is a pure
 *       function of element symbol but was resolved per atom via a {@code toLowerCase} + periodic-table
 *       lookup twice per build; both now go through {@link VdwRadiusCache}.</li>
 * </ul>
 *
 * <p>Everything else matches V11 - the same single-pass {@link PrunedSymmetricCellGridNeighborList}
 * construction, process-cached 256-bit dedup scan, last-occluder hint, cached tessellation, and
 * {@link ListSurfacePointStore} - so this isolates the copy-free-access + cached-VdW changes. Output is
 * bit-for-bit identical to {@link FasterNumericalSurface} (same neighbor set, same order, same radii;
 * only a per-query copy and a memoized lookup are removed). Scalar fallback when
 * {@code jdk.incubator.vector} is absent; {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV16DirectNbr extends DevSurfaceV1Soa {

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

    public DevSurfaceV16DirectNbr(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV16DirectNbr(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access (#2)
    }
}
