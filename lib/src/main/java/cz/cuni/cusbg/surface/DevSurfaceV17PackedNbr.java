package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step on top of {@link DevSurfaceV16DirectNbr}, targeting the cost the V16 profile
 * exposed. Once V16 removed the per-query neighbor copy, the construction's own edge buffering -
 * {@code edgeI.add(i); edgeJ.add(j)} into two HPPC {@code IntArrayList}s - became the dominant non-scan
 * hot method (~25% of CPU single-thread, ~34% at 16 threads, and most of the {@code int[]} allocation).
 *
 * <p>This variant swaps {@link PrunedSymmetricCellGridNeighborList} for
 * {@link PackedPrunedSymmetricCellGridNeighborList}: the same single distance pass, but each kept edge
 * is packed as an interleaved {@code (i, j)} pair into one cursor-managed {@code int[]} instead of two
 * {@code IntArrayList}s - one sequential write-stream, one capacity check per edge (not two), no per-add
 * wrapper overhead. Crucially it keeps the <em>single</em> distance pass, so it does not repeat the
 * mistake of {@link DevSurfaceV15LeanNbr} (which removed the buffers by doubling the distance pass and
 * regressed); the edge buffer stays, it is just made cheap.
 *
 * <p>Everything else matches V16 - copy-free CSR neighbor access ({@code directNeighbors}), cached VdW
 * radii (engine + source), process-cached 256-bit dedup scan, last-occluder hint, cached tessellation,
 * {@link ListSurfacePointStore}. Output is bit-for-bit identical to {@link FasterNumericalSurface}
 * (same neighbor set, same order, same radii). Scalar fallback when {@code jdk.incubator.vector} is
 * absent; {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV17PackedNbr extends DevSurfaceV1Soa {

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

    public DevSurfaceV17PackedNbr(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV17PackedNbr(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PackedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access (inherited from V16)
    }
}
