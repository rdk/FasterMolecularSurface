package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step on top of {@link DevSurfaceV17PackedNbr}, targeting the random coordinate gather the
 * V17 profile exposed. Once the edge buffer was made cheap (V17), the neighbor build's distance pass was
 * the top hot method at 16 threads (~47%), and ~87% of its samples were on the inner cell-scan loop's
 * candidate-coordinate read: {@code ax[j], ay[j], az[j], expandedR[j]} at the candidate's atom index
 * {@code j = cellAtoms[p]} - a 4-way random gather, since {@code j} runs in cell order, not atom order.
 * The build was memory-latency-bound on that gather.
 *
 * <p>This variant swaps {@link PackedPrunedSymmetricCellGridNeighborList} for
 * {@link CoordSortedPrunedSymmetricCellGridNeighborList}: it builds a cell-sorted interleaved coordinate
 * array once (a single gather pass) so the distance pass reads candidate coordinates <em>sequentially</em>
 * as it scans a cell, turning the per-candidate random gather into one streaming read. Same single
 * distance pass, same packed edge buffer.
 *
 * <p>Everything else matches V17 - copy-free CSR neighbor access ({@code directNeighbors}), cached VdW
 * radii, process-cached 256-bit dedup scan, last-occluder hint, cached tessellation,
 * {@link ListSurfacePointStore}. Output is bit-for-bit identical to {@link FasterNumericalSurface}
 * (the cell-sorted coordinates hold the same values, the comparisons and edge order are unchanged).
 * Scalar fallback when {@code jdk.incubator.vector} is absent; {@link #isVectorized()} reports the path.
 */
public class DevSurfaceV18SortedCoords extends DevSurfaceV1Soa {

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

    public DevSurfaceV18SortedCoords(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV18SortedCoords(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access (inherited from V16)
    }
}
