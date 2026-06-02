package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;

/**
 * Reusable per-thread scratch for {@link SoaNumericalSurface} (optimization A). Holds the engine's
 * transient per-build buffers - the extracted coordinate/radius arrays, the neighbor-query list, and
 * the per-neighbor {@code diff}/{@code thresh} arrays - so a worker thread reuses them across surface
 * builds instead of reallocating each time. Grown to the high-water mark and never shrunk.
 *
 * <p>Held in a {@code ThreadLocal} and used only when arena mode is enabled; a build runs {@code init()}
 * synchronously on one thread, so there is no concurrent access to a thread's scratch. The buffers are
 * fully overwritten before being read each build (coordinates/radii for {@code [0, n)}, the neighbor
 * list cleared, {@code diff}/{@code thresh} written for {@code [0, numNeighbors)}), so reuse cannot leak
 * stale values into a result. The coordinate arrays may be larger than the atom count; consumers that
 * need the count take it explicitly (see {@link CellGrid}, {@link PrunedSymmetricCellGridNeighborList}).
 *
 * <p>Does NOT hold the output point store (that is retained by each surface for its accessors, so it
 * cannot be reused across builds) nor the neighbor source's internal arrays (encapsulated in the
 * source).
 */
final class EngineScratch {

    private static final double[] EMPTY = new double[0];

    double[] ax = EMPTY, ay = EMPTY, az = EMPTY, atomRadius = EMPTY, atomRadius2 = EMPTY;
    final IntArrayList nbr = new IntArrayList(32);
    double[] diffX = new double[512], diffY = new double[512], diffZ = new double[512], thresh = new double[512];

    /** Ensure the per-atom arrays hold at least {@code n} entries. */
    void ensureAtoms(int n) {
        if (ax.length < n) {
            ax = new double[n]; ay = new double[n]; az = new double[n];
            atomRadius = new double[n]; atomRadius2 = new double[n];
        }
    }
}
