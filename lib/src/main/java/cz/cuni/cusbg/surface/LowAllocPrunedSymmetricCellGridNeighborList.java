package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.interfaces.IAtom;

import java.util.function.ToDoubleFunction;

/**
 * Low-allocation rebuild of {@link PrunedSymmetricCellGridNeighborList}. Same neighbor set (the exact
 * per-pair occlusion cutoff {@code d(i,j) < R_i + R_j}, both directions in a CSR adjacency, queries
 * served as array copies) but it eliminates the two intermediate {@code IntArrayList} edge buffers that
 * profiling flagged as the dominant allocation.
 *
 * <p>{@link PrunedSymmetricCellGridNeighborList} discovers every kept forward edge in a single distance
 * pass and buffers {@code (i, j)} into two growable {@code IntArrayList(n*8)} before assembling the CSR;
 * those two buffers were the bulk of the {@code int[]} allocation (61% of all allocation at 16 threads,
 * where the neighbor build is the hottest method). This variant trades that memory for arithmetic: it
 * makes <em>two</em> passes over the 13-cell forward half-stencil applying the same cutoff each time - a
 * counting pass that only accumulates per-atom degrees, then a fill pass that writes straight into the
 * CSR {@code adj} array - so nothing beyond {@code degree[]}, {@code expandedR[]}, {@code cursor[]} and
 * the final CSR is ever allocated. This mirrors what {@link LowAllocSymmetricCellGridNeighborList} does
 * for the unpruned source, now with the per-pair cutoff.
 *
 * <p>The fill pass rediscovers pairs in the same forward-stencil order as the counting pass, so each
 * atom's {@code adj} slice is filled in the identical order the single-pass variant produced; the
 * occlusion verdict is order-independent regardless (a point is buried iff <em>any</em> neighbor buries
 * it), so the surface is bit-for-bit identical to {@link PrunedSymmetricCellGridNeighborList} (hence to
 * {@link FasterNumericalSurface}). The neighbor SET is the same correctness-preserving subset, so like
 * the pruned source it does NOT match {@link CellGridNeighborList} / {@link NeighborList} as a set, only
 * the final surface output. Same dense-grid precondition as {@link CellGrid}.
 */
final class LowAllocPrunedSymmetricCellGridNeighborList implements NeighborSource {

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    LowAllocPrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                                double radius, double solventRadius) {
        this(atoms, ax, ay, az, radius, solventRadius, FasterNumericalSurface::getVdwRadius);
    }

    /**
     * @param vdwRadius van der Waals radius lookup (without solvent). Pass {@link VdwRadiusCache#get} to
     *        memoize it per element symbol (optimization D); must return the same value either way.
     */
    LowAllocPrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                                double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        int n = atoms.length;   // not ax.length: the coordinate arrays may be oversized (arena-reused)
        CellGrid g = new CellGrid(ax, ay, az, n, radius);

        // per-atom expanded radius (same definition the engine uses for thresh)
        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // pass 1: count each atom's degree (every kept forward pair contributes to both endpoints)
        int[] degree = new int[n];
        prunedForwardPairs(g, n, expandedR, (i, j) -> { degree[i]++; degree[j]++; });

        // CSR offsets
        adjStart = new int[n + 1];
        for (int i = 0; i < n; i++) adjStart[i + 1] = adjStart[i] + degree[i];
        adj = new int[adjStart[n]];

        // pass 2: fill the CSR adjacency directly (cursor[i] walks atom i's slice)
        int[] cursor = new int[n];
        System.arraycopy(adjStart, 0, cursor, 0, n);
        prunedForwardPairs(g, n, expandedR, (i, j) -> { adj[cursor[i]++] = j; adj[cursor[j]++] = i; });
    }

    /** Functional sink for a kept forward edge {@code (i, j)} with {@code i < j} discovery order. */
    private interface EdgeSink { void edge(int i, int j); }

    /**
     * Visit each kept forward pair once (13-cell forward half-stencil + same-cell j>i), applying the
     * per-pair cutoff {@code d2 < (R_i + R_j)^2}, calling {@code sink}. Identical pair-discovery logic to
     * {@link PrunedSymmetricCellGridNeighborList}'s single pass.
     */
    private static void prunedForwardPairs(CellGrid g, int n, double[] expandedR, EdgeSink sink) {
        for (int i = 0; i < n; i++) {
            double xi = g.ax[i], yi = g.ay[i], zi = g.az[i], ri = expandedR[i];
            int bx = g.cx[i], by = g.cy[i], bz = g.cz[i];
            for (int dx = -1; dx <= 1; dx++) {
                int ix = bx + dx; if (ix < g.minX || ix >= g.minX + g.nx) continue;
                for (int dy = -1; dy <= 1; dy++) {
                    int iy = by + dy; if (iy < g.minY || iy >= g.minY + g.ny) continue;
                    for (int dz = -1; dz <= 1; dz++) {
                        int lin = dx * 9 + dy * 3 + dz;
                        if (lin < 0) continue;   // backward cell: visited from the other side
                        int iz = bz + dz; if (iz < g.minZ || iz >= g.minZ + g.nz) continue;
                        int c = g.cellIndex(ix, iy, iz);
                        for (int p = g.cellStart[c], end = g.cellStart[c + 1]; p < end; p++) {
                            int j = g.cellAtoms[p];
                            if (lin == 0 && j <= i) continue;   // same cell: visit each pair once (j > i)
                            double x12 = g.ax[j] - xi, y12 = g.ay[j] - yi, z12 = g.az[j] - zi;
                            double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            double sumR = ri + expandedR[j];
                            if (d2 < sumR * sumR) sink.edge(i, j);   // expanded spheres overlap
                        }
                    }
                }
            }
        }
    }

    @Override
    public void getNeighborsInto(int i, IntArrayList out) {
        for (int p = adjStart[i], end = adjStart[i + 1]; p < end; p++) out.add(adj[p]);
    }
}
