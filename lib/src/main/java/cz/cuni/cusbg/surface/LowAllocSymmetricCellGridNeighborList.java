package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;

/**
 * Low-allocation rebuild of {@link SymmetricCellGridNeighborList}. Same result (each unordered
 * neighbor pair's distance is symmetric, both directions recorded in a CSR adjacency, queries served
 * as array copies, identical neighbor sets) but it eliminates the intermediate edge buffer that
 * profiling flagged.
 *
 * <p>{@link SymmetricCellGridNeighborList} discovers every forward edge in a single distance pass and
 * buffers them in two growable {@code IntArrayList}s before assembling the CSR; those buffers
 * (~16 bytes per edge, reallocated by doubling) were the {@code int[]} churn in the allocation profile.
 * This variant trades that memory for arithmetic: it makes <em>two</em> passes over the 13-cell forward
 * half-stencil - a counting pass that only accumulates per-atom degrees, then a fill pass that writes
 * straight into the CSR {@code adj} array - so nothing beyond {@code degree[]} and the final CSR is
 * ever allocated.
 *
 * <p><b>Tradeoff (measured on 4HHB, tess 2, GraalVM 25):</b> ~63% less allocation per build, but the
 * repeated distance pass makes it ~2-7% slower than {@link SymmetricCellGridNeighborList} (the gap is
 * largest at low tessellation, where the neighbor build is a bigger fraction of total time). Prefer
 * this only where GC pressure / throughput at scale matters more than single-run latency;
 * {@link SymmetricCellGridNeighborList} is the latency-optimal one.
 *
 * <p>Output is bit-for-bit identical to {@link SymmetricCellGridNeighborList} /
 * {@link CellGridNeighborList} / {@link NeighborList}. Same dense-grid precondition as {@link CellGrid}.
 */
final class LowAllocSymmetricCellGridNeighborList implements NeighborSource {

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    LowAllocSymmetricCellGridNeighborList(double[] ax, double[] ay, double[] az, double radius) {
        CellGrid g = new CellGrid(ax, ay, az, radius);
        int n = ax.length;

        // pass 1: count each atom's degree (every in-range forward pair contributes to both endpoints)
        int[] degree = new int[n];
        forwardPairs(g, n, (i, j) -> { degree[i]++; degree[j]++; });

        // CSR offsets
        adjStart = new int[n + 1];
        for (int i = 0; i < n; i++) adjStart[i + 1] = adjStart[i] + degree[i];
        adj = new int[adjStart[n]];

        // pass 2: fill the CSR adjacency directly (cursor[i] walks atom i's slice)
        int[] cursor = new int[n];
        System.arraycopy(adjStart, 0, cursor, 0, n);
        forwardPairs(g, n, (i, j) -> { adj[cursor[i]++] = j; adj[cursor[j]++] = i; });
    }

    /** Functional sink for an in-range forward edge {@code (i, j)} with {@code i < j} discovery order. */
    private interface EdgeSink { void edge(int i, int j); }

    /** Visit each in-range forward pair once (13-cell forward half-stencil + same-cell j>i), calling {@code sink}. */
    private static void forwardPairs(CellGrid g, int n, EdgeSink sink) {
        for (int i = 0; i < n; i++) {
            double xi = g.ax[i], yi = g.ay[i], zi = g.az[i];
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
                            if (x12 * x12 + y12 * y12 + z12 * z12 < g.cutoff2) sink.edge(i, j);
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
