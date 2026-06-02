package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;

/**
 * Cell-grid neighbor index that exploits the symmetry of the neighbor relation: j is within the
 * cutoff of i iff i is within the cutoff of j. {@link CellGridNeighborList} answers each per-atom
 * query independently and so evaluates every pair's distance <em>twice</em> (once from each endpoint);
 * profiling showed that query had become a large fraction of surface construction once the occlusion
 * scan was optimized. This variant instead computes each unordered pair's distance <em>once</em>, in a
 * single up-front pass over the {@link CellGrid}, and materializes a CSR adjacency
 * ({@link #adjStart} + {@link #adj}); queries are then a plain array copy.
 *
 * <p>The pass uses a half-stencil: for each atom it scans only the 13 "forward" neighbor cells
 * ({@code dx*9 + dy*3 + dz > 0}) plus, within its own cell, partners with a higher atom index. Each
 * cell-cell (and same-cell) pair is therefore visited from exactly one side; on a hit both directions
 * are recorded ({@code i->j} and {@code j->i}). The resulting neighbor <em>set</em> for each atom is
 * identical to {@link CellGridNeighborList} / {@link NeighborList} (only the within-list order differs,
 * which the occlusion test is invariant to), so the surface output is bit-for-bit unchanged.
 *
 * <p>Same dense-grid precondition as {@link CellGrid}: use the hash-box {@link NeighborList} for
 * sparse/outlier inputs.
 */
final class SymmetricCellGridNeighborList implements NeighborSource {

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    SymmetricCellGridNeighborList(double[] ax, double[] ay, double[] az, double radius) {
        CellGrid g = new CellGrid(ax, ay, az, radius);
        int n = ax.length;

        // single distance pass over forward pairs: record each in-range edge once, count both degrees
        int[] degree = new int[n];
        IntArrayList edgeI = new IntArrayList(n * 8);
        IntArrayList edgeJ = new IntArrayList(n * 8);
        for (int i = 0; i < n; i++) {
            double xi = g.ax[i], yi = g.ay[i], zi = g.az[i];
            int bx = g.cx[i], by = g.cy[i], bz = g.cz[i];
            for (int dx = -1; dx <= 1; dx++) {
                int ix = bx + dx; if (ix < g.minX || ix >= g.minX + g.nx) continue;
                for (int dy = -1; dy <= 1; dy++) {
                    int iy = by + dy; if (iy < g.minY || iy >= g.minY + g.ny) continue;
                    for (int dz = -1; dz <= 1; dz++) {
                        int lin = dx * 9 + dy * 3 + dz;
                        if (lin < 0) continue;   // backward cell: the pair is recorded from the other side
                        int iz = bz + dz; if (iz < g.minZ || iz >= g.minZ + g.nz) continue;
                        int c = g.cellIndex(ix, iy, iz);
                        for (int p = g.cellStart[c], end = g.cellStart[c + 1]; p < end; p++) {
                            int j = g.cellAtoms[p];
                            if (lin == 0 && j <= i) continue;   // same cell: visit each pair once (j > i)
                            double x12 = g.ax[j] - xi, y12 = g.ay[j] - yi, z12 = g.az[j] - zi;
                            if (x12 * x12 + y12 * y12 + z12 * z12 < g.cutoff2) {
                                edgeI.add(i); edgeJ.add(j);
                                degree[i]++; degree[j]++;
                            }
                        }
                    }
                }
            }
        }

        // build CSR adjacency from the recorded edges (both directions)
        adjStart = new int[n + 1];
        for (int i = 0; i < n; i++) adjStart[i + 1] = adjStart[i] + degree[i];
        adj = new int[adjStart[n]];
        int[] cursor = new int[n];
        System.arraycopy(adjStart, 0, cursor, 0, n);
        int[] ei = edgeI.buffer, ej = edgeJ.buffer;
        int edges = edgeI.size();
        for (int e = 0; e < edges; e++) {
            int i = ei[e], j = ej[e];
            adj[cursor[i]++] = j;
            adj[cursor[j]++] = i;
        }
    }

    @Override
    public void getNeighborsInto(int i, IntArrayList out) {
        for (int p = adjStart[i], end = adjStart[i + 1]; p < end; p++) out.add(adj[p]);
    }
}
