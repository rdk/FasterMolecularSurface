package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.interfaces.IAtom;

import java.util.function.ToDoubleFunction;

/**
 * Like {@link SymmetricCellGridNeighborList}, but applies the exact per-pair occlusion cutoff while
 * building the adjacency: a neighbor {@code j} is kept for atom {@code i} only if their expanded
 * spheres overlap, {@code d(i,j) < R_i + R_j} (with {@code R = vdwRadius + solventRadius}).
 *
 * <p>This is output-preserving for the surface. A neighbor buries a tessellation point {@code p} iff
 * {@code diff.p > thresh}; since {@code max(diff.p) = |diff| = d} and {@code thresh = (d^2 + R_i^2 -
 * R_j^2)/(2 R_i)}, one has {@code d > thresh} exactly when {@code d < R_i + R_j}. So when
 * {@code d >= R_i + R_j} the neighbor satisfies {@code diff.p <= d <= thresh} for every {@code p} and
 * can never bury anything: dropping it cannot change which points survive, their coordinates, or the
 * areas. The cell grid's global cutoff ({@code 2*(maxRadius+solvent)}) is looser than this per-pair
 * bound and admits many such never-firing neighbors; filtering them here shrinks the per-atom neighbor
 * list the engine then feeds to the occlusion scan (the dominant cost), with no change to the result.
 *
 * <p>The grid cell size is still the global cutoff, so the 3x3x3 forward half-stencil remains
 * exhaustive for the largest possible kept pair ({@code R_i + R_j <= 2*maxRadius}). The neighbor SET is
 * a (correctness-preserving) subset of {@link CellGridNeighborList} / {@link NeighborList}, so it does
 * NOT match them as a set - only the final surface output matches, bit-for-bit. Same dense-grid
 * precondition as {@link CellGrid}.
 */
final class PrunedSymmetricCellGridNeighborList implements NeighborSource {

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    PrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                        double radius, double solventRadius) {
        this(atoms, ax, ay, az, radius, solventRadius, FasterNumericalSurface::getVdwRadius);
    }

    /**
     * @param vdwRadius van der Waals radius lookup (without solvent). The default constructor uses
     *        {@link FasterNumericalSurface#getVdwRadius}; a variant may pass {@link VdwRadiusCache#get}
     *        to memoize it per element symbol. Must return the same value either way (bit-exact).
     */
    PrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                        double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        CellGrid g = new CellGrid(ax, ay, az, radius);
        int n = ax.length;

        // per-atom expanded radius (same definition the engine uses for thresh)
        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // single distance pass over forward pairs, keeping only pairs whose expanded spheres overlap
        int[] degree = new int[n];
        IntArrayList edgeI = new IntArrayList(n * 8);
        IntArrayList edgeJ = new IntArrayList(n * 8);
        for (int i = 0; i < n; i++) {
            double xi = g.ax[i], yi = g.ay[i], zi = g.az[i], ri = expandedR[i];
            int bx = g.cx[i], by = g.cy[i], bz = g.cz[i];
            for (int dx = -1; dx <= 1; dx++) {
                int ix = bx + dx; if (ix < g.minX || ix >= g.minX + g.nx) continue;
                for (int dy = -1; dy <= 1; dy++) {
                    int iy = by + dy; if (iy < g.minY || iy >= g.minY + g.ny) continue;
                    for (int dz = -1; dz <= 1; dz++) {
                        int lin = dx * 9 + dy * 3 + dz;
                        if (lin < 0) continue;   // backward cell: recorded from the other side
                        int iz = bz + dz; if (iz < g.minZ || iz >= g.minZ + g.nz) continue;
                        int c = g.cellIndex(ix, iy, iz);
                        for (int p = g.cellStart[c], end = g.cellStart[c + 1]; p < end; p++) {
                            int j = g.cellAtoms[p];
                            if (lin == 0 && j <= i) continue;   // same cell: visit each pair once (j > i)
                            double x12 = g.ax[j] - xi, y12 = g.ay[j] - yi, z12 = g.az[j] - zi;
                            double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            double sumR = ri + expandedR[j];
                            if (d2 < sumR * sumR) {   // expanded spheres overlap -> j can occlude i (and vice versa)
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
