package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.interfaces.IAtom;

import java.util.Arrays;
import java.util.function.ToDoubleFunction;

/**
 * Same neighbor set and CSR output as {@link PackedPrunedSymmetricCellGridNeighborList} (exact per-pair
 * occlusion cutoff, single distance pass, packed {@code int[]} edge buffer), but it removes the random
 * coordinate gather that profiling flagged as the distance pass's bottleneck.
 *
 * <p>Profiling V17 ({@link DevSurfaceV17PackedNbr}) showed the neighbor build was the top hot method at
 * 16 threads (~47%), with ~87% of its samples on the inner cell-scan loop: for each candidate {@code p}
 * in a cell it reads the candidate's coordinates at the candidate's <em>atom</em> index
 * {@code j = cellAtoms[p]} - {@code ax[j], ay[j], az[j], expandedR[j]} - a 4-way <em>random</em> gather,
 * since {@code j} is in cell order, not atom order. The build had become memory-latency-bound on that
 * gather, not on arithmetic or buffering.
 *
 * <p>This variant builds, once, a cell-sorted interleaved coordinate array
 * {@code coords[4p .. 4p+3] = (x, y, z, expandedR)} of {@code cellAtoms[p]} (a single gather pass over
 * the atoms), then the distance pass reads {@code coords} <em>sequentially</em> as {@code p} increments
 * within a cell - one streaming read (4 doubles per candidate, sequential) instead of four random loads.
 * The outer atom keeps its atom-order coordinate (sequential in {@code i}). Each candidate's coordinates
 * are visited once per appearance across neighbor cells (many times); the cell-sorted copy is built once
 * (n gathers), so the trade strongly favors the sorted layout.
 *
 * <p>Bit-for-bit identical to {@link PackedPrunedSymmetricCellGridNeighborList} (hence to
 * {@link FasterNumericalSurface}): {@code coords[4p] == ax[cellAtoms[p]]} exactly, the comparisons are
 * the same, and edges are recorded in the same outer-{@code i} / inner-{@code p} order. Implements
 * {@link DirectNeighborSource} for copy-free access. Same dense-grid precondition as {@link CellGrid}.
 */
final class CoordSortedPrunedSymmetricCellGridNeighborList implements DirectNeighborSource {

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    CoordSortedPrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                                   double radius, double solventRadius) {
        this(atoms, ax, ay, az, radius, solventRadius, FasterNumericalSurface::getVdwRadius);
    }

    /**
     * @param vdwRadius van der Waals radius lookup (without solvent). Pass {@link VdwRadiusCache#get} to
     *        memoize it per element symbol (optimization D); must return the same value either way.
     */
    CoordSortedPrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                                   double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        int n = atoms.length;   // not ax.length: the coordinate arrays may be oversized (arena-reused)
        CellGrid g = new CellGrid(ax, ay, az, n, radius);

        // per-atom expanded radius (atom order; used for the outer atom's r and to seed the sorted copy)
        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // cell-sorted, interleaved candidate coordinates: coords[4p..4p+3] = (x,y,z,expandedR) of the atom
        // at cell-order position p. Built once with a gather; the distance pass then reads it sequentially.
        double[] coords = new double[4 * n];
        for (int p = 0; p < n; p++) {
            int j = g.cellAtoms[p];
            int b = p << 2;
            coords[b] = g.ax[j]; coords[b + 1] = g.ay[j]; coords[b + 2] = g.az[j]; coords[b + 3] = expandedR[j];
        }

        // single distance pass over forward pairs, reading candidate coords sequentially from `coords`
        int[] degree = new int[n];
        int[] edges = new int[Math.max(16, n * 16)];
        int ec = 0;   // number of ints written (2 per edge)
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
                            int b = p << 2;
                            double x12 = coords[b] - xi, y12 = coords[b + 1] - yi, z12 = coords[b + 2] - zi;
                            double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            double sumR = ri + coords[b + 3];
                            if (d2 < sumR * sumR) {   // expanded spheres overlap -> j can occlude i (and vice versa)
                                if (ec + 2 > edges.length) edges = Arrays.copyOf(edges, edges.length * 2);
                                edges[ec++] = i; edges[ec++] = j;
                                degree[i]++; degree[j]++;
                            }
                        }
                    }
                }
            }
        }

        // build CSR adjacency from the packed edges (both directions)
        adjStart = new int[n + 1];
        for (int i = 0; i < n; i++) adjStart[i + 1] = adjStart[i] + degree[i];
        adj = new int[adjStart[n]];
        int[] cursor = new int[n];
        System.arraycopy(adjStart, 0, cursor, 0, n);
        for (int e = 0; e < ec; e += 2) {
            int i = edges[e], j = edges[e + 1];
            adj[cursor[i]++] = j;
            adj[cursor[j]++] = i;
        }
    }

    @Override
    public void getNeighborsInto(int i, IntArrayList out) {
        for (int p = adjStart[i], end = adjStart[i + 1]; p < end; p++) out.add(adj[p]);
    }

    // --- DirectNeighborSource: copy-free access to the CSR adjacency (used by DevSurfaceV18SortedCoords) ---

    @Override
    public int[] adjacency() { return adj; }

    @Override
    public int neighborStart(int i) { return adjStart[i]; }

    @Override
    public int neighborEnd(int i) { return adjStart[i + 1]; }
}
