package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.interfaces.IAtom;

import java.util.Arrays;
import java.util.function.ToDoubleFunction;

/**
 * Same neighbor set and CSR output as {@link PrunedSymmetricCellGridNeighborList} (exact per-pair
 * occlusion cutoff {@code d(i,j) < R_i + R_j}, both directions, single distance pass), but it records
 * the kept forward edges in one cursor-managed primitive {@code int[]} instead of two HPPC
 * {@link IntArrayList}s.
 *
 * <p>Profiling V16 ({@link DevSurfaceV16DirectNbr}) showed the edge buffering -
 * {@code edgeI.add(i); edgeJ.add(j)} - had become the dominant non-scan cost (~34% of CPU at 16
 * threads, and most of the {@code int[]} allocation). Two {@code IntArrayList}s mean two {@code add}
 * calls per edge, two bounds-checks, and two separate growing write-streams. This variant packs each
 * kept edge as an interleaved {@code (i, j)} pair into a single {@code int[]} grown by doubling: one
 * sequential write-stream, one capacity check per edge (not two), and no per-add wrapper overhead. It
 * keeps the <em>single</em> distance pass, so unlike {@link LowAllocPrunedSymmetricCellGridNeighborList}
 * (the two-pass {@code DevSurfaceV15LeanNbr} build, which regressed) it adds no extra arithmetic.
 *
 * <p>The CSR adjacency is assembled from the packed edges in the same forward-discovery order the
 * two-buffer variant produced, so each atom's {@code adj} slice is identical; the occlusion verdict is
 * order-independent regardless. Output is bit-for-bit identical to
 * {@link PrunedSymmetricCellGridNeighborList} (hence to {@link FasterNumericalSurface}). Implements
 * {@link DirectNeighborSource} for copy-free access. Same dense-grid precondition as {@link CellGrid}.
 */
final class PackedPrunedSymmetricCellGridNeighborList implements DirectNeighborSource {

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    PackedPrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                              double radius, double solventRadius) {
        this(atoms, ax, ay, az, radius, solventRadius, FasterNumericalSurface::getVdwRadius);
    }

    /**
     * @param vdwRadius van der Waals radius lookup (without solvent). Pass {@link VdwRadiusCache#get} to
     *        memoize it per element symbol (optimization D); must return the same value either way.
     */
    PackedPrunedSymmetricCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                              double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        int n = atoms.length;   // not ax.length: the coordinate arrays may be oversized (arena-reused)
        CellGrid g = new CellGrid(ax, ay, az, n, radius);

        // per-atom expanded radius (same definition the engine uses for thresh)
        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // single distance pass over forward pairs, packing each kept edge as interleaved (i, j) into one
        // int[] grown by doubling; degrees counted inline. n*16 ints == n*8 edges, same headroom the
        // two-IntArrayList variant pre-sized to.
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
                            double x12 = g.ax[j] - xi, y12 = g.ay[j] - yi, z12 = g.az[j] - zi;
                            double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            double sumR = ri + expandedR[j];
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

    // --- DirectNeighborSource: copy-free access to the CSR adjacency (used by DevSurfaceV17PackedNbr) ---

    @Override
    public int[] adjacency() { return adj; }

    @Override
    public int neighborStart(int i) { return adjStart[i]; }

    @Override
    public int neighborEnd(int i) { return adjStart[i + 1]; }
}
