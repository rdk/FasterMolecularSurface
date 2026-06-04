package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.interfaces.IAtom;

import java.util.Arrays;
import java.util.function.ToDoubleFunction;

/**
 * Backlog idea A2: the same neighbor set, CSR output, and cell-sorted sequential distance pass as
 * {@link CoordSortedPrunedSymmetricCellGridNeighborList}, but built over a {@link TightCellGrid} with
 * cells half the cutoff ({@code cellsPerCutoff = 2}) searched with a correspondingly larger ±2 (5×5×5)
 * stencil. Finer cells pack fewer atoms, so the distance pass — the build that profiling puts at ~47% of
 * CPU at 16 threads / tess 2, the regime p2rank runs in — visits a smaller candidate volume.
 *
 * <p><b>Bit-for-bit identical neighbor set</b> to the frozen ±1 build (hence to {@link FasterNumericalSurface}'s
 * areas): both the 3×3×3-over-cutoff-cells and 5×5×5-over-half-cells searches enclose every pair within the
 * cutoff, and both prune with the <em>identical</em> {@code d2 < (Ri+Rj)²} test on the identical coordinates,
 * so the surviving pairs are exactly the same. The symmetric forward-only pairing generalizes to any stencil
 * radius {@code K}: each offset {@code o} and its negation {@code -o} get opposite-sign linear indices
 * {@code lin = dx·S² + dy·S + dz} ({@code S = 2K+1}), so every cross-cell pair is recorded once (from the
 * {@code lin>0} side) and same-cell pairs once ({@code lin==0, j>i}). Implements {@link DirectNeighborSource}
 * for copy-free access. Same dense-grid precondition as {@link CellGrid}.
 */
final class TightCellGridNeighborList implements DirectNeighborSource {

    /** Cells per cutoff length: cell side = cutoff/2, searched with a ±2 (5×5×5) stencil. */
    private static final int CELLS_PER_CUTOFF = 2;

    private final int[] adjStart;   // prefix sums, length n+1
    private final int[] adj;        // neighbor indices grouped by atom, length 2*edgeCount

    TightCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                              double radius, double solventRadius) {
        this(atoms, ax, ay, az, radius, solventRadius, FasterNumericalSurface::getVdwRadius);
    }

    /**
     * @param vdwRadius van der Waals radius lookup (without solvent); pass {@link VdwRadiusCache#get} to
     *        memoize per element symbol. Must return the same value either way.
     */
    TightCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                              double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        int n = atoms.length;   // not ax.length: coordinate arrays may be oversized (arena-reused)
        final int K = CELLS_PER_CUTOFF;
        final int S = 2 * K + 1;
        TightCellGrid g = new TightCellGrid(ax, ay, az, n, radius, K);

        // per-atom expanded radius (atom order)
        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // cell-sorted interleaved candidate coords (V18's sequential-read win), built once with a gather
        double[] coords = new double[4 * n];
        for (int p = 0; p < n; p++) {
            int j = g.cellAtoms[p];
            int b = p << 2;
            coords[b] = g.ax[j]; coords[b + 1] = g.ay[j]; coords[b + 2] = g.az[j]; coords[b + 3] = expandedR[j];
        }

        // single distance pass over forward pairs in the ±K stencil, reading candidate coords sequentially
        int[] degree = new int[n];
        int[] edges = new int[Math.max(16, n * 16)];
        int ec = 0;
        for (int i = 0; i < n; i++) {
            double xi = g.ax[i], yi = g.ay[i], zi = g.az[i], ri = expandedR[i];
            int bx = g.cx[i], by = g.cy[i], bz = g.cz[i];
            for (int dx = -K; dx <= K; dx++) {
                int ix = bx + dx; if (ix < g.minX || ix >= g.minX + g.nx) continue;
                for (int dy = -K; dy <= K; dy++) {
                    int iy = by + dy; if (iy < g.minY || iy >= g.minY + g.ny) continue;
                    for (int dz = -K; dz <= K; dz++) {
                        int lin = dx * S * S + dy * S + dz;
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
                            if (d2 < sumR * sumR) {
                                if (ec + 2 > edges.length) edges = Arrays.copyOf(edges, edges.length * 2);
                                edges[ec++] = i; edges[ec++] = j;
                                degree[i]++; degree[j]++;
                            }
                        }
                    }
                }
            }
        }

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

    @Override public int[] adjacency() { return adj; }
    @Override public int neighborStart(int i) { return adjStart[i]; }
    @Override public int neighborEnd(int i) { return adjStart[i + 1]; }
}
