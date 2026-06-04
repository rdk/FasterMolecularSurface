package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.openscience.cdk.interfaces.IAtom;

import java.util.Arrays;
import java.util.function.ToDoubleFunction;

/**
 * Backlog idea A6: the same neighbor set and CSR output as
 * {@link CoordSortedPrunedSymmetricCellGridNeighborList} (same {@link CellGrid}, same ±1 stencil, same
 * per-pair {@code d2 < (Ri+Rj)²} prune), but the distance pass — the build that profiling puts at ~47% of
 * CPU at 16 threads / tess 2 — is <b>SIMD-vectorized</b>. The cell-sorted candidate coordinates are kept
 * as four SoA streams ({@code csx/csy/csz/csr}) so a lane-width of candidates is tested per iteration with
 * a masked compare; surviving lanes are drained scalar (low survivor density). The same-cell block is
 * handled scalar (its {@code j>i} guard does not vectorize cleanly), so the vectorized path covers the
 * cross-cell candidates, which dominate.
 *
 * <p><b>Bit-for-bit identical neighbor set</b> to the scalar build: the lane-wise
 * {@code dx·dx + dy·dy + dz·dz} (no FMA, same order as scalar) and the {@code < sumR²} compare reproduce
 * the scalar test exactly per lane, and the grid/stencil/prune are unchanged — so the surviving pairs are
 * the same. (Edge recording order differs, which does not affect areas or the surface point set: the scan
 * is an order-independent OR over neighbors and emission is per-atom/per-direction.) Implements
 * {@link DirectNeighborSource} for copy-free access. Same dense-grid precondition as {@link CellGrid}.
 *
 * <p>Uses {@code jdk.incubator.vector}; reference it only behind a Vector-API availability guard (the
 * surface falls back to the scalar build otherwise).
 */
final class SimdDistanceCellGridNeighborList implements DirectNeighborSource {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

    private final int[] adjStart;
    private final int[] adj;

    SimdDistanceCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                     double radius, double solventRadius) {
        this(atoms, ax, ay, az, radius, solventRadius, FasterNumericalSurface::getVdwRadius);
    }

    SimdDistanceCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                     double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        int n = atoms.length;
        CellGrid g = new CellGrid(ax, ay, az, n, radius);

        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // cell-sorted candidate coordinates as SoA streams (vs the AoS coords[4p] of the scalar build),
        // so a lane of candidates loads contiguously for the masked distance compare.
        double[] csx = new double[n], csy = new double[n], csz = new double[n], csr = new double[n];
        for (int p = 0; p < n; p++) {
            int j = g.cellAtoms[p];
            csx[p] = g.ax[j]; csy[p] = g.ay[j]; csz[p] = g.az[j]; csr[p] = expandedR[j];
        }

        int W = SPECIES.length();
        int[] degree = new int[n];
        int[] edges = new int[Math.max(16, n * 16)];
        int ec = 0;
        for (int i = 0; i < n; i++) {
            double xi = g.ax[i], yi = g.ay[i], zi = g.az[i], ri = expandedR[i];
            int bx = g.cx[i], by = g.cy[i], bz = g.cz[i];
            for (int dx = -1; dx <= 1; dx++) {
                int ix = bx + dx; if (ix < g.minX || ix >= g.minX + g.nx) continue;
                for (int dy = -1; dy <= 1; dy++) {
                    int iy = by + dy; if (iy < g.minY || iy >= g.minY + g.ny) continue;
                    for (int dz = -1; dz <= 1; dz++) {
                        int lin = dx * 9 + dy * 3 + dz;
                        if (lin < 0) continue;
                        int iz = bz + dz; if (iz < g.minZ || iz >= g.minZ + g.nz) continue;
                        int c = g.cellIndex(ix, iy, iz);
                        int start = g.cellStart[c], end = g.cellStart[c + 1];

                        if (lin == 0) {
                            // same cell: scalar, each pair once (j > i)
                            for (int p = start; p < end; p++) {
                                int j = g.cellAtoms[p];
                                if (j <= i) continue;
                                double x12 = csx[p] - xi, y12 = csy[p] - yi, z12 = csz[p] - zi;
                                double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                                double sumR = ri + csr[p];
                                if (d2 < sumR * sumR) {
                                    if (ec + 2 > edges.length) edges = Arrays.copyOf(edges, edges.length * 2);
                                    edges[ec++] = i; edges[ec++] = j; degree[i]++; degree[j]++;
                                }
                            }
                            continue;
                        }

                        // cross-cell: vectorized distance test over the candidate run, scalar tail
                        int p = start;
                        int vEnd = start + SPECIES.loopBound(end - start);
                        for (; p < vEnd; p += W) {
                            DoubleVector vx = DoubleVector.fromArray(SPECIES, csx, p);
                            DoubleVector vy = DoubleVector.fromArray(SPECIES, csy, p);
                            DoubleVector vz = DoubleVector.fromArray(SPECIES, csz, p);
                            DoubleVector vr = DoubleVector.fromArray(SPECIES, csr, p);
                            DoubleVector ddx = vx.sub(xi), ddy = vy.sub(yi), ddz = vz.sub(zi);
                            // lane-wise dx*dx + dy*dy + dz*dz, same order as scalar (no FMA)
                            DoubleVector d2 = ddx.mul(ddx).add(ddy.mul(ddy)).add(ddz.mul(ddz));
                            DoubleVector sumR = vr.add(ri);
                            VectorMask<Double> hit = d2.compare(VectorOperators.LT, sumR.mul(sumR));
                            long bits = hit.toLong();
                            while (bits != 0) {
                                int lane = Long.numberOfTrailingZeros(bits);
                                bits &= bits - 1;
                                int j = g.cellAtoms[p + lane];
                                if (ec + 2 > edges.length) edges = Arrays.copyOf(edges, edges.length * 2);
                                edges[ec++] = i; edges[ec++] = j; degree[i]++; degree[j]++;
                            }
                        }
                        for (; p < end; p++) {
                            int j = g.cellAtoms[p];
                            double x12 = csx[p] - xi, y12 = csy[p] - yi, z12 = csz[p] - zi;
                            double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            double sumR = ri + csr[p];
                            if (d2 < sumR * sumR) {
                                if (ec + 2 > edges.length) edges = Arrays.copyOf(edges, edges.length * 2);
                                edges[ec++] = i; edges[ec++] = j; degree[i]++; degree[j]++;
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
