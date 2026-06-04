package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.openscience.cdk.interfaces.IAtom;

import java.util.Arrays;
import java.util.function.ToDoubleFunction;

/**
 * Backlog idea C1: like {@link SimdDistanceCellGridNeighborList} (A6) but the distance test is done in
 * <b>single precision</b> — candidate coordinates narrowed to {@code float} and tested 8 lanes at a time.
 * A6 showed the neighbor build is bandwidth-bound at 16 threads; halving the candidate-coordinate traffic
 * (float vs double) attacks that ceiling directly.
 *
 * <p><b>Not bit-exact:</b> a pair whose squared distance is within {@code float} epsilon of {@code (Ri+Rj)²}
 * may be classified differently than the double build, so the neighbor set — and hence the surface — can
 * differ slightly. Such boundary pairs have near-empty occluding caps, so the area error is expected to be
 * tiny (tolerance-tested vs the double surface, like {@link FloatNumericalSurface}). The cell grid bucketing
 * stays in {@code double}; only the per-pair distance comparison is narrowed. Implements
 * {@link DirectNeighborSource}.
 */
final class FloatSimdDistanceCellGridNeighborList implements DirectNeighborSource {

    private static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;   // 8 lanes

    private final int[] adjStart;
    private final int[] adj;

    FloatSimdDistanceCellGridNeighborList(IAtom[] atoms, double[] ax, double[] ay, double[] az,
                                          double radius, double solventRadius, ToDoubleFunction<IAtom> vdwRadius) {
        int n = atoms.length;
        CellGrid g = new CellGrid(ax, ay, az, n, radius);

        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = vdwRadius.applyAsDouble(atoms[i]) + solventRadius;

        // cell-sorted candidate coords narrowed to float SoA (half the traffic of the double build)
        float[] csx = new float[n], csy = new float[n], csz = new float[n], csr = new float[n];
        for (int p = 0; p < n; p++) {
            int j = g.cellAtoms[p];
            csx[p] = (float) g.ax[j]; csy[p] = (float) g.ay[j]; csz[p] = (float) g.az[j]; csr[p] = (float) expandedR[j];
        }

        int W = SPECIES.length();
        int[] degree = new int[n];
        int[] edges = new int[Math.max(16, n * 16)];
        int ec = 0;
        for (int i = 0; i < n; i++) {
            float xi = (float) g.ax[i], yi = (float) g.ay[i], zi = (float) g.az[i], ri = (float) expandedR[i];
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
                            for (int p = start; p < end; p++) {
                                int j = g.cellAtoms[p];
                                if (j <= i) continue;
                                float x12 = csx[p] - xi, y12 = csy[p] - yi, z12 = csz[p] - zi;
                                float d2 = x12 * x12 + y12 * y12 + z12 * z12;
                                float sumR = ri + csr[p];
                                if (d2 < sumR * sumR) {
                                    if (ec + 2 > edges.length) edges = Arrays.copyOf(edges, edges.length * 2);
                                    edges[ec++] = i; edges[ec++] = j; degree[i]++; degree[j]++;
                                }
                            }
                            continue;
                        }

                        int p = start;
                        int vEnd = start + SPECIES.loopBound(end - start);
                        for (; p < vEnd; p += W) {
                            FloatVector vx = FloatVector.fromArray(SPECIES, csx, p);
                            FloatVector vy = FloatVector.fromArray(SPECIES, csy, p);
                            FloatVector vz = FloatVector.fromArray(SPECIES, csz, p);
                            FloatVector vr = FloatVector.fromArray(SPECIES, csr, p);
                            FloatVector ddx = vx.sub(xi), ddy = vy.sub(yi), ddz = vz.sub(zi);
                            FloatVector d2 = ddx.mul(ddx).add(ddy.mul(ddy)).add(ddz.mul(ddz));
                            FloatVector sumR = vr.add(ri);
                            VectorMask<Float> hit = d2.compare(VectorOperators.LT, sumR.mul(sumR));
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
                            float x12 = csx[p] - xi, y12 = csy[p] - yi, z12 = csz[p] - zi;
                            float d2 = x12 * x12 + y12 * y12 + z12 * z12;
                            float sumR = ri + csr[p];
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
