package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Arrays;

/**
 * Phase 5 (Lead: re-profile the build) instrumentation. A phase-decomposed CLONE of the V3 build
 * ({@link SimdDistanceCellGridNeighborList}) — the frozen class is NOT touched — with {@code nanoTime}
 * around each phase, run over the corpus, reporting where build time goes now that A6 has vectorized the
 * distance pass. Single-thread (the 16-thread bandwidth bottleneck is a separate effect; this localizes
 * the compute target). See {@code autoresearch/LOG.md} Phase 5.
 *
 * <p>Run: {@code ./gradlew scorecard --tests '*BuildProfileTest'}
 */
@Tag("scorecard")
class BuildProfileTest {

    private static final double SOLVENT = 1.4;
    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;
    private static final String[] PHASES = {"grid", "expandedR", "SoA repack", "distance pass", "edges->CSR"};

    @Test
    void profile() {
        // load the corpus once; the bit-exact engine uses radius = maxVdw + solvent for the grid cutoff
        IAtomContainer[] mols = TestStructures.loadAll();

        long[] acc = new long[5];
        long[] warm = new long[5];
        int totalAtoms = 0, totalEdges = 0;
        long[] candCount = new long[1];   // candidate distance tests performed (set on the first measured pass)

        // warmup (discard) then measured passes
        for (int pass = 0; pass < 8; pass++) {
            boolean measure = pass >= 3;
            long[] target = measure ? acc : warm;
            for (IAtomContainer mol : mols) {
                int[] counts = profileBuild(mol, target, (measure && pass == 3) ? candCount : null);
                if (measure && pass == 3) { totalAtoms += counts[0]; totalEdges += counts[1]; }
            }
        }
        int measuredPasses = 5;

        long sum = 0; for (long v : acc) sum += v;
        System.out.printf("%n=== Phase 5 build profile (V3 SIMD build clone, %d structures, %d measured passes, single-thread) ===%n",
                mols.length, measuredPasses);
        System.out.printf("corpus atoms=%d  undirected pairs=%d  total build %.2f ms/pass%n",
                totalAtoms, totalEdges, sum / 1e6 / measuredPasses);
        for (int i = 0; i < 5; i++) {
            System.out.printf("  %-14s %5.1f%%   (%.3f ms/pass)%n",
                    PHASES[i], sum == 0 ? 0 : 100.0 * acc[i] / sum, acc[i] / 1e6 / measuredPasses);
        }
        System.out.printf("  distance-pass candidate tests=%d for %d surviving pairs -> %.2f candidates/survivor (%.0f%% survive)%n",
                candCount[0], totalEdges,
                totalEdges == 0 ? 0 : (double) candCount[0] / totalEdges,
                candCount[0] == 0 ? 0 : 100.0 * totalEdges / candCount[0]);
        System.out.println();
    }

    /** Replays the V3 build with per-phase timers accumulated into {@code t}. Returns {n, undirectedPairs}. */
    private static int[] profileBuild(IAtomContainer mol, long[] t, long[] cand) {
        long candTests = 0;
        // mirror the engine's coordinate + radius extraction (see DevSurfaceV1Soa / FasterNumericalSurface)
        int n = mol.getAtomCount();
        IAtom[] atoms = new IAtom[n];
        double[] ax = new double[n], ay = new double[n], az = new double[n];
        double maxVdw = 0;
        for (int i = 0; i < n; i++) {
            IAtom a = mol.getAtom(i);
            atoms[i] = a;
            ax[i] = a.getPoint3d().x; ay[i] = a.getPoint3d().y; az[i] = a.getPoint3d().z;
            maxVdw = Math.max(maxVdw, VdwRadiusCache.get(a));
        }
        double radius = maxVdw + SOLVENT;

        long s;
        // (A) grid
        s = System.nanoTime();
        CellGrid g = new CellGrid(ax, ay, az, n, radius);
        t[0] += System.nanoTime() - s;

        // (B) expandedR
        s = System.nanoTime();
        double[] expandedR = new double[n];
        for (int i = 0; i < n; i++) expandedR[i] = VdwRadiusCache.get(atoms[i]) + SOLVENT;
        t[1] += System.nanoTime() - s;

        // (C) SoA repack
        s = System.nanoTime();
        double[] csx = new double[n], csy = new double[n], csz = new double[n], csr = new double[n];
        for (int p = 0; p < n; p++) {
            int j = g.cellAtoms[p];
            csx[p] = g.ax[j]; csy[p] = g.ay[j]; csz[p] = g.az[j]; csr[p] = expandedR[j];
        }
        t[2] += System.nanoTime() - s;

        // (D) distance pass
        s = System.nanoTime();
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
                            for (int p = start; p < end; p++) {
                                int j = g.cellAtoms[p];
                                if (j <= i) continue;
                                candTests++;
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
                        int p = start;
                        int vEnd = start + SPECIES.loopBound(end - start);
                        candTests += (end - start);   // every cross-cell candidate gets a distance test (vector lanes + scalar tail)
                        for (; p < vEnd; p += W) {
                            DoubleVector vx = DoubleVector.fromArray(SPECIES, csx, p);
                            DoubleVector vy = DoubleVector.fromArray(SPECIES, csy, p);
                            DoubleVector vz = DoubleVector.fromArray(SPECIES, csz, p);
                            DoubleVector vr = DoubleVector.fromArray(SPECIES, csr, p);
                            DoubleVector ddx = vx.sub(xi), ddy = vy.sub(yi), ddz = vz.sub(zi);
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
        t[3] += System.nanoTime() - s;

        // (E) edges -> CSR
        s = System.nanoTime();
        int[] adjStart = new int[n + 1];
        for (int i = 0; i < n; i++) adjStart[i + 1] = adjStart[i] + degree[i];
        int[] adj = new int[adjStart[n]];
        int[] cursor = new int[n];
        System.arraycopy(adjStart, 0, cursor, 0, n);
        for (int e = 0; e < ec; e += 2) {
            int i = edges[e], j = edges[e + 1];
            adj[cursor[i]++] = j;
            adj[cursor[j]++] = i;
        }
        t[4] += System.nanoTime() - s;

        if (cand != null) cand[0] += candTests;
        return new int[]{n, ec / 2};
    }
}
