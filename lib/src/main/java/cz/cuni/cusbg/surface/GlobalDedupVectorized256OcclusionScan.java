package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Like {@link DedupVectorized256OcclusionScan}, but the distinct-direction mapping is cached
 * <em>process-wide</em> (keyed by tessellation level) instead of rebuilt for every surface.
 *
 * <p>The mapping is a pure function of the tessellation, which is deterministic for a given level
 * ({@code "ico"}, identical points in identical order every build). Profiling the per-build dedup at
 * 16 threads showed its {@code String}-keyed {@code HashMap} rebuild was the largest allocation source
 * ({@code byte[]}/{@code String} churn), and the per-build construction had become a scaling tax once
 * the dedup made the scan cheap. Building it once and sharing an immutable {@link Mapping} removes that
 * churn entirely: every subsequent build at the same level is a cache hit. The mapping arrays are
 * never mutated, so the shared instance is safe to read from any number of threads concurrently.
 *
 * <p>The only per-build (per-instance) state is the {@code buried} verdict scratch, reused across the
 * build's atoms; this scan instance is created per surface, so it is never shared across threads.
 * Lane width is a {@code static final} 256-bit species (the Vector API needs a constant species to
 * intrinsify). Bit-for-bit identical to the scalar scans.
 */
final class GlobalDedupVectorized256OcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

    /** Process-wide cache of the immutable direction mapping, keyed by tessellation level. */
    private static final ConcurrentHashMap<Integer, Mapping> CACHE = new ConcurrentHashMap<>();

    private final int tesslevel;
    private Mapping mapping;     // resolved from the cache on the first collect() of this build
    private boolean[] buried;    // per-build verdict scratch, reused across atoms

    GlobalDedupVectorized256OcclusionScan(int tesslevel) {
        this.tesslevel = tesslevel;
    }

    @Override
    public void collect(double[] tx, double[] ty, double[] tz, int numTess,
                        int numNeighbors, double[] diffX, double[] diffY, double[] diffZ, double[] thresh,
                        double totalRadius, double atomX, double atomY, double atomZ,
                        SurfacePointSink sink) {
        Mapping m = this.mapping;
        if (m == null) {
            m = CACHE.get(tesslevel);
            if (m == null) {
                m = Mapping.build(tx, ty, tz, numTess);
                Mapping prev = CACHE.putIfAbsent(tesslevel, m);
                if (prev != null) m = prev;   // lost a benign race; everyone shares one immutable copy
            }
            this.mapping = m;
            this.buried = new boolean[m.numDir];
        }

        double[] ddx = m.ddx, ddy = m.ddy, ddz = m.ddz;
        int[] dirOf = m.dirOf;
        boolean[] buried = this.buried;
        int numDir = m.numDir;
        int laneWidth = SPECIES.length();
        int bound = SPECIES.loopBound(numNeighbors);
        int last = -1;

        // 1) verdict per distinct direction
        for (int d = 0; d < numDir; d++) {
            double px = ddx[d], py = ddy[d], pz = ddz[d];

            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                buried[d] = true;
                continue;
            }

            boolean b = false;
            int k = 0;
            for (; k < bound; k += laneWidth) {
                DoubleVector vx = DoubleVector.fromArray(SPECIES, diffX, k);
                DoubleVector vy = DoubleVector.fromArray(SPECIES, diffY, k);
                DoubleVector vz = DoubleVector.fromArray(SPECIES, diffZ, k);
                DoubleVector vth = DoubleVector.fromArray(SPECIES, thresh, k);
                DoubleVector dot = vx.mul(px).add(vy.mul(py)).add(vz.mul(pz));
                VectorMask<Double> over = dot.compare(VectorOperators.GT, vth);
                if (over.anyTrue()) {
                    last = k + over.firstTrue();
                    b = true;
                    break;
                }
            }
            if (!b) {
                for (; k < numNeighbors; k++) {
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                        last = k;
                        b = true;
                        break;
                    }
                }
            }
            buried[d] = b;
        }

        // 2) re-expand in original tessellation order (preserves multiset + order)
        for (int t = 0; t < numTess; t++) {
            if (!buried[dirOf[t]]) {
                sink.add(totalRadius * tx[t] + atomX, totalRadius * ty[t] + atomY, totalRadius * tz[t] + atomZ);
            }
        }
    }

    /** Immutable distinct-direction mapping for one tessellation level; built once, shared read-only. */
    static final class Mapping {
        final double[] ddx, ddy, ddz;
        final int[] dirOf;
        final int numDir;

        private Mapping(double[] ddx, double[] ddy, double[] ddz, int[] dirOf, int numDir) {
            this.ddx = ddx; this.ddy = ddy; this.ddz = ddz; this.dirOf = dirOf; this.numDir = numDir;
        }

        static Mapping build(double[] tx, double[] ty, double[] tz, int numTess) {
            Map<String, Integer> index = new HashMap<>(numTess * 2);
            int[] dirOf = new int[numTess];
            double[] dx = new double[numTess], dy = new double[numTess], dz = new double[numTess];
            int nd = 0;
            for (int t = 0; t < numTess; t++) {
                String key = tx[t] + "," + ty[t] + "," + tz[t];   // exact-value key (duplicates are bit-equal)
                Integer idx = index.get(key);
                if (idx == null) {
                    idx = nd;
                    index.put(key, idx);
                    dx[nd] = tx[t]; dy[nd] = ty[t]; dz[nd] = tz[t];
                    nd++;
                }
                dirOf[t] = idx;
            }
            return new Mapping(Arrays.copyOf(dx, nd), Arrays.copyOf(dy, nd), Arrays.copyOf(dz, nd), dirOf, nd);
        }
    }
}
