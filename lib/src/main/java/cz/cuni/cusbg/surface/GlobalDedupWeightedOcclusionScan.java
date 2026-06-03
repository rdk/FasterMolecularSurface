package cz.cuni.cusbg.surface;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Occlusion scan that emits ONE point per surviving distinct direction (weighted by its tessellation
 * multiplicity), instead of re-expanding to the full multiplicity like
 * {@link GlobalDedupVectorized256OcclusionScan}.
 *
 * <p>The icosahedral tessellation emits each direction with multiplicity >=5 (shared triangle
 * vertices), so the standard scans produce ~5.7x exact-coincident duplicate points that p2rank then
 * sparsifies away. This scan computes the same per-distinct-direction buried verdict, then for each
 * surviving direction calls {@link SurfacePointSink#addWeighted} once with that direction's
 * multiplicity. Paired with a deduplicating store ({@link DistinctFlatSurfacePointStore}) this yields a
 * surface with one point per distinct location and a <em>bit-exact area</em> (the weight preserves the
 * multiplicity-based count). Paired with a normal store the default {@code addWeighted} re-emits the
 * full multiplicity, reproducing the standard output exactly (used to validate this scan).
 *
 * <p>Scalar (no Vector API dependency): the dedup itself is the algorithmic win, and this is the path
 * that runs wherever {@code jdk.incubator.vector} is not enabled. The distinct-direction mapping
 * (coordinates + multiplicity) is a pure function of the tessellation level, cached process-wide.
 * The only per-instance state is the {@code buried} scratch; a scan instance is created per surface.
 */
final class GlobalDedupWeightedOcclusionScan implements OcclusionScan {

    /** Process-wide cache of the immutable distinct-direction mapping, keyed by tessellation level. */
    private static final ConcurrentHashMap<Integer, Mapping> CACHE = new ConcurrentHashMap<>();

    private final int tesslevel;
    private Mapping mapping;
    private boolean[] buried;

    GlobalDedupWeightedOcclusionScan(int tesslevel) {
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
                if (prev != null) m = prev;
            }
            this.mapping = m;
            this.buried = new boolean[m.numDir];
        }

        double[] ddx = m.ddx, ddy = m.ddy, ddz = m.ddz;
        int[] mult = m.mult;
        boolean[] buried = this.buried;
        int numDir = m.numDir;
        int last = -1;

        // verdict per distinct direction (scalar, with a last-occluder-first hint)
        for (int d = 0; d < numDir; d++) {
            double px = ddx[d], py = ddy[d], pz = ddz[d];
            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                buried[d] = true;
                continue;
            }
            boolean b = false;
            for (int k = 0; k < numNeighbors; k++) {
                if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                    last = k;
                    b = true;
                    break;
                }
            }
            buried[d] = b;
        }

        // emit one weighted point per surviving distinct direction (weight = multiplicity, for the area)
        for (int d = 0; d < numDir; d++) {
            if (!buried[d]) {
                sink.addWeighted(totalRadius * ddx[d] + atomX, totalRadius * ddy[d] + atomY, totalRadius * ddz[d] + atomZ, mult[d]);
            }
        }
    }

    /** Immutable distinct-direction mapping (coordinates + per-direction multiplicity) for one level. */
    static final class Mapping {
        final double[] ddx, ddy, ddz;
        final int[] mult;     // mult[d] = number of tessellation points mapping to distinct direction d
        final int numDir;

        private Mapping(double[] ddx, double[] ddy, double[] ddz, int[] mult, int numDir) {
            this.ddx = ddx; this.ddy = ddy; this.ddz = ddz; this.mult = mult; this.numDir = numDir;
        }

        static Mapping build(double[] tx, double[] ty, double[] tz, int numTess) {
            Map<String, Integer> index = new HashMap<>(numTess * 2);
            double[] dx = new double[numTess], dy = new double[numTess], dz = new double[numTess];
            int[] mult = new int[numTess];
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
                mult[idx]++;
            }
            return new Mapping(Arrays.copyOf(dx, nd), Arrays.copyOf(dy, nd), Arrays.copyOf(dz, nd), Arrays.copyOf(mult, nd), nd);
        }
    }
}
