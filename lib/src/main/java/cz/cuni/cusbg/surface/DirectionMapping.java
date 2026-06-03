package cz.cuni.cusbg.surface;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Immutable distinct-direction mapping for one tessellation level, shared by the weighted dedup
 * occlusion scans ({@link WeightedDedupOcclusionScan} and
 * {@link Vectorized256WeightedDedupOcclusionScan}).
 *
 * <p>The icosahedral tessellation emits each direction with multiplicity &gt;=5 (shared triangle
 * vertices), so the raw point list holds ~5.7x exact-coincident duplicates. This mapping collapses
 * them to one entry per distinct direction ({@code ddx/ddy/ddz}) plus that direction's multiplicity
 * ({@code mult}); a scan evaluates the buried test once per distinct direction and emits one weighted
 * point, preserving the multiplicity-based area exactly.
 *
 * <p>The mapping is a pure function of the tessellation level (deterministic points in a deterministic
 * order), so it is built once and cached process-wide; the arrays are never mutated, so the shared
 * instance is safe to read from any number of threads concurrently. This is the single place the
 * distinct-direction collapse lives for the V2 stack - the older dedup scans
 * ({@link GlobalDedupWeightedOcclusionScan}, {@link GlobalDedupVectorized256OcclusionScan},
 * {@link DedupVectorized256OcclusionScan}) each carried their own copy of this build, which this
 * class removes for the new scans.
 */
final class DirectionMapping {

    /** Process-wide cache of the immutable mapping, keyed by tessellation level. */
    private static final ConcurrentHashMap<Integer, DirectionMapping> CACHE = new ConcurrentHashMap<>();

    final double[] ddx, ddy, ddz;   // distinct direction unit vectors, length numDir
    final int[] mult;               // mult[d] = number of tessellation points mapping to direction d
    final int numDir;

    private DirectionMapping(double[] ddx, double[] ddy, double[] ddz, int[] mult, int numDir) {
        this.ddx = ddx; this.ddy = ddy; this.ddz = ddz; this.mult = mult; this.numDir = numDir;
    }

    /** Resolve the mapping for {@code tesslevel} from the process-wide cache, building it once on a miss. */
    static DirectionMapping forLevel(int tesslevel, double[] tx, double[] ty, double[] tz, int numTess) {
        DirectionMapping m = CACHE.get(tesslevel);
        if (m == null) {
            m = build(tx, ty, tz, numTess);
            DirectionMapping prev = CACHE.putIfAbsent(tesslevel, m);
            if (prev != null) m = prev;   // lost a benign race; everyone shares one immutable copy
        }
        return m;
    }

    private static DirectionMapping build(double[] tx, double[] ty, double[] tz, int numTess) {
        Map<Dir, Integer> index = new HashMap<>(numTess * 2);
        double[] dx = new double[numTess], dy = new double[numTess], dz = new double[numTess];
        int[] mult = new int[numTess];
        int nd = 0;
        for (int t = 0; t < numTess; t++) {
            Dir key = new Dir(tx[t], ty[t], tz[t]);   // bit-exact duplicates collapse (record equals over doubles)
            Integer idx = index.get(key);
            if (idx == null) {
                idx = nd;
                index.put(key, idx);
                dx[nd] = tx[t]; dy[nd] = ty[t]; dz[nd] = tz[t];
                nd++;
            }
            mult[idx]++;
        }
        return new DirectionMapping(Arrays.copyOf(dx, nd), Arrays.copyOf(dy, nd), Arrays.copyOf(dz, nd),
                Arrays.copyOf(mult, nd), nd);
    }

    /**
     * Emit one weighted surface point per surviving distinct direction (weight = multiplicity, so the
     * area stays multiplicity-exact). Shared by the scalar and vectorized weighted scans, which differ
     * only in how they compute the {@code buried} verdict.
     */
    void emitWeighted(boolean[] buried, double totalRadius, double atomX, double atomY, double atomZ,
                      SurfacePointSink sink) {
        for (int d = 0; d < numDir; d++) {
            if (!buried[d]) {
                sink.addWeighted(totalRadius * ddx[d] + atomX, totalRadius * ddy[d] + atomY,
                        totalRadius * ddz[d] + atomZ, mult[d]);
            }
        }
    }

    /**
     * Map key: a bit-exact 3D direction. Record equality/hashCode is over the raw {@code double}s
     * ({@code Double.equals} semantics), so only exact-coincident tessellation points collapse - the
     * same criterion the older {@code String}-keyed builds used, without the per-point string. Built
     * once per level and cached, so this is off the hot path.
     */
    private record Dir(double x, double y, double z) { }
}
