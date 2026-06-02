package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Direction-deduplicated 256-bit SIMD occlusion scan.
 *
 * <p>CDK's icosahedral tessellation emits each surface direction multiple times (≈5.7x at every level:
 * e.g. 240 points but only 42 distinct directions at tess 2), because it returns the three vertices of
 * every triangle and vertices are shared. The buried test {@code diff.p > thresh} depends only on the
 * direction {@code p}, so all duplicates of a direction get the same verdict. This scan therefore
 * evaluates each <em>distinct</em> direction once (≈5.7x fewer occlusion scans) and then re-expands:
 * it walks the original tessellation points in order and emits each one whose direction survived. The
 * emitted multiset and its order are identical to scanning every point, so the result (surface points,
 * their order, and per-atom areas) is bit-for-bit unchanged.
 *
 * <p>The distinct-direction mapping is a pure function of the tessellation arrays, which the engine
 * reuses (same array instances) for every atom of a build, so it is computed once per build and
 * memoized on this instance. The instance therefore carries per-build state and must NOT be shared
 * across builds/threads; the owning variant creates a fresh one per surface (see
 * {@link DedupVectorizedSymmetricHintedGridSoaNumericalSurface}).
 *
 * <p>The lane width is a {@code static final} 256-bit species (the Vector API requires a constant
 * species for JIT intrinsification; see {@link Vectorized256OcclusionScan}). Bit-for-bit identical to
 * the scalar scans (lane-wise {@code mul}/{@code add}, no FMA, same order).
 */
final class DedupVectorized256OcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

    // per-build memo of the tessellation's distinct directions
    private double[] mapTx;            // identity of the tx array the mapping was built for
    private double[] ddx, ddy, ddz;    // distinct directions
    private int[] dirOf;               // dirOf[t] = distinct-direction index of original point t
    private int numDir;
    private boolean[] buried;          // scratch verdict per distinct direction

    @Override
    public void collect(double[] tx, double[] ty, double[] tz, int numTess,
                        int numNeighbors, double[] diffX, double[] diffY, double[] diffZ, double[] thresh,
                        double totalRadius, double atomX, double atomY, double atomZ,
                        SurfacePointSink sink) {
        if (mapTx != tx) buildMapping(tx, ty, tz, numTess);   // once per build (same tx for every atom)

        double[] ddx = this.ddx, ddy = this.ddy, ddz = this.ddz;
        boolean[] buried = this.buried;
        int numDir = this.numDir;
        int laneWidth = SPECIES.length();
        int bound = SPECIES.loopBound(numNeighbors);
        int last = -1;   // last-occluder hint across distinct directions

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

        // 2) re-expand: emit surviving points in original tessellation order (preserves multiset + order)
        int[] dirOf = this.dirOf;
        for (int t = 0; t < numTess; t++) {
            if (!buried[dirOf[t]]) {
                sink.add(totalRadius * tx[t] + atomX, totalRadius * ty[t] + atomY, totalRadius * tz[t] + atomZ);
            }
        }
    }

    /** Build the distinct-direction set and the original-point -> direction map (once per build). */
    private void buildMapping(double[] tx, double[] ty, double[] tz, int numTess) {
        Map<String, Integer> index = new HashMap<>(numTess * 2);
        int[] dir = new int[numTess];
        double[] dx = new double[numTess], dy = new double[numTess], dz = new double[numTess];
        int nd = 0;
        for (int t = 0; t < numTess; t++) {
            String key = tx[t] + "," + ty[t] + "," + tz[t];   // exact-value key: duplicates have identical coords
            Integer idx = index.get(key);
            if (idx == null) {
                idx = nd;
                index.put(key, idx);
                dx[nd] = tx[t]; dy[nd] = ty[t]; dz[nd] = tz[t];
                nd++;
            }
            dir[t] = idx;
        }
        this.numDir = nd;
        this.ddx = Arrays.copyOf(dx, nd);
        this.ddy = Arrays.copyOf(dy, nd);
        this.ddz = Arrays.copyOf(dz, nd);
        this.dirOf = dir;
        this.buried = new boolean[nd];
        this.mapTx = tx;
    }
}
