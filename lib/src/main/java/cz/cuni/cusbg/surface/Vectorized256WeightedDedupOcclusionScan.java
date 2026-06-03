package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Vectorized (256-bit) counterpart of {@link WeightedDedupOcclusionScan}: same distinct-direction
 * weighted dedup over the shared {@link DirectionMapping}, but the per-direction buried test runs a
 * lane-width of neighbors at a time via the {@code jdk.incubator.vector} API. It is to
 * {@link WeightedDedupOcclusionScan} what {@link GlobalDedupVectorized256OcclusionScan} is to
 * {@link GlobalDedupWeightedOcclusionScan}, and the verdict is bit-for-bit identical to the scalar
 * scan (same per-lane multiply/add order, same {@code >} test, same last-occluder-first hint).
 *
 * <p>Lane width is a {@code static final} 256-bit species - the Vector API needs a constant species to
 * intrinsify, so the width is baked into the class (mirroring
 * {@link GlobalDedupVectorized256OcclusionScan}). This class is referenced only behind a Vector-API
 * availability probe (see {@link DistinctPackedNumericalSurfaceV2}); the scalar
 * {@link WeightedDedupOcclusionScan} is the fallback wherever the incubator module is absent.
 *
 * <p>The only per-instance state is the {@code buried} verdict scratch, reused across the build's
 * atoms; a scan instance is created per surface, so it is never shared across threads.
 */
final class Vectorized256WeightedDedupOcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

    private final int tesslevel;
    private DirectionMapping mapping;   // resolved from the process-wide cache on the first collect()
    private boolean[] buried;           // per-build verdict scratch, reused across atoms

    Vectorized256WeightedDedupOcclusionScan(int tesslevel) {
        this.tesslevel = tesslevel;
    }

    @Override
    public void collect(double[] tx, double[] ty, double[] tz, int numTess,
                        int numNeighbors, double[] diffX, double[] diffY, double[] diffZ, double[] thresh,
                        double totalRadius, double atomX, double atomY, double atomZ,
                        SurfacePointSink sink) {
        DirectionMapping m = mapping;
        if (m == null) {
            m = DirectionMapping.forLevel(tesslevel, tx, ty, tz, numTess);
            mapping = m;
            buried = new boolean[m.numDir];
        }

        double[] ddx = m.ddx, ddy = m.ddy, ddz = m.ddz;
        boolean[] buried = this.buried;
        int numDir = m.numDir;
        int laneWidth = SPECIES.length();
        int bound = SPECIES.loopBound(numNeighbors);
        int last = -1;

        // verdict per distinct direction (a lane-width of neighbors at a time, scalar tail + hint)
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

        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
