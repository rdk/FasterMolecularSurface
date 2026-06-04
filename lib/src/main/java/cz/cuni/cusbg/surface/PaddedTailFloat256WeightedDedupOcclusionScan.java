package cz.cuni.cusbg.surface;

import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Backlog idea C2: {@link Float256WeightedDedupOcclusionScan} with the scalar remainder loop eliminated by
 * padding (idea A7). A7 lost on the 4-wide double scan (rung 22) because the tail it removed was cheap, but
 * the tail is proportionally larger on this 8-wide float scan, so it may pay here. Unlike the double case,
 * the float scan <em>owns</em> its narrowed neighbor scratch ({@code fdx/fdy/fdz/fth}), so it can always
 * size it to a lane-width multiple and pad with non-burying sentinels ({@code thresh=+∞}) — no dependency
 * on the engine's scratch capacity.
 *
 * <p><b>Bit-for-bit identical to {@link Float256WeightedDedupOcclusionScan}</b> (and hence to
 * {@link FloatNumericalSurface}): a sentinel gives {@code dot = 0 > +∞ = false}, so it can never bury a
 * direction nor become the hint.
 */
final class PaddedTailFloat256WeightedDedupOcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;

    private final int tesslevel;
    private DirectionMapping mapping;
    private boolean[] buried;
    private float[] fdx, fdy, fdz, fth;

    PaddedTailFloat256WeightedDedupOcclusionScan(int tesslevel) {
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

        int laneWidth = SPECIES.length();
        int rem = numNeighbors % laneWidth;
        int padded = (rem == 0) ? numNeighbors : numNeighbors + (laneWidth - rem);

        // narrow this atom's neighbor scratch to float, sized to a lane multiple, pad with sentinels
        if (fdx == null || fdx.length < padded) {
            fdx = new float[padded]; fdy = new float[padded]; fdz = new float[padded]; fth = new float[padded];
        }
        float[] fdx = this.fdx, fdy = this.fdy, fdz = this.fdz, fth = this.fth;
        for (int k = 0; k < numNeighbors; k++) {
            fdx[k] = (float) diffX[k]; fdy[k] = (float) diffY[k];
            fdz[k] = (float) diffZ[k]; fth[k] = (float) thresh[k];
        }
        for (int k = numNeighbors; k < padded; k++) {
            fdx[k] = 0f; fdy[k] = 0f; fdz[k] = 0f; fth[k] = Float.POSITIVE_INFINITY;
        }

        float[] fddx = m.fddx, fddy = m.fddy, fddz = m.fddz;
        boolean[] buried = this.buried;
        int numDir = m.numDir;
        int last = -1;

        for (int d = 0; d < numDir; d++) {
            float px = fddx[d], py = fddy[d], pz = fddz[d];

            if (last >= 0 && fdx[last] * px + fdy[last] * py + fdz[last] * pz > fth[last]) {
                buried[d] = true;
                continue;
            }

            boolean b = false;
            for (int k = 0; k < padded; k += laneWidth) {   // padded is a lane multiple -> no scalar tail
                FloatVector vx = FloatVector.fromArray(SPECIES, fdx, k);
                FloatVector vy = FloatVector.fromArray(SPECIES, fdy, k);
                FloatVector vz = FloatVector.fromArray(SPECIES, fdz, k);
                FloatVector vth = FloatVector.fromArray(SPECIES, fth, k);
                FloatVector dot = vx.mul(px).add(vy.mul(py)).add(vz.mul(pz));
                VectorMask<Float> over = dot.compare(VectorOperators.GT, vth);
                if (over.anyTrue()) {
                    last = k + over.firstTrue();
                    b = true;
                    break;
                }
            }
            buried[d] = b;
        }

        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
