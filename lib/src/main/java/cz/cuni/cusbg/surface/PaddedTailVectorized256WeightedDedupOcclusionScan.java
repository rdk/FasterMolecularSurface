package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * {@link Vectorized256WeightedDedupOcclusionScan} with backlog idea A7: the scalar remainder loop is
 * eliminated by padding the neighbor scratch to a multiple of the lane width with <em>non-burying
 * sentinels</em> ({@code diff = 0}, {@code thresh = +∞}), so every distinct direction's neighbor scan is
 * pure-vector with no tail. A sentinel gives {@code dot = 0 > +∞ = false}, so it can never bury a
 * direction nor become the last-occluder hint — the verdict is unchanged, hence <b>bit-for-bit identical
 * to the scalar / un-padded scans</b>.
 *
 * <p>The padding is written once per atom into the engine's reused scratch beyond {@code numNeighbors};
 * the next atom's pre-pass overwrites {@code [0, numNeighbors)} and re-pads, so stale sentinels are never
 * read. It is applied only when the scratch has spare capacity ({@code diffX.length ≥} the padded count);
 * for the one at-capacity atom (the engine grows the scratch to exactly {@code numNeighbors} for the
 * largest neighbor list) it falls back to the scalar tail, so this is never slower than the un-padded scan.
 */
final class PaddedTailVectorized256WeightedDedupOcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

    private final int tesslevel;
    private DirectionMapping mapping;
    private boolean[] buried;

    PaddedTailVectorized256WeightedDedupOcclusionScan(int tesslevel) {
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

        // A7: pad the neighbor scratch with non-burying sentinels so the per-direction scan is all-vector.
        // Done once per atom; falls back to a scalar tail when there is no spare scratch capacity.
        int rem = numNeighbors % laneWidth;
        int scanCount = numNeighbors;
        if (rem != 0) {
            int pad = laneWidth - rem;
            if (diffX.length >= numNeighbors + pad) {
                for (int k = numNeighbors; k < numNeighbors + pad; k++) {
                    diffX[k] = 0.0; diffY[k] = 0.0; diffZ[k] = 0.0; thresh[k] = Double.POSITIVE_INFINITY;
                }
                scanCount = numNeighbors + pad;   // multiple of laneWidth -> no scalar tail
            }
        }
        int bound = SPECIES.loopBound(scanCount);
        int last = -1;

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
                // scalar tail only when padding could not be applied (scanCount == numNeighbors, rem != 0)
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
