package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Backlog idea C3: {@link Vectorized256WeightedDedupOcclusionScan} with the single last-occluder hint
 * generalized to <b>one hint per sphere octant</b>. The instrumentation showed the single hint already
 * resolves ~50% of directions in one test; since a neighbor's occluding cap covers a contiguous patch of
 * the sphere, caching the last occluder <em>per octant</em> of the direction (sign bits of the unit
 * direction) should fire more often as the scan moves across octants. Bit-for-bit identical to the
 * single-hint scan: the hint only changes which neighbor is tested first; the buried verdict is the
 * order-independent OR over all neighbors.
 */
final class RegionHintedVectorized256WeightedDedupOcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

    private final int tesslevel;
    private DirectionMapping mapping;
    private boolean[] buried;
    private int[] region;            // per distinct direction: its octant 0..7 (sign bits of x,y,z)
    private final int[] last = new int[8];   // per-octant last-occluder hint

    RegionHintedVectorized256WeightedDedupOcclusionScan(int tesslevel) {
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
            region = new int[m.numDir];
            for (int d = 0; d < m.numDir; d++) {
                region[d] = (m.ddx[d] < 0 ? 1 : 0) | (m.ddy[d] < 0 ? 2 : 0) | (m.ddz[d] < 0 ? 4 : 0);
            }
        }

        double[] ddx = m.ddx, ddy = m.ddy, ddz = m.ddz;
        boolean[] buried = this.buried;
        int[] region = this.region;
        int[] last = this.last;
        for (int r = 0; r < 8; r++) last[r] = -1;
        int numDir = m.numDir;
        int laneWidth = SPECIES.length();
        int bound = SPECIES.loopBound(numNeighbors);

        for (int d = 0; d < numDir; d++) {
            double px = ddx[d], py = ddy[d], pz = ddz[d];
            int r = region[d];
            int h = last[r];

            if (h >= 0 && diffX[h] * px + diffY[h] * py + diffZ[h] * pz > thresh[h]) {
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
                    last[r] = k + over.firstTrue();
                    b = true;
                    break;
                }
            }
            if (!b) {
                for (; k < numNeighbors; k++) {
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                        last[r] = k;
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
