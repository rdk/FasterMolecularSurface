package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import javax.vecmath.Point3d;
import java.util.List;

/**
 * SIMD occlusion scan using the {@code jdk.incubator.vector} API: the buried test is evaluated for a
 * whole {@link VectorSpecies#length() lane-width} of neighbors at once. Same buried predicate, scan
 * order, and last-occluder-first hint as {@link OcclusionScan#LAST_OCCLUDER_FIRST}; only the inner
 * neighbor loop is vectorized.
 *
 * <p>For a tessellation point {@code (px,py,pz)} the per-neighbor dot product
 * {@code diffX*px + diffY*py + diffZ*pz} is computed with lane-wise {@code mul}/{@code add} (NOT a
 * fused multiply-add), in the same operation order as the scalar scan, so every lane reproduces the
 * scalar result bit-for-bit. A point is buried iff some lane exceeds its threshold ({@code anyTrue}).
 * The result (surviving set, coordinates, areas) is therefore identical to the scalar scans; only the
 * arithmetic throughput differs.
 *
 * <p>This class is referenced only behind a guard (see
 * {@link VectorizedSymmetricHintedGridSoaNumericalSurface}); on a JVM without the incubator module its
 * static initializer fails and the caller falls back to the scalar scan.
 */
final class VectorizedOcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    @Override
    public void collect(double[] tx, double[] ty, double[] tz, int numTess,
                        int numNeighbors, double[] diffX, double[] diffY, double[] diffZ, double[] thresh,
                        double totalRadius, double atomX, double atomY, double atomZ,
                        List<Point3d> points) {
        int laneWidth = SPECIES.length();
        int bound = SPECIES.loopBound(numNeighbors);   // largest multiple of laneWidth <= numNeighbors
        int last = -1;   // last-occluder hint (per-atom local; reset here)

        for (int t = 0; t < numTess; t++) {
            double px = tx[t], py = ty[t], pz = tz[t];

            // hint: test the cached last occluder first (scalar, single comparison)
            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                continue;
            }

            boolean buried = false;
            int k = 0;
            for (; k < bound; k += laneWidth) {
                DoubleVector vx = DoubleVector.fromArray(SPECIES, diffX, k);
                DoubleVector vy = DoubleVector.fromArray(SPECIES, diffY, k);
                DoubleVector vz = DoubleVector.fromArray(SPECIES, diffZ, k);
                DoubleVector vth = DoubleVector.fromArray(SPECIES, thresh, k);
                // lane-wise diffX*px + diffY*py + diffZ*pz, same order as the scalar scan (no FMA)
                DoubleVector dot = vx.mul(px).add(vy.mul(py)).add(vz.mul(pz));
                VectorMask<Double> over = dot.compare(VectorOperators.GT, vth);
                if (over.anyTrue()) {
                    last = k + over.firstTrue();   // first burying neighbor in this lane block
                    buried = true;
                    break;
                }
            }
            if (!buried) {
                // scalar tail for the remaining (numNeighbors % laneWidth) neighbors
                for (; k < numNeighbors; k++) {
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                        last = k;
                        buried = true;
                        break;
                    }
                }
            }
            if (!buried) {
                points.add(new Point3d(totalRadius * px + atomX, totalRadius * py + atomY, totalRadius * pz + atomZ));
            }
        }
    }
}
