package cz.cuni.cusbg.surface;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

import javax.vecmath.Point3d;
import java.util.List;

/**
 * Identical to {@link VectorizedOcclusionScan} but pinned to 256-bit lanes (4 doubles) instead of
 * {@link DoubleVector#SPECIES_PREFERRED}.
 *
 * <p>Why a whole separate class rather than a configurable species: the {@code jdk.incubator.vector}
 * JIT intrinsics require the {@link VectorSpecies} to be a compile-time constant (a {@code static final}
 * field) to generate efficient SIMD; passing the species as an instance field or method parameter
 * defeats intrinsification and collapses throughput by ~10x (measured). So the species must be baked in
 * as a constant per class. This one bakes in {@code SPECIES_256}.
 *
 * <p>Rationale for 256-bit: on AVX-512 hardware {@code SPECIES_PREFERRED} is 512-bit, but for the short
 * per-atom neighbor lists this kernel sees (especially after the per-pair occlusion cutoff in
 * {@link PrunedSymmetricCellGridNeighborList}), the 512-bit bursts are too brief to amortize AVX-512's
 * frequency-license downclock and the scalar<->512-bit transition penalties, which drags the kernel
 * below even scalar speed. 256-bit avoids the downclock while still vectorizing 4-wide.
 *
 * <p>Bit-for-bit identical to the scalar scans (lane-wise {@code mul}/{@code add}, no FMA, same order).
 * Referenced only behind a guard; its static initializer touches the incubator module, so on a JVM
 * without it the caller falls back to the scalar scan.
 */
final class Vectorized256OcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_256;

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
