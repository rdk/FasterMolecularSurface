package cz.cuni.cusbg.surface;

import jdk.incubator.vector.FloatVector;
import jdk.incubator.vector.VectorMask;
import jdk.incubator.vector.VectorOperators;
import jdk.incubator.vector.VectorSpecies;

/**
 * Single-precision (256-bit, 8-lane) weighted dedup occlusion scan - the {@code float} analog of
 * {@link Vectorized256WeightedDedupOcclusionScan}. The buried VERDICT is computed in {@code float}:
 * 8 neighbors per vector instead of 4, and half the memory traffic when the per-atom neighbor scratch
 * is re-read across all distinct directions. Profiling identified this scan as ~72% of
 * {@link DistinctPackedNumericalSurfaceV2}'s CPU at tessellation level 4, and double precision had it
 * bandwidth/throughput bound - this is the direct lever on that ceiling.
 *
 * <p>Only the verdict is single precision. The emitted point POSITIONS and the area-defining weights
 * stay in {@code double}: {@link DirectionMapping#emitWeighted} uses the double distinct directions,
 * the double atom centre and the double radius, so a surviving direction's coordinates are bit-identical
 * to the double scans. Likewise the coordinate subtraction stays in {@code double} (done by the engine
 * pre-pass) to avoid cancellation; this scan only narrows the already-differenced {@code diff/thresh}
 * scratch to {@code float} for the comparison.
 *
 * <p><b>Not bit-exact:</b> a tessellation point within {@code float} epsilon of the occlusion boundary
 * may flip survival, so the surviving SET (and hence the area) can differ slightly from the double
 * scans - within SAS tessellation discretization error. It therefore ships behind a tolerance-tested
 * surface ({@link FloatNumericalSurface}), not the bit-exact equivalence harness.
 *
 * <p>Per-instance state only (the float neighbor scratch and the {@code buried} verdict array); a scan
 * instance is created per surface, so it is never shared across threads.
 */
final class Float256WeightedDedupOcclusionScan implements OcclusionScan {

    private static final VectorSpecies<Float> SPECIES = FloatVector.SPECIES_256;

    private final int tesslevel;
    private DirectionMapping mapping;   // resolved from the process-wide cache on the first collect()
    private boolean[] buried;           // per-build verdict scratch, reused across atoms
    private float[] fdx, fdy, fdz, fth; // per-atom neighbor scratch narrowed to float, reused across atoms

    Float256WeightedDedupOcclusionScan(int tesslevel) {
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

        // narrow this atom's neighbor scratch to float once; reused across all distinct directions below
        if (fdx == null || fdx.length < numNeighbors) {
            fdx = new float[numNeighbors]; fdy = new float[numNeighbors];
            fdz = new float[numNeighbors]; fth = new float[numNeighbors];
        }
        float[] fdx = this.fdx, fdy = this.fdy, fdz = this.fdz, fth = this.fth;
        for (int k = 0; k < numNeighbors; k++) {
            fdx[k] = (float) diffX[k]; fdy[k] = (float) diffY[k];
            fdz[k] = (float) diffZ[k]; fth[k] = (float) thresh[k];
        }

        float[] fddx = m.fddx, fddy = m.fddy, fddz = m.fddz;
        boolean[] buried = this.buried;
        int numDir = m.numDir;
        int laneWidth = SPECIES.length();
        int bound = SPECIES.loopBound(numNeighbors);
        int last = -1;

        // verdict per distinct direction (8 float neighbors at a time, scalar tail + last-occluder hint)
        for (int d = 0; d < numDir; d++) {
            float px = fddx[d], py = fddy[d], pz = fddz[d];

            if (last >= 0 && fdx[last] * px + fdy[last] * py + fdz[last] * pz > fth[last]) {
                buried[d] = true;
                continue;
            }

            boolean b = false;
            int k = 0;
            for (; k < bound; k += laneWidth) {
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
            if (!b) {
                for (; k < numNeighbors; k++) {
                    if (fdx[k] * px + fdy[k] * py + fdz[k] * pz > fth[k]) {
                        last = k;
                        b = true;
                        break;
                    }
                }
            }
            buried[d] = b;
        }

        // emit surviving directions in DOUBLE (positions + area weights bit-identical to the double scans)
        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
