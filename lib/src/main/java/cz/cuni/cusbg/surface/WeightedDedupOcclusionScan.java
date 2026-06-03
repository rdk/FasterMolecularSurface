package cz.cuni.cusbg.surface;

/**
 * Scalar occlusion scan that emits ONE point per surviving distinct direction, weighted by that
 * direction's tessellation multiplicity, using the shared {@link DirectionMapping}. It is the V2
 * counterpart of {@link GlobalDedupWeightedOcclusionScan} and the scalar fallback for
 * {@link Vectorized256WeightedDedupOcclusionScan} (the path that runs wherever
 * {@code jdk.incubator.vector} is not enabled).
 *
 * <p>Same buried predicate as the reference scan, with a last-occluder-first hint: it computes the
 * per-distinct-direction verdict, then calls {@link DirectionMapping#emitWeighted} to add one weighted
 * point per surviving direction. Paired with a deduplicating store
 * ({@link DistinctFlatSurfacePointStoreV2}) this yields one point per distinct location with a
 * <em>bit-exact area</em> (the weight preserves the multiplicity-based count); paired with a normal
 * store the default {@code addWeighted} re-emits the full multiplicity.
 *
 * <p>The only per-instance state is the {@code buried} verdict scratch, reused across the build's
 * atoms; a scan instance is created per surface, so it is never shared across threads.
 */
final class WeightedDedupOcclusionScan implements OcclusionScan {

    private final int tesslevel;
    private DirectionMapping mapping;   // resolved from the process-wide cache on the first collect()
    private boolean[] buried;           // per-build verdict scratch, reused across atoms

    WeightedDedupOcclusionScan(int tesslevel) {
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

        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
