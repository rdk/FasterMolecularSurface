package cz.cuni.cusbg.surface;

/**
 * Strategy for the occlusion scan in {@link DevSurfaceV1Soa}: given one atom's tessellation
 * points and its per-neighbor scratch (the parallel {@code diffX/diffY/diffZ/thresh} arrays, first
 * {@code numNeighbors} entries), decide which tessellation points survive and append the surviving
 * surface points to {@code sink}.
 *
 * <p>A point at unit-sphere direction {@code (px,py,pz)} is <em>buried</em> iff some neighbor
 * satisfies {@code diffX*px + diffY*py + diffZ*pz > thresh}; surviving points are mapped onto the
 * atom's expanded sphere as {@code (totalRadius*p + atomCenter)}. Every strategy must implement
 * exactly that predicate, so the surviving set, the point coordinates, and the areas are identical
 * regardless of strategy — only <em>how</em> the scan reaches the verdict (scan order, early-exit
 * hints) may differ. See {@link #STANDARD}.
 *
 * <p>Like {@link NeighborOrdering}, the strategy is supplied to {@link DevSurfaceV1Soa} through
 * its constructor (not an overridable method), so it is fixed before any computation runs and never
 * dispatches on a partially-constructed subclass. Implementations must therefore be stateless across
 * atoms: any per-atom state (e.g. a last-occluder hint) must be a local of {@link #collect}, reset on
 * each call.
 */
@FunctionalInterface
interface OcclusionScan {

    /**
     * Append the surviving surface points for one atom to {@code sink}. Reads only entries
     * {@code [0, numNeighbors)} of the four parallel arrays.
     */
    void collect(double[] tx, double[] ty, double[] tz, int numTess,
                 int numNeighbors, double[] diffX, double[] diffY, double[] diffZ, double[] thresh,
                 double totalRadius, double atomX, double atomY, double atomZ,
                 SurfacePointSink sink);

    /** Reference scan: test every neighbor in array order, break on the first that buries the point. */
    OcclusionScan STANDARD = (tx, ty, tz, numTess, numNeighbors, diffX, diffY, diffZ, thresh,
                              totalRadius, atomX, atomY, atomZ, sink) -> {
        for (int t = 0; t < numTess; t++) {
            double px = tx[t], py = ty[t], pz = tz[t];
            boolean buried = false;
            for (int k = 0; k < numNeighbors; k++) {
                if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                    buried = true;
                    break;
                }
            }
            if (!buried) {
                sink.add(totalRadius * px + atomX, totalRadius * py + atomY, totalRadius * pz + atomZ);
            }
        }
    };

    /**
     * Same buried predicate as {@link #STANDARD}, with a last-occluder-first hint: the neighbor that
     * buried the previous tessellation point is tested first for the next one. Consecutive tessellation
     * points are spatially coherent, so that neighbor usually buries the next point too, exiting after a
     * single comparison; on a miss it falls back to the full in-order scan (which still tests the hinted
     * neighbor, so the result is unchanged). The hint is one {@code int} local per atom, reset on entry,
     * so the scan is stateless across atoms. See {@link DevSurfaceV4Hinted}.
     */
    OcclusionScan LAST_OCCLUDER_FIRST = (tx, ty, tz, numTess, numNeighbors, diffX, diffY, diffZ, thresh,
                                         totalRadius, atomX, atomY, atomZ, sink) -> {
        int last = -1;   // neighbor that buried the previous point, or -1 (per-atom local; reset here)
        for (int t = 0; t < numTess; t++) {
            double px = tx[t], py = ty[t], pz = tz[t];

            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                continue;   // buried by the cached occluder, last unchanged
            }

            boolean buried = false;
            for (int k = 0; k < numNeighbors; k++) {
                if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                    last = k;   // remember this occluder for the next point
                    buried = true;
                    break;
                }
            }
            if (!buried) {
                sink.add(totalRadius * px + atomX, totalRadius * py + atomY, totalRadius * pz + atomZ);
            }
        }
    };
}
