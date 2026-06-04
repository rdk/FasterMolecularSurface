package cz.cuni.cusbg.surface;

/**
 * A scalar weighted-dedup occlusion scan that also accumulates {@link ScanStats} — used by the
 * feasibility-instrumentation test to measure the headroom of ideas A1/A3/A4 on the real engine
 * (the engine computes diff/thresh exactly as for V2; this scan just observes them). It produces a
 * correct surface (emits via the shared {@link DirectionMapping}) so the build runs normally.
 *
 * <p>Per direction it runs the real scan (last-occluder hint + early exit) to measure A1's early-exit
 * effectiveness, then runs a no-exit pass over all neighbors to mark which neighbors bury at least one
 * direction (A4) and whether the atom is fully buried (A3). The second pass is analysis-only overhead.
 */
final class InstrumentingOcclusionScan implements OcclusionScan {

    private final int tesslevel;
    private final ScanStats stats;
    private DirectionMapping mapping;
    private boolean[] buried;
    private boolean[] neighborBuries;

    InstrumentingOcclusionScan(int tesslevel, ScanStats stats) {
        this.tesslevel = tesslevel;
        this.stats = stats;
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
        int numDir = m.numDir;
        boolean[] buried = this.buried;
        if (neighborBuries == null || neighborBuries.length < numNeighbors) {
            neighborBuries = new boolean[Math.max(numNeighbors, 64)];
        }
        boolean[] nb = neighborBuries;
        for (int k = 0; k < numNeighbors; k++) nb[k] = false;

        int last = -1;
        int surviving = 0;
        long testsThisAtom = 0;

        for (int d = 0; d < numDir; d++) {
            double px = ddx[d], py = ddy[d], pz = ddz[d];

            // (1) real scan: last-occluder hint then early-exit linear scan
            boolean b = false;
            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                b = true;
                stats.hintFires++;
                testsThisAtom++;   // the hint test
            } else {
                for (int k = 0; k < numNeighbors; k++) {
                    testsThisAtom++;
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                        b = true;
                        last = k;
                        break;
                    }
                }
            }
            buried[d] = b;
            stats.directions++;
            if (b) stats.buriedDirections++; else surviving++;

            // (2) analysis pass (no early exit): which neighbors bury this direction
            for (int k = 0; k < numNeighbors; k++) {
                if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) nb[k] = true;
            }
        }

        // per-atom accumulation
        stats.atoms++;
        stats.neighborTests += testsThisAtom;
        stats.noExitTests += (long) numDir * numNeighbors;
        stats.neighborSlots += numNeighbors;
        if (surviving == 0) {
            stats.fullyBuriedAtoms++;
            stats.neighborTestsFullyBuried += testsThisAtom;
        }
        for (int k = 0; k < numNeighbors; k++) if (!nb[k]) stats.zeroBuryNeighbors++;

        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
