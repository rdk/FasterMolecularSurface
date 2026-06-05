package cz.cuni.cusbg.surface;

/**
 * Phase 3 (Lead A5 / DCLM) feasibility instrumentation. A5 is the same neighbor-major spatial-narrowing
 * family Phase 1 ({@link BitmaskFeasibilityScan}) closed at tess 3; its one distinct, untested angle is
 * the cost of <em>surviving</em> directions — those the current direction-major scan must test against
 * ALL neighbors because no occluder lets it early-exit. DCLM's spatial bin could test a survivor against
 * only the neighbors whose cap is near it. This scan measures, over the corpus, how big that survivor
 * cost is and how much a (perfectly-binned) DCLM could shrink it.
 *
 * <p>Reuses {@link CapDirectionLut} as the spatial index (same role a DCLM lattice plays): a survivor's
 * "candidate neighbors" are the neighbors whose cap-LUT mask contains that direction — i.e. the neighbors
 * a sound DCLM bin would still have to exact-test. Prints only.
 */
final class DclmFeasibilityScan implements OcclusionScan {

    static final class Stats {
        int lutGrid, lutBins;
        long atoms;
        long directions;                 // total distinct directions scanned
        long survivingDirections;        // directions that survive (no neighbor buries them)

        long currentTests;               // direction-major hint+early-exit (the cost baseline)
        long survivorCurrentTests;       // current tests spent on surviving directions (= sum numNeighbors over survivors)
        long survivorDclmCandidates;     // DCLM tests on survivors (= neighbors whose cap-LUT includes the survivor)

        void report(String title) {
            System.out.printf("%n=== %s (LUT grid=%d bins=%d) ===%n", title, lutGrid, lutBins);
            System.out.printf("atoms=%d  directions=%d  surviving=%d (%.1f%%)%n",
                    atoms, directions, survivingDirections, pct(survivingDirections, directions));
            System.out.printf("  surviving-direction tests are %.1f%% of all current scan tests (%d / %d)%n",
                    pct(survivorCurrentTests, currentTests), survivorCurrentTests, currentTests);
            System.out.printf("  per surviving direction: current tests %.1f (all neighbors) vs DCLM candidates %.1f%n",
                    survivingDirections == 0 ? 0 : (double) survivorCurrentTests / survivingDirections,
                    survivingDirections == 0 ? 0 : (double) survivorDclmCandidates / survivingDirections);
            System.out.printf("  DCLM would do %.1f%% of the survivor tests (%d / %d) -> survivor saving %.1f%% of ALL current tests%n",
                    pct(survivorDclmCandidates, survivorCurrentTests), survivorDclmCandidates, survivorCurrentTests,
                    pct(survivorCurrentTests - survivorDclmCandidates, currentTests));
            System.out.println();
        }

        private static double pct(long a, long b) { return b == 0 ? 0 : 100.0 * a / b; }
    }

    private final int tesslevel, lutGrid, lutBins;
    private final Stats stats;
    private DirectionMapping mapping;
    private CapDirectionLut lut;
    private long[] maskScratch;
    private boolean[] buried;
    private int[] candidateCount;

    DclmFeasibilityScan(int tesslevel, int lutGrid, int lutBins, Stats stats) {
        this.tesslevel = tesslevel; this.lutGrid = lutGrid; this.lutBins = lutBins; this.stats = stats;
        stats.lutGrid = lutGrid; stats.lutBins = lutBins;
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
            candidateCount = new int[m.numDir];
            lut = new CapDirectionLut(m.ddx, m.ddy, m.ddz, m.numDir, lutGrid, lutBins);
            maskScratch = new long[lut.words];
        }
        double[] ddx = m.ddx, ddy = m.ddy, ddz = m.ddz;
        int numDir = m.numDir;
        boolean[] buried = this.buried;
        int[] candidateCount = this.candidateCount;
        for (int d = 0; d < numDir; d++) candidateCount[d] = 0;

        // (1) production direction-major scan (hint + early-exit) for the cost baseline + buried verdict
        int last = -1;
        long currentTests = 0;
        for (int d = 0; d < numDir; d++) {
            double px = ddx[d], py = ddy[d], pz = ddz[d];
            boolean b = false;
            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                b = true; currentTests++;
            } else {
                for (int k = 0; k < numNeighbors; k++) {
                    currentTests++;
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) { b = true; last = k; break; }
                }
            }
            buried[d] = b;
        }

        // (2) DCLM candidate count per direction = #neighbors whose cap-LUT mask includes the direction
        long[] mask = maskScratch;
        for (int k = 0; k < numNeighbors; k++) {
            lut.maskFor(diffX[k], diffY[k], diffZ[k], thresh[k], mask);
            for (int wi = 0; wi < mask.length; wi++) {
                long bitsW = mask[wi];
                int dbase = wi << 6;
                while (bitsW != 0) {
                    int d = dbase + Long.numberOfTrailingZeros(bitsW);
                    bitsW &= bitsW - 1;
                    if (d < numDir) candidateCount[d]++;
                }
            }
        }

        // (3) accumulate, focusing on surviving directions
        stats.atoms++;
        stats.directions += numDir;
        stats.currentTests += currentTests;
        for (int d = 0; d < numDir; d++) {
            if (!buried[d]) {
                stats.survivingDirections++;
                stats.survivorCurrentTests += numNeighbors;      // current scan tests ALL neighbors for a survivor
                stats.survivorDclmCandidates += candidateCount[d]; // DCLM tests only candidate neighbors
            }
        }

        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
