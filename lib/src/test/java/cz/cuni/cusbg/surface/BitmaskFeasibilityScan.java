package cz.cuni.cusbg.surface;

/**
 * Feasibility instrumentation for Lead #1 (the LUT-bitmask / fully-buried-atom prize, backlog §1b).
 * Wired into V2's exact pipeline like {@link InstrumentingOcclusionScan}, it produces a correct surface
 * (so the build runs normally) while measuring the cost a <em>neighbor-major</em> cap→direction bitmask
 * scan would pay, in directly comparable exact-dot-test units, against today's direction-major scan.
 *
 * <p>Per atom it computes three numbers:
 * <ul>
 *   <li><b>current</b> — the production direction-major scan (last-occluder hint + per-direction
 *       early exit). This is what the bitmask scan must beat.</li>
 *   <li><b>perfect-mask floor</b> — with a free, exact cap mask, neighbor-major does exactly one
 *       confirming exact-test per buried direction (survivors never tested, confirmed directions never
 *       re-tested), so the floor equals {@code buriedDirections}. This is the best the idea could ever do.</li>
 *   <li><b>no-mask</b> — neighbor-major with no mask at all: each processed neighbor exact-tests every
 *       still-unconfirmed direction, stopping once all directions are confirmed (the atom-level
 *       early-exit). The upper bound if the cap mask were useless.</li>
 * </ul>
 * It also records, for the fully-buried atoms (the ones the atom-level early-exit fires on), how many
 * neighbors are processed before all directions are confirmed — the per-atom budget a real cap-mask
 * build must fit inside. Prints only; asserts nothing.
 */
final class BitmaskFeasibilityScan implements OcclusionScan {

    /** Mutable corpus-wide counters; split into the fully-buried subset and the whole population. */
    static final class Stats {
        long atoms, fullyBuriedAtoms;
        long neighborSlots;              // sum of numNeighbors over all atoms
        long neighborSlotsFullyBuried;   // ... over fully-buried atoms only

        long currentTests;               // direction-major (hint + early exit) exact dot-tests
        long currentTestsFullyBuried;
        long buriedDirections;           // == perfect-mask neighbor-major floor (exact dot-tests)
        long buriedDirectionsFullyBuried;
        long noMaskTests;                // neighbor-major, no mask, atom-level early-exit
        long noMaskTestsFullyBuried;

        long neighborsProcessedFullyBuried;   // neighbors touched before all dirs confirmed (mask budget)

        // real cap-mask LUT (neighbor-major, mask narrows candidates, atom-level early-exit)
        int lutGrid, lutBins;
        long realMaskTests;              // exact dot-tests done with the real LUT mask
        long realMaskTestsFullyBuried;
        long realMaskCandidates;         // sum of set bits in the masks of processed neighbors (false positives = these - confirmations)
        long realMaskConfirmations;      // candidates that actually buried (== buriedDirections if no early-exit waste)
        long realMaskNeighborsProcessed; // neighbors whose mask was fetched (mask-build invocations)
        long soundnessViolations;        // truly-buried (neighbor,direction) pairs MISSING from the mask (must be 0)

        void report(String title) {
            System.out.printf("%n=== %s ===%n", title);
            System.out.printf("atoms=%d  fully-buried=%d (%.1f%%)  avg neighbors/atom=%.1f%n",
                    atoms, fullyBuriedAtoms, pct(fullyBuriedAtoms, atoms),
                    atoms == 0 ? 0 : (double) neighborSlots / atoms);

            System.out.printf("%n[exact dot-tests, whole corpus]%n");
            System.out.printf("  current (direction-major, hint+early-exit): %d%n", currentTests);
            System.out.printf("  perfect-mask floor (= buriedDirections):    %d  (%.1f%% of current)%n",
                    buriedDirections, pct(buriedDirections, currentTests));
            System.out.printf("  no-mask neighbor-major:                      %d  (%.1f%% of current)%n",
                    noMaskTests, pct(noMaskTests, currentTests));
            System.out.printf("  --> headroom gap current-floor:              %d exact-tests (%.1f%% of current)%n",
                    currentTests - buriedDirections, pct(currentTests - buriedDirections, currentTests));

            System.out.printf("%n[fully-buried atoms only — where the atom-level early-exit pays]%n");
            System.out.printf("  current tests on them:        %d  (%.1f%% of all current tests)%n",
                    currentTestsFullyBuried, pct(currentTestsFullyBuried, currentTests));
            System.out.printf("  perfect-mask floor on them:   %d  (%.1f%% of their current)%n",
                    buriedDirectionsFullyBuried, pct(buriedDirectionsFullyBuried, currentTestsFullyBuried));
            System.out.printf("  no-mask on them:              %d  (%.1f%% of their current)%n",
                    noMaskTestsFullyBuried, pct(noMaskTestsFullyBuried, currentTestsFullyBuried));
            System.out.printf("  neighbors processed until fully confirmed: %d of %d slots (%.1f%%)%n",
                    neighborsProcessedFullyBuried, neighborSlotsFullyBuried,
                    pct(neighborsProcessedFullyBuried, neighborSlotsFullyBuried));
            System.out.printf("  --> mask-build budget: %.2f exact-test-equivalents per processed neighbor%n",
                    neighborsProcessedFullyBuried == 0 ? 0
                            : (double) (currentTestsFullyBuried - buriedDirectionsFullyBuried)
                              / neighborsProcessedFullyBuried);

            System.out.printf("%n[real cap-mask LUT, grid=%d bins=%d]  soundness violations=%d (MUST be 0)%n",
                    lutGrid, lutBins, soundnessViolations);
            System.out.printf("  real-mask exact dot-tests:    %d  (%.1f%% of current, floor is %.1f%%)%n",
                    realMaskTests, pct(realMaskTests, currentTests), pct(buriedDirections, currentTests));
            System.out.printf("  real-mask on fully-buried:    %d  (%.1f%% of their current)%n",
                    realMaskTestsFullyBuried, pct(realMaskTestsFullyBuried, currentTestsFullyBuried));
            System.out.printf("  mask candidates examined:     %d  (false-positive rate %.1f%%: %d wasted of %d)%n",
                    realMaskCandidates, pct(realMaskCandidates - realMaskConfirmations, realMaskCandidates),
                    realMaskCandidates - realMaskConfirmations, realMaskCandidates);
            System.out.printf("  avg mask popcount per processed neighbor: %.1f directions — the per-neighbor scan work%n",
                    realMaskNeighborsProcessed == 0 ? 0 : (double) realMaskCandidates / realMaskNeighborsProcessed);
            System.out.println();
        }

        private static double pct(long a, long b) { return b == 0 ? 0 : 100.0 * a / b; }
    }

    private final int tesslevel;
    private final int lutGrid, lutBins;
    private final Stats stats;
    private DirectionMapping mapping;
    private CapDirectionLut lut;
    private long[] maskScratch;
    private boolean[] buried;
    private boolean[] confirmed;
    private boolean[] confirmed2;

    BitmaskFeasibilityScan(int tesslevel, int lutGrid, int lutBins, Stats stats) {
        this.tesslevel = tesslevel;
        this.lutGrid = lutGrid;
        this.lutBins = lutBins;
        this.stats = stats;
        stats.lutGrid = lutGrid;
        stats.lutBins = lutBins;
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
            confirmed = new boolean[m.numDir];
            confirmed2 = new boolean[m.numDir];
            lut = new CapDirectionLut(m.ddx, m.ddy, m.ddz, m.numDir, lutGrid, lutBins);
            maskScratch = new long[lut.words];
        }
        double[] ddx = m.ddx, ddy = m.ddy, ddz = m.ddz;
        int numDir = m.numDir;
        boolean[] buried = this.buried;

        // (1) current production scan: direction-major, last-occluder hint + early exit.
        int last = -1;
        long currentTests = 0;
        long buriedDirs = 0;
        for (int d = 0; d < numDir; d++) {
            double px = ddx[d], py = ddy[d], pz = ddz[d];
            boolean b = false;
            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                b = true;
                currentTests++;
            } else {
                for (int k = 0; k < numNeighbors; k++) {
                    currentTests++;
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                        b = true; last = k; break;
                    }
                }
            }
            buried[d] = b;
            if (b) buriedDirs++;
        }
        boolean fullyBuried = (buriedDirs == numDir);

        // (2) no-mask neighbor-major: each neighbor exact-tests every still-unconfirmed direction,
        //     stop once all directions confirmed (the atom-level early-exit the idea is built on).
        boolean[] confirmed = this.confirmed;
        for (int d = 0; d < numDir; d++) confirmed[d] = false;
        int remaining = numDir;
        long noMaskTests = 0;
        int neighborsProcessed = 0;
        for (int k = 0; k < numNeighbors && remaining > 0; k++) {
            neighborsProcessed++;
            double nx = diffX[k], ny = diffY[k], nz = diffZ[k], th = thresh[k];
            for (int d = 0; d < numDir; d++) {
                if (confirmed[d]) continue;
                noMaskTests++;
                if (nx * ddx[d] + ny * ddy[d] + nz * ddz[d] > th) {
                    confirmed[d] = true;
                    remaining--;
                }
            }
        }

        // (3) real cap-mask LUT neighbor-major: mask narrows candidates, atom-level early-exit,
        //     plus a full soundness sweep (every truly-buried pair must be in the mask).
        boolean[] confirmed2 = this.confirmed2;
        for (int d = 0; d < numDir; d++) confirmed2[d] = false;
        int remaining2 = numDir;
        long realTests = 0, realCandidates = 0, realConfirmations = 0;
        int realNeighborsProcessed = 0;
        long[] mask = maskScratch;
        for (int k = 0; k < numNeighbors; k++) {
            double nx = diffX[k], ny = diffY[k], nz = diffZ[k], th = thresh[k];
            lut.maskFor(nx, ny, nz, th, mask);

            // soundness: every direction this neighbor actually buries must be present in the mask
            for (int d = 0; d < numDir; d++) {
                boolean exact = nx * ddx[d] + ny * ddy[d] + nz * ddz[d] > th;
                if (exact && (mask[d >>> 6] & (1L << (d & 63))) == 0) stats.soundnessViolations++;
            }

            if (remaining2 > 0) {   // candidate scan + atom-level early-exit (measured cost)
                realNeighborsProcessed++;
                for (int wi = 0; wi < mask.length; wi++) {
                    long bitsW = mask[wi];
                    int dbase = wi << 6;
                    while (bitsW != 0) {
                        int d = dbase + Long.numberOfTrailingZeros(bitsW);
                        bitsW &= bitsW - 1;
                        if (d >= numDir) break;
                        realCandidates++;
                        if (confirmed2[d]) continue;
                        realTests++;
                        if (nx * ddx[d] + ny * ddy[d] + nz * ddz[d] > th) {
                            confirmed2[d] = true; remaining2--; realConfirmations++;
                        }
                    }
                }
            }
        }

        // accumulate
        stats.realMaskTests += realTests;
        stats.realMaskCandidates += realCandidates;
        stats.realMaskConfirmations += realConfirmations;
        stats.realMaskNeighborsProcessed += realNeighborsProcessed;
        if (fullyBuried) stats.realMaskTestsFullyBuried += realTests;

        stats.atoms++;
        stats.neighborSlots += numNeighbors;
        stats.currentTests += currentTests;
        stats.buriedDirections += buriedDirs;   // == perfect-mask floor exact-tests
        stats.noMaskTests += noMaskTests;
        if (fullyBuried) {
            stats.fullyBuriedAtoms++;
            stats.neighborSlotsFullyBuried += numNeighbors;
            stats.currentTestsFullyBuried += currentTests;
            stats.buriedDirectionsFullyBuried += buriedDirs;
            stats.noMaskTestsFullyBuried += noMaskTests;
            stats.neighborsProcessedFullyBuried += neighborsProcessed;
        }

        m.emitWeighted(buried, totalRadius, atomX, atomY, atomZ, sink);
    }
}
