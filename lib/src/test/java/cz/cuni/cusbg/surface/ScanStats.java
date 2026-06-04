package cz.cuni.cusbg.surface;

/**
 * Mutable counters accumulated by {@link InstrumentingOcclusionScan} across a corpus, to decide whether
 * the open Tier-A scan/build ideas (A1 bitmask scan, A3 whole-atom reject, A4 coverage-based neighbor
 * prune) have real headroom — measured, not assumed.
 */
final class ScanStats {

    // A1: early-exit / hint effectiveness
    long directions;          // total distinct directions scanned (across all atoms)
    long buriedDirections;
    long hintFires;           // directions resolved by the last-occluder hint (1 test)
    long neighborTests;       // neighbor dot-tests actually performed (with hint + early exit)
    long noExitTests;         // neighbor dot-tests if there were NO early exit (numDir * numNeighbors)

    // A3: whole-atom reject
    long atoms;
    long fullyBuriedAtoms;            // atoms with zero surviving directions (emit no points)
    long neighborTestsFullyBuried;   // neighbor-tests spent on fully-buried atoms (their scan-time share)

    // A4: coverage-based neighbor prune
    long neighborSlots;       // sum of numNeighbors over atoms
    long zeroBuryNeighbors;   // neighbors that bury NO direction of their atom

    void report(String title) {
        System.out.printf("%n=== %s ===%n", title);
        System.out.printf("atoms=%d  directions(distinct)=%d  avg neighbors/atom=%.1f%n",
                atoms, directions, atoms == 0 ? 0 : (double) neighborSlots / atoms);

        System.out.printf("%n[A1] early-exit / hint effectiveness (scan = %.0f%% of CPU at tess 2 / 16t):%n", 34.0);
        System.out.printf("  buried directions:        %.1f%%%n", pct(buriedDirections, directions));
        System.out.printf("  resolved by hint (1 test): %.1f%% of directions%n", pct(hintFires, directions));
        System.out.printf("  avg neighbor-tests / direction (with early exit): %.2f%n",
                directions == 0 ? 0 : (double) neighborTests / directions);
        System.out.printf("  early exit saves: %.1f%% of neighbor-tests vs no-exit (%d vs %d)%n",
                pct(noExitTests - neighborTests, noExitTests), neighborTests, noExitTests);

        System.out.printf("%n[A3] whole-atom reject potential:%n");
        System.out.printf("  fully-buried atoms:        %.1f%% (%d / %d)%n",
                pct(fullyBuriedAtoms, atoms), fullyBuriedAtoms, atoms);
        System.out.printf("  their share of neighbor-tests: %.1f%% (the scan time a perfect skip would save)%n",
                pct(neighborTestsFullyBuried, neighborTests));

        System.out.printf("%n[A4] coverage-based neighbor-prune potential:%n");
        System.out.printf("  neighbors that bury NO direction: %.1f%% (%d / %d)%n",
                pct(zeroBuryNeighbors, neighborSlots), zeroBuryNeighbors, neighborSlots);
        System.out.println();
    }

    private static double pct(long a, long b) {
        return b == 0 ? 0 : 100.0 * a / b;
    }
}
