package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Lead #1 (LUT-bitmask / fully-buried-atom prize) feasibility instrumentation. Runs the
 * {@link BitmaskFeasibilityScan} over the corpus at tess 2 (p2rank's operating point) and tess 3,
 * printing the idealized floor / budget for a neighbor-major cap→direction bitmask scan vs today's
 * direction-major early-exit scan. Count-based and load-immune; prints only. See {@code autoresearch/LOG.md}
 * Phase 1.
 *
 * <p>Run: {@code ./gradlew scorecard --tests '*BitmaskFeasibilityTest'}
 */
@Tag("scorecard")
class BitmaskFeasibilityTest {

    private static final double SOLVENT = 1.4;

    private static final class InstrumentedSurface extends DevSurfaceV1Soa {
        InstrumentedSurface(IAtomContainer mol, double solvent, int tess, int grid, int bins,
                            BitmaskFeasibilityScan.Stats stats) {
            super(mol, solvent, tess,
                    (atoms, ax, ay, az, radius) ->
                            new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solvent, VdwRadiusCache::get),
                    NeighborOrdering.NONE,
                    new BitmaskFeasibilityScan(tess, grid, bins, stats),
                    TessellationProvider.CACHED,
                    DistinctFlatSurfacePointStoreV2::new,
                    false,
                    VdwRadiusCache::get,
                    true);
        }
    }

    /** Sweep a few LUT resolutions per tess level; coarser = cheaper lookup but more false positives. */
    private static final int[][] CONFIGS = {{6, 8}, {8, 16}, {12, 24}, {16, 32}};

    @Test
    void instrument() {
        for (int tess : new int[]{2, 3}) {
            for (int[] cfg : CONFIGS) {
                int grid = cfg[0], bins = cfg[1];
                BitmaskFeasibilityScan.Stats stats = new BitmaskFeasibilityScan.Stats();
                for (TestStructures.Structure s : TestStructures.Structure.values()) {
                    new InstrumentedSurface(s.load(), SOLVENT, tess, grid, bins, stats);
                }
                stats.report("Lead #1 bitmask feasibility @ tess " + tess
                        + " (LUT grid=" + grid + " bins=" + bins + ")");
            }
        }
    }
}
