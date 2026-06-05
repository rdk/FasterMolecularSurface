package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Phase 3 (Lead A5 / DCLM) feasibility: measures the surviving-direction cost share of the scan and how
 * much a perfectly-binned DCLM could shrink it, at tess 2 and 3. Count-based, load-immune; prints only.
 * See {@code autoresearch/LOG.md} Phase 3.
 *
 * <p>Run: {@code ./gradlew scorecard --tests '*DclmFeasibilityTest'}
 */
@Tag("scorecard")
class DclmFeasibilityTest {

    private static final double SOLVENT = 1.4;

    private static final class InstrumentedSurface extends DevSurfaceV1Soa {
        InstrumentedSurface(IAtomContainer mol, double solvent, int tess, int grid, int bins,
                            DclmFeasibilityScan.Stats stats) {
            super(mol, solvent, tess,
                    (atoms, ax, ay, az, radius) ->
                            new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solvent, VdwRadiusCache::get),
                    NeighborOrdering.NONE,
                    new DclmFeasibilityScan(tess, grid, bins, stats),
                    TessellationProvider.CACHED,
                    DistinctFlatSurfacePointStoreV2::new,
                    false,
                    VdwRadiusCache::get,
                    true);
        }
    }

    @Test
    void instrument() {
        for (int tess : new int[]{2, 3}) {
            DclmFeasibilityScan.Stats stats = new DclmFeasibilityScan.Stats();
            for (TestStructures.Structure s : TestStructures.Structure.values()) {
                new InstrumentedSurface(s.load(), SOLVENT, tess, 16, 32, stats);   // finest LUT (tightest candidate sets)
            }
            stats.report("Phase 3 A5/DCLM survivor feasibility @ tess " + tess);
        }
    }
}
