package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Feasibility instrumentation for the open Tier-A scan/build ideas — run before scaffolding them, so the
 * decision is measured rather than assumed (the lesson from A2 and A7, both of which lost on unverified
 * premises). Wires {@link InstrumentingOcclusionScan} into V2's exact pipeline over the corpus at tess 2
 * (p2rank's operating point) and prints the headroom for A1 (early-exit effectiveness), A3 (fully-buried
 * atoms), and A4 (neighbors that bury nothing). Prints only; asserts nothing.
 *
 * <p>Run: {@code ./gradlew scorecard --tests '*ScanInstrumentationTest'}
 */
@Tag("scorecard")
class ScanInstrumentationTest {

    private static final double SOLVENT = 1.4;
    private static final int TESS = 2;

    /** V2's pipeline with the instrumenting scan substituted (diff/thresh computed identically to V2). */
    private static final class InstrumentedSurface extends DevSurfaceV1Soa {
        InstrumentedSurface(IAtomContainer mol, double solvent, int tess, ScanStats stats) {
            super(mol, solvent, tess,
                    (atoms, ax, ay, az, radius) ->
                            new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solvent, VdwRadiusCache::get),
                    NeighborOrdering.NONE,
                    new InstrumentingOcclusionScan(tess, stats),
                    TessellationProvider.CACHED,
                    DistinctFlatSurfacePointStoreV2::new,
                    false,
                    VdwRadiusCache::get,
                    true);
        }
    }

    @Test
    void instrument() {
        ScanStats stats = new ScanStats();
        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            // building the surface runs the instrumenting scan over every atom (eager init)
            new InstrumentedSurface(s.load(), SOLVENT, TESS, stats);
        }
        stats.report("Scan/build feasibility @ tess 2 over the corpus");
    }
}
