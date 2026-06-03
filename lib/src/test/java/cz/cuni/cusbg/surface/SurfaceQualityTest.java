package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Family-relative invariants on the intrinsic {@link SurfaceQuality} metrics — a regression guard on the
 * point distribution that needs no oracle, so it also covers (future) sampling surfaces. Runs on the
 * smallest structure at tess 2 to keep the O(points²) measurement cheap.
 */
class SurfaceQualityTest {

    private static final double SOLVENT = 1.4;
    private static final int TESS = 2;

    @Test
    void distinctSurfaceHasEssentiallyNoDuplicates() {
        SurfaceQuality.Report q = SurfaceQuality.measure(
                new DistinctPackedNumericalSurfaceV2(TestStructures.Structure.CRAMBIN.load(), SOLVENT, TESS));
        assertTrue(q.totalPoints() > 0, "non-empty surface");
        assertTrue(q.duplicateRatio() < 1e-3,
                "distinct surface should emit ~no coincident duplicates, got " + q.duplicateRatio());
        assertTrue(Double.isFinite(q.evennessCV()) && q.evennessCV() >= 0, "evenness CV is a finite, non-negative number");
        assertTrue(q.minMeanNNRatio() > 0 && q.minMeanNNRatio() <= 1.0001, "min/mean NN ratio in (0,1]");
    }

    @Test
    void fullMultiplicitySurfaceIsMostlyDuplicates() {
        SurfaceQuality.Report q = SurfaceQuality.measure(
                new PackedNumericalSurface(TestStructures.Structure.CRAMBIN.load(), SOLVENT, TESS));
        // CDK's icosahedral tessellation repeats each direction ~5.7x at the same position, so most
        // points are coincident with another — the metric that distinguishes the full-multiplicity
        // family from the distinct one.
        assertTrue(q.duplicateRatio() > 0.5,
                "full-multiplicity surface should be mostly coincident points, got " + q.duplicateRatio());
    }
}
