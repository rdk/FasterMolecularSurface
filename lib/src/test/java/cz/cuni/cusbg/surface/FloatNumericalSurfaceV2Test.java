package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tolerance gate for {@link FloatNumericalSurfaceV2} (float build + float scan stacked). Held to the same
 * bounds {@link FloatNumericalSurface} documents, against the exact distinct surface
 * {@link DistinctPackedNumericalSurfaceV3}: total-area ≤ 1e-4 relative, per-atom area ≤ 2%, distinct
 * point-set symmetric difference ≤ 1e-3 of the point count. Stacking the two float effects must not push
 * the error past the established single-precision envelope.
 */
class FloatNumericalSurfaceV2Test {

    private static final double AREA_TOTAL_REL_TOL = 1e-4;
    private static final double AREA_ATOM_REL_TOL = 2e-2;
    private static final double SET_DIFF_REL_TOL = 1e-3;

    static Stream<Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void withinFloatToleranceVsExactDistinct(TestStructures.Structure s, double solvent, int tess) {
        SurfaceAccuracy.Report r = SurfaceAccuracy.compare(
                new FloatNumericalSurfaceV2(s.load(), solvent, tess),
                new DistinctPackedNumericalSurfaceV3(s.load(), solvent, tess));
        assertTrue(r.areaTotalRelErr() <= AREA_TOTAL_REL_TOL,
                () -> s + " total-area rel err " + r.areaTotalRelErr());
        assertTrue(r.areaAtomMaxRelErr() <= AREA_ATOM_REL_TOL,
                () -> s + " per-atom area max rel err " + r.areaAtomMaxRelErr());
        assertTrue(r.pointSetSymDiffFrac() <= SET_DIFF_REL_TOL,
                () -> s + " point-set symmetric difference " + r.pointSetSymDiffFrac());
    }
}
