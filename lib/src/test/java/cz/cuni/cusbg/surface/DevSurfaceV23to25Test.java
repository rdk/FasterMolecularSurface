package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Correctness gates for the three "other-ideas" Dev rungs (C1–C3), each by its fidelity tier:
 * <ul>
 *   <li>V23 (region-binned hint, C3) — <b>bit-for-bit</b> identical to {@code DistinctPackedNumericalSurfaceV3}.</li>
 *   <li>V25 (padded float scan, C2) — <b>bit-for-bit</b> identical to {@code FloatNumericalSurface}.</li>
 *   <li>V24 (float-precision build, C1) — <b>tolerance</b> vs V3: total-area and per-atom area within the
 *       same bounds as the float surface; the float build only reclassifies near-cutoff (near-empty-cap) pairs.</li>
 * </ul>
 */
class DevSurfaceV23to25Test {

    static Stream<Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void v23BitForBitVsV3(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new DistinctPackedNumericalSurfaceV3(s.load(), solvent, tess),
                new DevSurfaceV23RegionHint(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void v25BitForBitVsFloat(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FloatNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV25FloatPaddedScan(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void v24FloatBuildWithinTolerance(TestStructures.Structure s, double solvent, int tess) {
        SurfaceAccuracy.Report r = SurfaceAccuracy.compare(
                new DevSurfaceV24FloatBuild(s.load(), solvent, tess),
                new DistinctPackedNumericalSurfaceV3(s.load(), solvent, tess));
        assertTrue(r.areaTotalRelErr() <= 1e-3,
                () -> s + " total-area rel err " + r.areaTotalRelErr() + " exceeds 1e-3");
        assertTrue(r.areaAtomMaxRelErr() <= 5e-2,
                () -> s + " per-atom area max rel err " + r.areaAtomMaxRelErr() + " exceeds 5e-2");
    }
}
