package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

/**
 * The recommended {@link DistinctPackedNumericalSurfaceV3} (V2 + the SIMD neighbor build) must be
 * <b>bit-for-bit identical</b> to {@link DistinctPackedNumericalSurfaceV2}: it changes only how the
 * build's distance test is evaluated, so areas, the surface-point set, and point order are unchanged.
 * Full corpus × (solvent, tess) matrix, including the Van der Waals path (solvent=0).
 */
class DistinctPackedNumericalSurfaceV3Test {

    static Stream<Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void bitForBitVsV2(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess),
                new DistinctPackedNumericalSurfaceV3(s.load(), solvent, tess));
    }
}
