package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

/**
 * {@link DevSurfaceV22PaddedTail} (backlog A6 + A7 stacked) must be <b>bit-for-bit identical</b> to
 * {@link DistinctPackedNumericalSurfaceV2}: the SIMD build yields the same neighbor set, and the padded
 * scan's non-burying sentinels cannot change any direction's verdict. Full corpus × (solvent, tess)
 * matrix, including the Van der Waals path (solvent=0).
 */
class DevSurfaceV22PaddedTailTest {

    static Stream<Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void bitForBitVsV2(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess),
                new DevSurfaceV22PaddedTail(s.load(), solvent, tess));
    }
}
