package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

/**
 * {@link DevSurfaceV21SimdBuild} (backlog A6: SIMD-vectorized neighbor-build distance pass) must be
 * <b>bit-for-bit identical</b> to {@link DistinctPackedNumericalSurfaceV2}, the surface it re-bases: it
 * changes only how the distance test is evaluated (a lane at a time vs scalar), yielding the same neighbor
 * set, so per-atom areas, the surface-point set, and point order are unchanged. Covers the full corpus ×
 * (solvent, tess) matrix, including the Van der Waals path (solvent=0).
 */
class DevSurfaceV21SimdBuildTest {

    static Stream<Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void bitForBitVsV2(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess),
                new DevSurfaceV21SimdBuild(s.load(), solvent, tess));
    }
}
