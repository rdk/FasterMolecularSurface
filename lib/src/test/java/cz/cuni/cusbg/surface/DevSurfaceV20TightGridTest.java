package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

/**
 * {@link DevSurfaceV20TightGrid} (backlog A2: tight-cell neighbor build) must be <b>bit-for-bit identical</b>
 * to {@link DistinctPackedNumericalSurfaceV2}, the surface it re-bases: it changes only the neighbor build,
 * which yields the same neighbor set, so per-atom areas, the surface-point set, and point order are all
 * unchanged. Covers the full corpus × (solvent, tess) matrix — including the Van der Waals path (solvent=0)
 * where boundary ties are densest.
 */
class DevSurfaceV20TightGridTest {

    static Stream<Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void bitForBitVsV2(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess),
                new DevSurfaceV20TightGrid(s.load(), solvent, tess));
    }
}
