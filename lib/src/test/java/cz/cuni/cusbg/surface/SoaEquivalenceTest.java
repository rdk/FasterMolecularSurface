package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * Direct equivalence: {@link DevSurfaceV1Soa} must reproduce {@link FasterNumericalSurface}
 * bit-for-bit (it is a pure data-layout change). Asserts identical per-atom areas, identical total
 * point count, and identical surface-point coordinates across the test corpus and a range of
 * (solvent radius, tessellation level) configurations — including tess=2, p2rank's operating point.
 */
class SoaEquivalenceTest {

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void soaMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        FasterNumericalSurface ref = new FasterNumericalSurface(s.load(), solvent, tess);
        DevSurfaceV1Soa    soa = new DevSurfaceV1Soa(s.load(), solvent, tess);
        VariantEquivalence.assertBitForBit(s, solvent, tess, ref, soa);
    }
}
