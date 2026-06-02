package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import javax.vecmath.Point3d;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Direct equivalence: {@link SoaNumericalSurface} must reproduce {@link FasterNumericalSurface}
 * bit-for-bit (it is a pure data-layout change). Asserts identical per-atom areas, identical total
 * point count, and identical surface-point coordinates across the test corpus.
 */
class SoaEquivalenceTest {

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structures() {
        return TestStructures.all();
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void soaMatchesFasterExactly(TestStructures.Structure s) {
        FasterNumericalSurface ref = new FasterNumericalSurface(s.load());
        SoaNumericalSurface    soa = new SoaNumericalSurface(s.load());

        double[] ra = ref.getAllSurfaceAreas();
        double[] sa = soa.getAllSurfaceAreas();
        assertEquals(ra.length, sa.length, () -> s + " area array length");
        for (int i = 0; i < ra.length; i++) {
            assertEquals(ra[i], sa[i], 0.0, "per-atom area[" + i + "] for " + s);
        }

        Point3d[] rp = ref.getAllSurfacePoints();
        Point3d[] sp = soa.getAllSurfacePoints();
        assertEquals(rp.length, sp.length, () -> s + " total surface-point count");
        for (int i = 0; i < rp.length; i++) {
            assertEquals(rp[i].x, sp[i].x, 0.0, "point[" + i + "].x for " + s);
            assertEquals(rp[i].y, sp[i].y, 0.0, "point[" + i + "].y for " + s);
            assertEquals(rp[i].z, sp[i].z, 0.0, "point[" + i + "].z for " + s);
        }
    }
}
