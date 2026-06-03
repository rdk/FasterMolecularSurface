package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import javax.vecmath.Point3d;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tolerance tests for {@link FloatNumericalSurface}: the float-verdict distinct surface is NOT
 * bit-exact - a tessellation point within float epsilon of the occlusion boundary may flip survival
 * in single precision, which both drops points the exact surface keeps and keeps points it drops - so
 * it is checked against {@link FasterNumericalSurface} / {@link DistinctPackedNumericalSurfaceV2} with
 * tolerances rather than the bit-exact equivalence harness.
 *
 * <p>Tolerances are set with wide margin over the worst observed deviation across the full config
 * matrix (tess 2-4 at solvent 1.4, plus the degenerate van der Waals surface at solvent 0.0, where the
 * most boundary ties occur):
 * <ul>
 *   <li>total surface area: observed &le; 1.4e-5 relative -&gt; assert 1e-4;</li>
 *   <li>per-atom area: observed &le; 0.9% relative (one flip is ~1/atomPointCount) -&gt; assert 2%;</li>
 *   <li>distinct point-set symmetric difference vs the exact V2 surface: observed &le; 1.4e-5 of the
 *       point count (single-digit points out of hundreds of thousands) -&gt; assert 1e-3.</li>
 * </ul>
 * On a JVM without the Vector API the surface falls back to the exact scalar double scan, so deviation
 * is zero and every tolerance passes trivially.
 */
class FloatNumericalSurfaceTest {

    private static final double AREA_TOTAL_REL_TOL = 1e-4;
    private static final double AREA_ATOM_REL_TOL = 2e-2;
    private static final double SET_DIFF_REL_TOL = 1e-3;

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    /** Total and per-atom areas match FasterNumericalSurface within tolerance (verdict-only float error). */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void areasApproxMatchFaster(TestStructures.Structure s, double solvent, int tess) {
        FasterNumericalSurface ref = new FasterNumericalSurface(s.load(), solvent, tess);
        FloatNumericalSurface flt = new FloatNumericalSurface(s.load(), solvent, tess);

        double rt = ref.getTotalSurfaceArea();
        assertEquals(rt, flt.getTotalSurfaceArea(), AREA_TOTAL_REL_TOL * rt, "total area");
        double[] ra = ref.getAllSurfaceAreas(), fa = flt.getAllSurfaceAreas();
        assertEquals(ra.length, fa.length, "atom count");
        for (int i = 0; i < ra.length; i++) {
            assertEquals(ra[i], fa[i], AREA_ATOM_REL_TOL * Math.max(1.0, ra[i]), "area[" + i + "]");
        }
    }

    /**
     * Emitted points are distinct (no coincident duplicates), and the distinct point SET differs from
     * the exact V2 surface by only a tiny boundary-flip fraction. Surviving points have double-computed
     * positions identical to the exact surface; only set membership can differ near the boundary.
     */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void distinctPointSetCloseToExactV2(TestStructures.Structure s, double solvent, int tess) {
        Point3d[] floatPts = new FloatNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();
        Set<String> floatSet = new HashSet<>(floatPts.length * 2);
        for (Point3d p : floatPts) {
            assertTrue(floatSet.add(key(p)), "distinct points must be unique: " + p);
        }
        assertTrue(floatPts.length > 0, "non-empty");

        Set<String> v2Set = keys(new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess).getAllSurfacePoints());
        int d = 0;
        for (String k : floatSet) if (!v2Set.contains(k)) d++;   // points float kept that V2 dropped
        for (String k : v2Set) if (!floatSet.contains(k)) d++;   // points V2 kept that float dropped
        final int diff = d;
        assertTrue(diff <= Math.ceil(SET_DIFF_REL_TOL * v2Set.size()),
                () -> s + " float/V2 point-set symmetric difference " + diff + " exceeds "
                        + (SET_DIFF_REL_TOL * 100) + "% of " + v2Set.size());
    }

    private static String key(Point3d p) { return p.x + "," + p.y + "," + p.z; }

    private static Set<String> keys(Point3d[] pts) {
        Set<String> set = new HashSet<>(pts.length * 2);
        for (Point3d p : pts) set.add(key(p));
        return set;
    }
}
