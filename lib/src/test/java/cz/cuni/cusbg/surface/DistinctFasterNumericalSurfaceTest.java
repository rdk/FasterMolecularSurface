package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import javax.vecmath.Point3d;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for {@link DistinctFasterNumericalSurface}: areas bit-identical to {@link FasterNumericalSurface},
 * distinct (de-duplicated) point set, and identical output to {@link DistinctPackedNumericalSurface}.
 */
class DistinctFasterNumericalSurfaceTest {

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void areasMatchFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        FasterNumericalSurface faster = new FasterNumericalSurface(s.load(), solvent, tess);
        DistinctFasterNumericalSurface distinct = new DistinctFasterNumericalSurface(s.load(), solvent, tess);

        assertEquals(faster.getTotalSurfaceArea(), distinct.getTotalSurfaceArea(), 0.0d, "total area");
        double[] fa = faster.getAllSurfaceAreas(), da = distinct.getAllSurfaceAreas();
        assertEquals(fa.length, da.length, "atom count");
        for (int i = 0; i < fa.length; i++) assertEquals(fa[i], da[i], 0.0d, "area[" + i + "]");
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void emitsFewerUniquePointsSubsetOfFaster(TestStructures.Structure s, double solvent, int tess) {
        Point3d[] fasterPts = new FasterNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();
        Point3d[] distinctPts = new DistinctFasterNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();

        assertTrue(distinctPts.length > 0, "non-empty");
        assertTrue(distinctPts.length < fasterPts.length, "distinct must be fewer than full multiplicity");
        assertTrue((double) fasterPts.length / distinctPts.length > 3.0, "expected large dedup ratio");

        Set<String> fasterSet = keys(fasterPts);
        Set<String> distinctSet = new HashSet<>(distinctPts.length * 2);
        for (Point3d p : distinctPts) {
            assertTrue(fasterSet.contains(key(p)), "distinct point must exist in the full surface: " + p);
            assertTrue(distinctSet.add(key(p)), "distinct points must be unique: " + p);
        }
    }

    /**
     * Strong per-atom equivalence to its OWN original ({@link FasterNumericalSurface}): identical total
     * area, identical per-atom area, and FOR EACH ATOM the distinct point set equals the de-duplicated
     * original point set exactly (both directions), with no coincident duplicates left.
     */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void perAtomEquivalentToFasterOriginal(TestStructures.Structure s, double solvent, int tess) throws Exception {
        FasterNumericalSurface orig = new FasterNumericalSurface(s.load(), solvent, tess);
        DistinctFasterNumericalSurface dist = new DistinctFasterNumericalSurface(s.load(), solvent, tess);

        assertEquals(orig.getTotalSurfaceArea(), dist.getTotalSurfaceArea(), 0.0d, "total area");
        double[] oa = orig.getAllSurfaceAreas(), da = dist.getAllSurfaceAreas();
        assertEquals(oa.length, da.length, "atom count");

        for (int i = 0; i < oa.length; i++) {
            assertEquals(oa[i], da[i], 0.0d, "area[" + i + "]");

            Set<String> origUnique = keys(orig.getSurfacePoints(i));   // de-duplicated original locations
            Point3d[] distPts = dist.getSurfacePoints(i);
            Set<String> distSet = new HashSet<>(distPts.length * 2);
            for (Point3d p : distPts) {
                assertTrue(distSet.add(key(p)), "atom " + i + ": no coincident duplicates allowed: " + p);
            }
            assertEquals(origUnique, distSet,
                    "atom " + i + ": distinct point set must equal the unique original locations");
        }
    }

    /** The two distinct implementations (this and the Packed-pipeline based one) must agree exactly. */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void matchesDistinctPackedNumericalSurface(TestStructures.Structure s, double solvent, int tess) {
        DistinctFasterNumericalSurface a = new DistinctFasterNumericalSurface(s.load(), solvent, tess);
        DistinctPackedNumericalSurface b = new DistinctPackedNumericalSurface(s.load(), solvent, tess);

        assertEquals(a.getTotalSurfaceArea(), b.getTotalSurfaceArea(), 0.0d, "total area");
        Point3d[] pa = a.getAllSurfacePoints(), pb = b.getAllSurfacePoints();
        assertEquals(pa.length, pb.length, "point count");
        assertEquals(keys(pa), keys(pb), "distinct point sets must be identical");
    }

    private static String key(Point3d p) { return p.x + "," + p.y + "," + p.z; }

    private static Set<String> keys(Point3d[] pts) {
        Set<String> set = new HashSet<>(pts.length * 2);
        for (Point3d p : pts) set.add(key(p));
        return set;
    }
}
