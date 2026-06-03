package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import javax.vecmath.Point3d;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for {@link DistinctPackedNumericalSurface}: it emits one point per distinct surviving direction
 * (no coincident duplicates) yet computes the SAME area as the full-multiplicity surfaces.
 */
class DistinctPackedNumericalSurfaceTest {

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    /** Areas (per-atom and total) must be bit-for-bit identical to FasterNumericalSurface. */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void areasMatchFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        FasterNumericalSurface faster = new FasterNumericalSurface(s.load(), solvent, tess);
        DistinctPackedNumericalSurface distinct = new DistinctPackedNumericalSurface(s.load(), solvent, tess);

        assertEquals(faster.getTotalSurfaceArea(), distinct.getTotalSurfaceArea(), 0.0d, "total area");
        double[] fa = faster.getAllSurfaceAreas(), da = distinct.getAllSurfaceAreas();
        assertEquals(fa.length, da.length, "atom count");
        for (int i = 0; i < fa.length; i++) assertEquals(fa[i], da[i], 0.0d, "area[" + i + "]");
    }

    /** The distinct surface has far fewer points, each of which is one of FasterNumericalSurface's points. */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void emitsFewerPointsAllPresentInFaster(TestStructures.Structure s, double solvent, int tess) {
        Point3d[] fasterPts = new FasterNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();
        Point3d[] distinctPts = new DistinctPackedNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();

        assertTrue(distinctPts.length > 0, "non-empty");
        assertTrue(distinctPts.length < fasterPts.length, "distinct must be fewer than full multiplicity");
        // tessellation emits each direction with multiplicity >=5, so dedup should remove most points
        assertTrue((double) fasterPts.length / distinctPts.length > 3.0, "expected large dedup ratio");

        Set<String> fasterSet = new HashSet<>(fasterPts.length * 2);
        for (Point3d p : fasterPts) fasterSet.add(p.x + "," + p.y + "," + p.z);
        for (Point3d p : distinctPts) {
            assertTrue(fasterSet.contains(p.x + "," + p.y + "," + p.z),
                    "distinct point must exist in the full surface: " + p);
        }
        // distinct points themselves must be unique (no coincident duplicates left)
        Set<String> distinctSet = new HashSet<>(distinctPts.length * 2);
        for (Point3d p : distinctPts) {
            assertTrue(distinctSet.add(p.x + "," + p.y + "," + p.z), "distinct points must be unique: " + p);
        }
    }

    /**
     * Strong per-atom equivalence to its OWN original ({@link PackedNumericalSurface}): identical total
     * area, identical per-atom area, and FOR EACH ATOM the distinct point set equals the de-duplicated
     * original point set exactly - both directions, so nothing the original has is missing and nothing
     * extra is added - with no coincident duplicates left in the distinct surface.
     */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void perAtomEquivalentToPackedOriginal(TestStructures.Structure s, double solvent, int tess) throws Exception {
        PackedNumericalSurface orig = new PackedNumericalSurface(s.load(), solvent, tess);
        DistinctPackedNumericalSurface dist = new DistinctPackedNumericalSurface(s.load(), solvent, tess);

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

    /**
     * The weighted scan itself is correct: paired with a normal (full-multiplicity) store its default
     * {@code addWeighted} re-emits each direction's multiplicity, producing the SAME MULTISET of points
     * as FasterNumericalSurface (same set of locations, each with the same multiplicity). The ORDER
     * differs (the weighted scan groups a direction's copies, vs the original interleaved tessellation
     * order), which is irrelevant for the surface, so we compare as a multiset rather than bit-for-bit.
     */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void weightedScanWithFullStoreMatchesFasterAsMultiset(TestStructures.Structure s, double solvent, int tess) {
        DevSurfaceV1Soa fullViaWeighted = new DevSurfaceV1Soa(s.load(), solvent, tess,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solvent, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                new GlobalDedupWeightedOcclusionScan(tess),
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);
        Point3d[] faster = new FasterNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();
        Point3d[] weighted = fullViaWeighted.getAllSurfacePoints();
        assertEquals(multiset(faster), multiset(weighted), "weighted-scan full-multiplicity multiset must equal faster's");
    }

    private static String key(Point3d p) { return p.x + "," + p.y + "," + p.z; }

    private static Set<String> keys(Point3d[] pts) {
        Set<String> set = new HashSet<>(pts.length * 2);
        for (Point3d p : pts) set.add(key(p));
        return set;
    }

    private static java.util.Map<String, Integer> multiset(Point3d[] pts) {
        java.util.Map<String, Integer> m = new java.util.HashMap<>(pts.length * 2);
        for (Point3d p : pts) m.merge(p.x + "," + p.y + "," + p.z, 1, Integer::sum);
        return m;
    }
}
