package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import javax.vecmath.Point3d;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for {@link DistinctPackedNumericalSurfaceV2}: it must produce the SAME distinct surface as
 * {@link DistinctPackedNumericalSurface} (bit-exact areas, identical per-atom distinct point sets) via
 * the shared {@link DirectionMapping} and the weighted dedup scan (vectorized when the Vector API is
 * available, scalar otherwise), while keeping the bit-exact-area contract against
 * {@link FasterNumericalSurface}.
 */
class DistinctPackedNumericalSurfaceV2Test {

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    /** Areas (per-atom and total) must be bit-for-bit identical to FasterNumericalSurface. */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void areasMatchFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        FasterNumericalSurface faster = new FasterNumericalSurface(s.load(), solvent, tess);
        DistinctPackedNumericalSurfaceV2 distinct = new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess);

        assertEquals(faster.getTotalSurfaceArea(), distinct.getTotalSurfaceArea(), 0.0d, "total area");
        double[] fa = faster.getAllSurfaceAreas(), da = distinct.getAllSurfaceAreas();
        assertEquals(fa.length, da.length, "atom count");
        for (int i = 0; i < fa.length; i++) assertEquals(fa[i], da[i], 0.0d, "area[" + i + "]");
    }

    /** V2 must be equivalent to V1: identical total/per-atom areas and identical per-atom distinct point sets. */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void matchesV1Exactly(TestStructures.Structure s, double solvent, int tess) throws Exception {
        DistinctPackedNumericalSurface v1 = new DistinctPackedNumericalSurface(s.load(), solvent, tess);
        DistinctPackedNumericalSurfaceV2 v2 = new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess);

        assertEquals(v1.getTotalSurfaceArea(), v2.getTotalSurfaceArea(), 0.0d, "total area");
        double[] a1 = v1.getAllSurfaceAreas(), a2 = v2.getAllSurfaceAreas();
        assertEquals(a1.length, a2.length, "atom count");
        for (int i = 0; i < a1.length; i++) {
            assertEquals(a1[i], a2[i], 0.0d, "area[" + i + "]");
            assertEquals(keys(v1.getSurfacePoints(i)), keys(v2.getSurfacePoints(i)),
                    "atom " + i + ": distinct point set must match V1");
        }
    }

    /** The distinct surface has far fewer points, each of which is one of FasterNumericalSurface's points. */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void emitsFewerPointsAllPresentInFaster(TestStructures.Structure s, double solvent, int tess) {
        Point3d[] fasterPts = new FasterNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints();
        Point3d[] distinctPts = new DistinctPackedNumericalSurfaceV2(s.load(), solvent, tess).getAllSurfacePoints();

        assertTrue(distinctPts.length > 0, "non-empty");
        assertTrue(distinctPts.length < fasterPts.length, "distinct must be fewer than full multiplicity");
        assertTrue((double) fasterPts.length / distinctPts.length > 3.0, "expected large dedup ratio");

        Set<String> fasterSet = keys(fasterPts);
        for (Point3d p : distinctPts) {
            assertTrue(fasterSet.contains(key(p)), "distinct point must exist in the full surface: " + p);
        }
        Set<String> distinctSet = new HashSet<>(distinctPts.length * 2);
        for (Point3d p : distinctPts) {
            assertTrue(distinctSet.add(key(p)), "distinct points must be unique: " + p);
        }
    }

    /**
     * Validates the new weighted dedup scans <em>independently of the dedup store</em>: paired with a
     * normal full-multiplicity store ({@link ListSurfacePointStore}), each scan's default
     * {@code addWeighted} re-emits every direction's multiplicity, so the result must be the SAME
     * MULTISET of points as {@link FasterNumericalSurface} (same locations, same multiplicities). The
     * order differs (the weighted scan groups a direction's copies), which is irrelevant, so we compare
     * as a multiset. Exercises the scalar {@link WeightedDedupOcclusionScan} always, and the vectorized
     * {@link Vectorized256WeightedDedupOcclusionScan} when the Vector API is available - the same scans
     * {@link DistinctPackedNumericalSurfaceV2} selects between.
     */
    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void weightedScansWithFullStoreMatchFasterAsMultiset(TestStructures.Structure s, double solvent, int tess) {
        java.util.Map<String, Integer> faster = multiset(new FasterNumericalSurface(s.load(), solvent, tess).getAllSurfacePoints());

        assertEquals(faster, multiset(fullViaScan(s, solvent, tess, new WeightedDedupOcclusionScan(tess)).getAllSurfacePoints()),
                "scalar weighted scan full-multiplicity multiset must equal faster's");

        if (DistinctPackedNumericalSurfaceV2.isVectorized()) {
            assertEquals(faster, multiset(fullViaScan(s, solvent, tess, new Vectorized256WeightedDedupOcclusionScan(tess)).getAllSurfacePoints()),
                    "vectorized weighted scan full-multiplicity multiset must equal faster's");
        }
    }

    /** The V2 compute pipeline with the given scan but a normal (full-multiplicity) list store. */
    private static MolecularSurface fullViaScan(TestStructures.Structure s, double solvent, int tess, OcclusionScan scan) {
        return new DevSurfaceV1Soa(s.load(), solvent, tess,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solvent, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                scan,
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);
    }

    private static String key(Point3d p) { return p.x + "," + p.y + "," + p.z; }

    private static Set<String> keys(Point3d[] pts) {
        Set<String> set = new HashSet<>(pts.length * 2);
        for (Point3d p : pts) set.add(key(p));
        return set;
    }

    private static java.util.Map<String, Integer> multiset(Point3d[] pts) {
        java.util.Map<String, Integer> m = new java.util.HashMap<>(pts.length * 2);
        for (Point3d p : pts) m.merge(key(p), 1, Integer::sum);
        return m;
    }
}
