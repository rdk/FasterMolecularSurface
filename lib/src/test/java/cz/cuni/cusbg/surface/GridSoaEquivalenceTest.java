package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point3d;
import java.util.HashSet;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Equivalence checks for the flat-grid variant (step 2):
 * <ol>
 *   <li>{@link CellGridNeighborList} returns the same neighbor SET as the hash-box
 *       {@link NeighborList} for every atom (order may differ, set must not).</li>
 *   <li>{@link DevSurfaceV2Grid} reproduces {@link FasterNumericalSurface} bit-for-bit
 *       (areas and surface-point coordinates).</li>
 * </ol>
 */
class GridSoaEquivalenceTest {

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structures() {
        return TestStructures.all();
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void gridNeighborSetsMatchNeighborList(TestStructures.Structure s) {
        assertSetsMatchNeighborList(s, (atoms, ax, ay, az, r) -> new CellGridNeighborList(ax, ay, az, r));
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void symmetricGridNeighborSetsMatchNeighborList(TestStructures.Structure s) {
        assertSetsMatchNeighborList(s, (atoms, ax, ay, az, r) -> new SymmetricCellGridNeighborList(ax, ay, az, r));
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void lowAllocSymmetricGridNeighborSetsMatchNeighborList(TestStructures.Structure s) {
        assertSetsMatchNeighborList(s, (atoms, ax, ay, az, r) -> new LowAllocSymmetricCellGridNeighborList(ax, ay, az, r));
    }

    /** Every atom's neighbor SET from {@code factory}'s source must equal the hash-box NeighborList's. */
    private static void assertSetsMatchNeighborList(TestStructures.Structure s, NeighborSourceFactory factory) {
        IAtomContainer mol = s.load();
        org.openscience.cdk.interfaces.IAtom[] atoms = AtomContainerManipulator.getAtomArray(mol);
        int n = atoms.length;

        double solvent = 1.4, maxR = 0;
        double[] ax = new double[n], ay = new double[n], az = new double[n];
        for (int i = 0; i < n; i++) {
            Point3d p = atoms[i].getPoint3d();
            ax[i] = p.x; ay[i] = p.y; az[i] = p.z;
            double r = FasterNumericalSurface.getVdwRadius(atoms[i]) + solvent;
            if (r > maxR) maxR = r;
        }
        double radius = maxR + solvent;

        NeighborList ref = new NeighborList(atoms, radius);
        NeighborSource src = factory.create(atoms, ax, ay, az, radius);

        IntArrayList a = new IntArrayList(), b = new IntArrayList();
        for (int i = 0; i < n; i++) {
            final int fi = i;
            a.clear(); b.clear();
            ref.getNeighborsInto(fi, a);
            src.getNeighborsInto(fi, b);
            assertEquals(toSet(a), toSet(b), () -> s + ": neighbor set mismatch at atom " + fi);
        }
    }

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structureConfigs() {
        return TestStructures.structureConfigs();
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV2GridMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV2Grid(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV3SortedMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // neighbor reordering must not change any result, at any tessellation level
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV3Sorted(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV4HintedMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // the last-occluder-first hint changes only the scan order, never the surviving set
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV4Hinted(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV5SymmetricMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // the symmetric neighbor precompute yields the same neighbor sets, so the result is unchanged
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV5Symmetric(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV6LowAllocMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // the lower-allocation two-pass symmetric precompute still yields the same neighbor sets
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV6LowAlloc(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV7SimdMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // the SIMD scan computes the dot product lane-for-lane (no FMA), so it must match bit-for-bit
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV7Simd(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV8PrunedMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // the per-pair occlusion cutoff only drops neighbors that provably never bury, so dropping them
        // (and the 256-bit scan) cannot change the surface: must match bit-for-bit
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV8Pruned(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV9DedupMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // dedup scans each distinct direction once then re-expands in original point order: the surviving
        // multiset and its order are unchanged, so it must match the per-point scan bit-for-bit
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV9Dedup(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV10CachedMapMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // the process-wide-cached direction mapping is the same mapping, just computed once; result unchanged
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV10CachedMap(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV11CachedTessMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // caching the tessellation arrays holds the same values the engine would build; result unchanged
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV11CachedTess(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV12FlatMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // flat double[] point storage + cached VdW radii change only storage and a memoized lookup,
        // not the coordinates/areas, so the surface must match bit-for-bit
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV12Flat(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV13ArenaMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        // reusing per-thread scratch across builds changes only allocation; buffers are fully overwritten
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV13Arena(s.load(), solvent, tess));
    }

    @ParameterizedTest(name = "{0} solvent={1} tess={2}")
    @MethodSource("structureConfigs")
    void devSurfaceV14ArenaFlatMatchesFasterExactly(TestStructures.Structure s, double solvent, int tess) {
        VariantEquivalence.assertBitForBit(s, solvent, tess,
                new FasterNumericalSurface(s.load(), solvent, tess),
                new DevSurfaceV14ArenaFlat(s.load(), solvent, tess));
    }

    private static Set<Integer> toSet(IntArrayList list) {
        Set<Integer> set = new HashSet<>();
        for (int k = 0; k < list.size(); k++) set.add(list.get(k));
        return set;
    }
}
