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
 *   <li>{@link GridSoaNumericalSurface} reproduces {@link FasterNumericalSurface} bit-for-bit
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
        CellGridNeighborList grid = new CellGridNeighborList(ax, ay, az, radius);

        IntArrayList a = new IntArrayList(), b = new IntArrayList();
        for (int i = 0; i < n; i++) {
            final int fi = i;
            a.clear(); b.clear();
            ref.getNeighborsInto(fi, a);
            grid.getNeighborsInto(fi, b);
            assertEquals(toSet(a), toSet(b), () -> s + ": neighbor set mismatch at atom " + fi);
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void gridMatchesFasterExactly(TestStructures.Structure s) {
        assertMatchesFaster(s, new GridSoaNumericalSurface(s.load()));
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void orderedGridMatchesFasterExactly(TestStructures.Structure s) {
        // neighbor reordering must not change any result
        assertMatchesFaster(s, new OrderedGridSoaNumericalSurface(s.load()));
    }

    private static void assertMatchesFaster(TestStructures.Structure s, MolecularSurface variant) {
        FasterNumericalSurface ref = new FasterNumericalSurface(s.load());

        double[] ra = ref.getAllSurfaceAreas(), ga = variant.getAllSurfaceAreas();
        assertEquals(ra.length, ga.length, () -> s + " area array length");
        for (int i = 0; i < ra.length; i++) assertEquals(ra[i], ga[i], 0.0, "per-atom area[" + i + "] for " + s);

        Point3d[] rp = ref.getAllSurfacePoints(), gp = variant.getAllSurfacePoints();
        assertEquals(rp.length, gp.length, () -> s + " total surface-point count");
        for (int i = 0; i < rp.length; i++) {
            assertEquals(rp[i].x, gp[i].x, 0.0, "point[" + i + "].x for " + s);
            assertEquals(rp[i].y, gp[i].y, 0.0, "point[" + i + "].y for " + s);
            assertEquals(rp[i].z, gp[i].z, 0.0, "point[" + i + "].z for " + s);
        }
    }

    private static Set<Integer> toSet(IntArrayList list) {
        Set<Integer> set = new HashSet<>();
        for (int k = 0; k < list.size(); k++) set.add(list.get(k));
        return set;
    }
}
