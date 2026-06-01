package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point3d;
import java.util.Arrays;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Direct unit tests for {@link NeighborList}.
 *
 * <p>{@code NeighborList} is the spatial-acceleration structure under the surface algorithm and a
 * prime target for future optimization, so it is pinned independently here. The neighbour relation
 * is defined as: atoms whose centres are within {@code 2*radius} (the box size) of atom {@code i},
 * excluding {@code i} itself. The authoritative check is equivalence to a brute-force O(n^2) scan.
 */
class NeighborListTest {

    /** Brute-force reference: all j != i with distance^2 < (2*radius)^2. */
    private static int[] bruteForce(IAtom[] atoms, int i, double radius) {
        double cutoff2 = (2 * radius) * (2 * radius);
        TreeSet<Integer> result = new TreeSet<>();
        Point3d pi = atoms[i].getPoint3d();
        for (int j = 0; j < atoms.length; j++) {
            if (j == i) continue;
            if (atoms[j].getPoint3d().distanceSquared(pi) < cutoff2) {
                result.add(j);
            }
        }
        return result.stream().mapToInt(Integer::intValue).toArray();
    }

    private static int[] sorted(int[] a) {
        int[] c = a.clone();
        Arrays.sort(c);
        return c;
    }

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structures() {
        // Use the smaller structures so the O(n^2) brute-force reference stays fast.
        return java.util.stream.Stream.of(
                org.junit.jupiter.params.provider.Arguments.of(TestStructures.Structure.CRAMBIN),
                org.junit.jupiter.params.provider.Arguments.of(TestStructures.Structure.BPTI),
                org.junit.jupiter.params.provider.Arguments.of(TestStructures.Structure.UBIQUITIN));
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void matchesBruteForceOverAllAtoms(TestStructures.Structure s) {
        IAtom[] atoms = AtomContainerManipulator.getAtomArray(s.load());
        for (double radius : new double[]{1.5, 3.0, 5.0}) {
            NeighborList nl = new NeighborList(atoms, radius);
            for (int i = 0; i < atoms.length; i++) {
                int[] expected = bruteForce(atoms, i, radius);
                int[] actual = sorted(nl.getNeighbors(i));
                assertArrayEqualsMsg(expected, actual, s + " atom " + i + " radius " + radius);
            }
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void getNeighborsIntoMatchesGetNeighbors(TestStructures.Structure s) {
        IAtom[] atoms = AtomContainerManipulator.getAtomArray(s.load());
        NeighborList nl = new NeighborList(atoms, 3.0);
        IntArrayList buf = new IntArrayList();
        for (int i = 0; i < atoms.length; i++) {
            int[] viaArray = sorted(nl.getNeighbors(i));
            buf.clear();
            nl.getNeighborsInto(i, buf);
            int[] viaBuffer = sorted(buf.toArray());
            assertArrayEqualsMsg(viaArray, viaBuffer, s + " getNeighborsInto disagrees at atom " + i);
            assertEquals(viaArray.length, nl.getNumberOfNeighbors(i),
                    s + " getNumberOfNeighbors disagrees at atom " + i);
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void neighbourRelationIsSymmetricAndExcludesSelf(TestStructures.Structure s) {
        IAtom[] atoms = AtomContainerManipulator.getAtomArray(s.load());
        NeighborList nl = new NeighborList(atoms, 3.0);
        for (int i = 0; i < atoms.length; i++) {
            int[] ni = nl.getNeighbors(i);
            for (int j : ni) {
                assertFalse(j == i, "atom must not be its own neighbour");
                int[] nj = nl.getNeighbors(j);
                assertTrue(contains(nj, i),
                        s + " neighbour relation not symmetric: " + j + " in N(" + i + ") but not vice versa");
            }
        }
    }

    @Test
    void explicitLineFixture() {
        // Five atoms on the x-axis at 0,1,2,3,4. radius 0.75 -> cutoff distance 1.5.
        AtomContainer c = new AtomContainer(5, 0, 0, 0);
        for (int k = 0; k < 5; k++) c.addAtom(new Atom("C", new Point3d(k, 0, 0)));
        IAtom[] atoms = AtomContainerManipulator.getAtomArray(c);
        NeighborList nl = new NeighborList(atoms, 0.75);

        assertArrayEqualsMsg(new int[]{1}, sorted(nl.getNeighbors(0)), "endpoint neighbours");
        assertArrayEqualsMsg(new int[]{1, 3}, sorted(nl.getNeighbors(2)), "middle neighbours");
        assertArrayEqualsMsg(new int[]{3}, sorted(nl.getNeighbors(4)), "endpoint neighbours");
    }

    private static boolean contains(int[] a, int v) {
        for (int x : a) if (x == v) return true;
        return false;
    }

    private static void assertArrayEqualsMsg(int[] expected, int[] actual, String msg) {
        assertEquals(Arrays.toString(expected), Arrays.toString(actual), msg);
    }
}
