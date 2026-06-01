package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.AtomContainer;

import javax.vecmath.Point3d;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Tests for the public API contract of {@link FasterNumericalSurface}: index-bounds behaviour,
 * constructor parameter handling, immutability of returned collections, and input validation.
 *
 * <p>Uses the small 1CRN structure so these run fast.
 */
class ApiContractTest {

    private static FasterNumericalSurface crambin() {
        return new FasterNumericalSurface(TestStructures.Structure.CRAMBIN.load());
    }

    // --- Index bounds: must raise CDKException, never AIOOBE ----------------------------

    @Test
    void getSurfaceAreaRejectsIndexAtSize() {
        FasterNumericalSurface s = crambin();
        assertThrows(CDKException.class, () -> s.getSurfaceArea(s.getAllSurfaceAreas().length));
    }

    @Test
    void getSurfaceAreaRejectsNegativeIndex() {
        // Regression: a negative index used to throw ArrayIndexOutOfBoundsException, not CDKException.
        FasterNumericalSurface s = crambin();
        assertThrows(CDKException.class, () -> s.getSurfaceArea(-1));
    }

    @Test
    void getSurfacePointsRejectsIndexAtSize() {
        FasterNumericalSurface s = crambin();
        assertThrows(CDKException.class, () -> s.getSurfacePoints(s.getAllSurfaceAreas().length));
    }

    @Test
    void getSurfacePointsRejectsNegativeIndex() {
        // Regression: negative index must map to CDKException, consistent with the upper bound.
        FasterNumericalSurface s = crambin();
        assertThrows(CDKException.class, () -> s.getSurfacePoints(-1));
    }

    @Test
    void validIndicesDoNotThrow() {
        FasterNumericalSurface s = crambin();
        int last = s.getAllSurfaceAreas().length - 1;
        assertDoesNotThrow(() -> s.getSurfaceArea(0));
        assertDoesNotThrow(() -> s.getSurfaceArea(last));
        assertDoesNotThrow(() -> s.getSurfacePoints(0));
        assertDoesNotThrow(() -> s.getSurfacePoints(last));
    }

    // --- Returned collections are not back-doors into internal state --------------------

    @Test
    void atomSurfaceMapValuesAreUnmodifiable() {
        Map<IAtom, List<Point3d>> map = crambin().getAtomSurfaceMap();
        List<Point3d> anyExposed = map.values().iterator().next();
        assertThrows(UnsupportedOperationException.class,
                () -> anyExposed.add(new Point3d(0, 0, 0)));
    }

    // --- Constructor parameters are honoured --------------------------------------------

    @Test
    void zeroSolventRadiusGivesVanDerWaalsSurface() {
        // solventRadius = 0 must be accepted and yield a valid (positive) Van der Waals surface.
        FasterNumericalSurface s = new FasterNumericalSurface(
                TestStructures.Structure.CRAMBIN.load(), 0.0, 4);
        assertTrue(s.getTotalSurfaceArea() > 0);
        assertTrue(s.getAllSurfacePoints().length > 0);
    }

    @Test
    void tessellationLevelChangesPointDensity() {
        IAtomContainer mol = TestStructures.Structure.CRAMBIN.load();
        int coarse = new FasterNumericalSurface(mol, 1.4, 2).getAllSurfacePoints().length;
        int fine = new FasterNumericalSurface(TestStructures.Structure.CRAMBIN.load(), 1.4, 4)
                .getAllSurfacePoints().length;
        assertTrue(fine > coarse,
                "Higher tessellation level should produce more surface points (coarse=" + coarse
                        + ", fine=" + fine + ")");
    }

    @Test
    void isolatedAtomAreaIsCloseToAnalyticSphere() {
        // A single isolated carbon has no occluders, so its area must approximate 4*pi*r^2
        // (r = VdW 1.7 + solvent 1.4) to within the tessellation discretisation (~1% at level 4).
        double r = 1.7 + 1.4;
        double analytic = 4.0 * Math.PI * r * r;
        double actual = oneAtomArea("C", 4);
        assertEquals(analytic, actual, analytic * 0.01,
                "Isolated-atom area should be within 1% of 4*pi*r^2");
    }

    private static double oneAtomArea(String symbol, int tessLevel) {
        AtomContainer c = new AtomContainer(1, 0, 0, 0);
        c.addAtom(new Atom(symbol, new Point3d(0, 0, 0)));
        return new FasterNumericalSurface(c, 1.4, tessLevel).getTotalSurfaceArea();
    }

    // --- Input validation ----------------------------------------------------------------

    @Test
    void atomWithoutCoordinatesIsRejected() {
        AtomContainer c = new AtomContainer(2, 0, 0, 0);
        c.addAtom(new Atom("C", new Point3d(0, 0, 0)));
        c.addAtom(new Atom("C")); // no 3D coordinate
        assertThrows(IllegalArgumentException.class, () -> new FasterNumericalSurface(c));
    }

    @Test
    void twoFarApartAtomsAreNearlyTheSumOfIsolatedAtoms() {
        // Sanity on the area accounting: atoms far beyond any interaction range do not occlude.
        double iso = oneAtomArea("C", 4);
        AtomContainer c = new AtomContainer(2, 0, 0, 0);
        c.addAtom(new Atom("C", new Point3d(0, 0, 0)));
        c.addAtom(new Atom("C", new Point3d(1000, 0, 0)));
        double pair = new FasterNumericalSurface(c, 1.4, 4).getTotalSurfaceArea();
        assertEquals(2 * iso, pair, 1e-9 * 2 * iso,
                "Two far-apart identical atoms should equal twice an isolated atom");
        assertNotEquals(0.0, pair);
    }
}
