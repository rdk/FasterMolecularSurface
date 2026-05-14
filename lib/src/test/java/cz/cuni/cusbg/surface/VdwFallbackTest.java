package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.silent.Atom;
import org.openscience.cdk.silent.AtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import javax.vecmath.Point3d;
import java.io.InputStream;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Regression tests for the VdW-radius fallback added to {@link FasterNumericalSurface}.
 *
 * <p>Without the fallback, surface computation on atoms whose CDK {@code Elements} enum
 * entry has {@code null} for VdW radius (Co, Ni, Cu, Rh, Os, Ir, plus radioactive /
 * synthetic elements) throws NPE in {@code init()}. See {@code radii-vdw.txt} bundled
 * with this jar and the writeup in p2rank-private/local/cdk-vdw-radius-gap.md.
 */
class VdwFallbackTest {

    private static IAtomContainer singleAtomContainer(String symbol, double x, double y, double z) {
        AtomContainer c = new AtomContainer(1, 0, 0, 0);
        c.addAtom(new Atom(symbol, new Point3d(x, y, z)));
        return c;
    }

    private static IAtomContainer twoAtomContainer(String symbol1, String symbol2) {
        AtomContainer c = new AtomContainer(2, 0, 0, 0);
        c.addAtom(new Atom(symbol1, new Point3d(0.0, 0.0, 0.0)));
        c.addAtom(new Atom(symbol2, new Point3d(3.0, 0.0, 0.0)));
        return c;
    }

    // --- Single-atom regressions: must not NPE ---

    @Test
    public void cobaltAtomDoesNotCrash() {
        FasterNumericalSurface surface = new FasterNumericalSurface(singleAtomContainer("Co", 0, 0, 0));
        assertTrue(surface.getTotalSurfaceArea() > 0,
                "Cobalt atom should produce a positive surface area, got " + surface.getTotalSurfaceArea());
    }

    @Test
    public void nickelAtomDoesNotCrash() {
        FasterNumericalSurface surface = new FasterNumericalSurface(singleAtomContainer("Ni", 0, 0, 0));
        assertTrue(surface.getTotalSurfaceArea() > 0);
    }

    @Test
    public void copperAtomDoesNotCrash() {
        FasterNumericalSurface surface = new FasterNumericalSurface(singleAtomContainer("Cu", 0, 0, 0));
        assertTrue(surface.getTotalSurfaceArea() > 0);
    }

    @Test
    public void iridiumAtomDoesNotCrash() {
        FasterNumericalSurface surface = new FasterNumericalSurface(singleAtomContainer("Ir", 0, 0, 0));
        assertTrue(surface.getTotalSurfaceArea() > 0);
    }

    // --- Fallback gives the expected radius (within tolerance) ---

    @Test
    public void cobaltSurfaceMatchesExpectedRadius() {
        // VdW from CDK's radii-vdw.txt for Co (Z=27) is 2.0 A.
        // With solvent 1.4 A, total radius = 3.4 A. Area = 4 * pi * r^2 ~= 145.27 A^2.
        // Allow 1% tolerance for tessellation discretisation.
        FasterNumericalSurface surface = new FasterNumericalSurface(singleAtomContainer("Co", 0, 0, 0));
        double expected = 4.0 * Math.PI * 3.4 * 3.4;
        double actual = surface.getTotalSurfaceArea();
        assertEquals(expected, actual, expected * 0.01,
                "Co single-atom surface should be ~4*pi*(2.0+1.4)^2 A^2, got " + actual);
    }

    @Test
    public void cobaltAndKryptonSurfacesAreClose() {
        // Krypton has CDK enum VdW = 2.02 A. Cobalt falls back to 2.0 A.
        // Single-atom surfaces should differ by < 2.5% (proportional to r^2).
        double co = new FasterNumericalSurface(singleAtomContainer("Co", 0, 0, 0)).getTotalSurfaceArea();
        double kr = new FasterNumericalSurface(singleAtomContainer("Kr", 0, 0, 0)).getTotalSurfaceArea();
        double rel = Math.abs(co - kr) / kr;
        assertTrue(rel < 0.025,
                "Co and Kr single-atom surfaces should be within 2.5% (got Co=" + co + ", Kr=" + kr + ", rel=" + rel + ")");
    }

    // --- Two-atom mixed: cobalt in proximity to carbon should still compute cleanly ---

    @Test
    public void carbonPlusCobaltDoesNotCrash() {
        FasterNumericalSurface surface = new FasterNumericalSurface(twoAtomContainer("C", "Co"));
        assertTrue(surface.getTotalSurfaceArea() > 0);
        // Should be less than the sum of isolated atoms (they occlude each other).
        double cIsolated = new FasterNumericalSurface(singleAtomContainer("C", 0, 0, 0)).getTotalSurfaceArea();
        double coIsolated = new FasterNumericalSurface(singleAtomContainer("Co", 0, 0, 0)).getTotalSurfaceArea();
        assertTrue(surface.getTotalSurfaceArea() < cIsolated + coIsolated,
                "Mixed C+Co surface should show some occlusion vs. isolated atoms");
    }

    // --- Unknown / synthetic element: must fall back to 2.0 A default ---

    @Test
    public void syntheticElementDoesNotCrash() {
        // Tennessine (Ts, Z=117) is null in CDK's enum and likely absent from radii-vdw.txt.
        // Should hit the 2.0 A default fallback. Just verify no NPE.
        FasterNumericalSurface surface = new FasterNumericalSurface(singleAtomContainer("Ts", 0, 0, 0));
        assertTrue(surface.getTotalSurfaceArea() > 0);
    }

    // --- Real-world PDB: cobalamin (vitamin B12) structure 3CI3 ---

    @Test
    public void cobalaminStructureSurfaceComputes() throws Exception {
        // 3CI3 is methylmalonyl-CoA mutase. Without the VdW fallback, computing the
        // surface over ALL atoms (including the cobalamin cobalt) throws NPE.
        IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
        IChemFile chemFile;
        try (InputStream in = getClass().getResourceAsStream("/data/pdb/3CI3.pdb");
             PDBReader pdbr = new PDBReader(in)) {
            chemFile = pdbr.read(bldr.newInstance(IChemFile.class));
        }
        IAtomContainer mol = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);

        // Sanity: structure must contain at least one cobalt atom (otherwise the test
        // doesn't exercise the fallback path).
        boolean hasCobalt = false;
        for (IAtom a : mol.atoms()) {
            if ("Co".equals(a.getSymbol())) {
                hasCobalt = true;
                break;
            }
        }
        assertTrue(hasCobalt, "Sanity: 3CI3 must contain at least one Co atom");

        FasterNumericalSurface surface = new FasterNumericalSurface(mol);
        assertTrue(surface.getTotalSurfaceArea() > 0,
                "3CI3 (cobalamin-bearing) surface area should be positive");
        assertTrue(surface.getAllSurfacePoints().length > 0,
                "3CI3 should produce at least one surface point");
    }
}
