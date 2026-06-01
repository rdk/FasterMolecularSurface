package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Establishes {@link FasterNumericalSurface} as a faithful optimization of CDK's reference
 * {@code org.openscience.cdk.geometry.surface.NumericalSurface}.
 *
 * <p>This is the <em>oracle</em>: it proves the fast implementation reproduces CDK's results
 * (within floating-point tolerance) across a wide range of structures. The golden-value test
 * then pins those results so future variants can be checked cheaply without re-running CDK.
 *
 * <p>3CI3 (cobalamin) is deliberately excluded — CDK's reference throws NPE on its cobalt atom
 * (the VdW-radius gap that motivated {@link FasterNumericalSurface}'s fallback), so there is no
 * reference value to compare against. That structure is covered by {@link VdwFallbackTest}.
 */
class CdkEquivalenceTest {

    /** All structures except the ones CDK's reference implementation cannot process. */
    static Stream<Arguments> cdkCompatibleStructures() {
        return Stream.of(TestStructures.Structure.values())
                .filter(s -> s != TestStructures.Structure.COBALAMIN)
                .map(Arguments::of);
    }

    private static NumericalSurface reference(IAtomContainer mol, double solventRadius, int tessLevel) {
        NumericalSurface ns = new NumericalSurface(mol, solventRadius, tessLevel);
        ns.calculateSurface();
        return ns;
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("cdkCompatibleStructures")
    void totalAreaMatchesCdk(TestStructures.Structure s) {
        IAtomContainer mol = s.load();
        double fast = new FasterNumericalSurface(mol).getTotalSurfaceArea();
        double ref = reference(s.load(), 1.4, 4).getTotalSurfaceArea();
        assertEquals(ref, fast, 1e-6 * ref,
                () -> s + " total surface area must match CDK reference (ref=" + ref + ", fast=" + fast + ")");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("cdkCompatibleStructures")
    void perAtomAreasMatchCdk(TestStructures.Structure s) {
        IAtomContainer mol = s.load();
        double[] fast = new FasterNumericalSurface(mol).getAllSurfaceAreas();
        double[] ref = reference(s.load(), 1.4, 4).getAllSurfaceAreas();
        assertEquals(ref.length, fast.length, () -> s + " atom count mismatch vs CDK");
        for (int i = 0; i < ref.length; i++) {
            final int idx = i;
            assertEquals(ref[i], fast[i], 1e-6 * Math.max(1.0, ref[i]),
                    () -> s + " per-atom area mismatch at atom " + idx);
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("cdkCompatibleStructures")
    void surfacePointCountMatchesCdk(TestStructures.Structure s) {
        IAtomContainer mol = s.load();
        int fast = new FasterNumericalSurface(mol).getAllSurfacePoints().length;
        int ref = reference(s.load(), 1.4, 4).getAllSurfacePoints().length;
        assertEquals(ref, fast, () -> s + " surface point count must match CDK reference");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("cdkCompatibleStructures")
    void matchesCdkAtVanDerWaalsRadius(TestStructures.Structure s) {
        // solvent radius 0 (Van der Waals surface), coarser tessellation for speed.
        IAtomContainer mol = s.load();
        double fast = new FasterNumericalSurface(mol, 0.0, 3).getTotalSurfaceArea();
        double ref = reference(s.load(), 0.0, 3).getTotalSurfaceArea();
        assertTrue(ref > 0, () -> s + " reference VdW area should be positive");
        assertEquals(ref, fast, 1e-6 * ref,
                () -> s + " VdW-surface total area must match CDK reference");
    }
}
