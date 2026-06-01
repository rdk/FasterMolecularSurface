package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import javax.vecmath.Point3d;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Implementation-independent behavioural contract for {@link MolecularSurface}.
 *
 * <p>Every variant of the surface algorithm must satisfy these properties. A new variant
 * gets the entire battery for free by subclassing this test and overriding the two factory
 * methods. The properties checked here are internal-consistency invariants that hold for any
 * correct implementation, independent of the exact numbers produced (those are pinned
 * separately by the golden-value test and validated against CDK by the equivalence test).
 */
abstract class AbstractMolecularSurfaceContractTest {

    /** Default parameters (solvent radius 1.4 A, tessellation level 4). */
    protected abstract MolecularSurface create(IAtomContainer mol);

    /** User-specified parameters. */
    protected abstract MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel);

    static java.util.stream.Stream<org.junit.jupiter.params.provider.Arguments> structures() {
        return TestStructures.all();
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void totalAreaIsPositive(TestStructures.Structure s) {
        MolecularSurface surface = create(s.load());
        assertTrue(surface.getTotalSurfaceArea() > 0,
                () -> s + " should have positive total surface area");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void perAtomAreasAreNonNegativeAndSumToTotal(TestStructures.Structure s) {
        MolecularSurface surface = create(s.load());
        double[] areas = surface.getAllSurfaceAreas();
        assertEquals(s.atomCount, areas.length, () -> s + " area array length must equal atom count");
        double sum = 0;
        for (double a : areas) {
            assertTrue(a >= 0, () -> s + " per-atom area must be non-negative");
            sum += a;
        }
        // getTotalSurfaceArea sums the same array, so this is exact.
        assertEquals(surface.getTotalSurfaceArea(), sum, 1e-9 * Math.max(1, Math.abs(sum)),
                () -> s + " sum of per-atom areas must equal total area");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void getSurfaceAreaMatchesAllSurfaceAreas(TestStructures.Structure s) throws Exception {
        MolecularSurface surface = create(s.load());
        double[] areas = surface.getAllSurfaceAreas();
        for (int i = 0; i < areas.length; i++) {
            assertEquals(areas[i], surface.getSurfaceArea(i), 0.0,
                    "getSurfaceArea(" + i + ") must match getAllSurfaceAreas()[" + i + "]");
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void allSurfacePointsEqualsSumOfPerAtomPoints(TestStructures.Structure s) throws Exception {
        MolecularSurface surface = create(s.load());
        int total = surface.getAllSurfacePoints().length;

        int sumPerAtom = 0;
        for (int i = 0; i < s.atomCount; i++) {
            sumPerAtom += surface.getSurfacePoints(i).length;
        }
        assertEquals(total, sumPerAtom,
                () -> s + " getAllSurfacePoints() must equal the sum of per-atom point counts");

        int sumMap = 0;
        for (List<Point3d> pts : surface.getAtomSurfaceMap().values()) {
            sumMap += pts.size();
        }
        assertEquals(total, sumMap,
                () -> s + " atom-surface map point counts must equal getAllSurfacePoints()");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void atomSurfaceMapContainsExactlyTheExposedAtoms(TestStructures.Structure s) throws Exception {
        MolecularSurface surface = create(s.load());
        Map<IAtom, List<Point3d>> map = surface.getAtomSurfaceMap();

        // No mapped atom may have an empty point list (buried atoms are excluded entirely).
        for (List<Point3d> pts : map.values()) {
            assertTrue(!pts.isEmpty(), () -> s + " mapped atoms must have at least one surface point");
        }
        assertTrue(map.size() <= s.atomCount, () -> s + " map cannot have more atoms than the molecule");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void surfacePointsLieOnTheExpandedSphereOfTheirAtom(TestStructures.Structure s) throws Exception {
        // Every surface point of an atom must sit at distance (vdw + solvent) from that atom's
        // centre. We don't know the per-element radius here, but for each atom the distance must
        // be constant across all its points (the tessellation sphere is rigidly scaled+translated).
        IAtomContainer mol = s.load();
        MolecularSurface surface = create(mol);
        Map<IAtom, List<Point3d>> map = surface.getAtomSurfaceMap();
        for (Map.Entry<IAtom, List<Point3d>> e : map.entrySet()) {
            Point3d centre = e.getKey().getPoint3d();
            List<Point3d> pts = e.getValue();
            double r0 = centre.distance(pts.get(0));
            assertTrue(r0 > 0, () -> s + " surface radius must be positive");
            for (Point3d p : pts) {
                assertEquals(r0, centre.distance(p), 1e-6 * r0,
                        () -> s + " all surface points of an atom must be equidistant from its centre");
            }
        }
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void computationIsDeterministic(TestStructures.Structure s) {
        // Two independent computations of the same input must agree exactly.
        IAtomContainer mol = s.load();
        double a1 = create(mol).getTotalSurfaceArea();
        double a2 = create(s.load()).getTotalSurfaceArea();
        assertEquals(a1, a2, 0.0, () -> s + " total area must be deterministic across runs");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void customParametersProduceAPositiveDistinctSurface(TestStructures.Structure s) {
        // The custom-parameter constructor must work and a different solvent radius must change
        // the result. (The direction of change is not monotonic for packed molecules: a larger
        // probe expands every atom sphere but also occludes more of each neighbour, so total SASA
        // may rise or fall - hence we only assert "positive and different", not an ordering.)
        IAtomContainer mol = s.load();
        double vdw = create(mol, 0.0, 4).getTotalSurfaceArea();
        double sas = create(mol, 1.4, 4).getTotalSurfaceArea();
        assertTrue(vdw > 0, () -> s + " Van der Waals area must be positive");
        assertTrue(sas > 0, () -> s + " solvent-accessible area must be positive");
        assertTrue(Math.abs(sas - vdw) > 1e-6 * vdw,
                () -> s + " changing the solvent radius must change the surface area");
    }

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void translationInvariance(TestStructures.Structure s) {
        // Rigidly translating the whole molecule must not change any surface area.
        IAtomContainer mol = s.load();
        double before = create(mol).getTotalSurfaceArea();

        IAtomContainer shifted = s.load();
        for (IAtom a : shifted.atoms()) {
            Point3d p = a.getPoint3d();
            a.setPoint3d(new Point3d(p.x + 137.0, p.y - 71.5, p.z + 42.25));
        }
        double after = create(shifted).getTotalSurfaceArea();
        assertEquals(before, after, 1e-6 * before,
                () -> s + " total area must be invariant under rigid translation");
    }

    // --- Golden values: pinned baseline every variant must reproduce ---------------------

    @ParameterizedTest(name = "{0}")
    @MethodSource("structures")
    void matchesGoldenValues(TestStructures.Structure s) throws Exception {
        GoldenValues.Row g = GoldenValues.get(s.pdbId);
        MolecularSurface surface = create(s.load());

        assertEquals(g.atomCount, surface.getAllSurfaceAreas().length,
                () -> s + " atom count drifted from golden baseline");
        assertEquals(g.surfaceAtomCount, surface.getAtomSurfaceMap().size(),
                () -> s + " exposed-atom count drifted from golden baseline");
        assertEquals(g.totalPoints, surface.getAllSurfacePoints().length,
                () -> s + " total surface-point count drifted from golden baseline");
        assertEquals(g.totalArea, surface.getTotalSurfaceArea(), 1e-6 * g.totalArea,
                () -> s + " total surface area drifted from golden baseline (golden=" + g.totalArea
                        + ", got=" + surface.getTotalSurfaceArea() + ")");
    }
}
