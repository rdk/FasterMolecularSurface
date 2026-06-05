package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Phase 6 gate: V24 (`DevSurfaceV24FloatBuild` = float-precision neighbor build + double scan, idea C1)
 * is non-bit-exact (a pair at the cutoff boundary can flip when the {@code d²<sumR²} test is done in
 * float), so before any tess-3 speed claim it must pass the tolerance envelope vs the bit-exact double
 * reference V3, over the corpus at tess 2 AND tess 3. Envelope: total area ≤ 1e-4 rel, per-atom ≤ 2%,
 * point-set symmetric difference ≤ 1e-3. Load-immune. See {@code autoresearch/LOG.md} Phase 6.
 *
 * <p>Run: {@code ./gradlew scorecard --tests '*FloatBuildToleranceGateTest'}
 */
@Tag("scorecard")
class FloatBuildToleranceGateTest {

    private static final double SOLVENT = 1.4;

    @Test
    void gate() {
        for (int tess : new int[]{2, 3}) {
            double maxTotal = 0, maxAtom = 0, maxSym = 0;
            String worst = "";
            for (TestStructures.Structure s : TestStructures.Structure.values()) {
                IAtomContainer mol = s.load();
                MolecularSurface cand = SurfaceCatalog.byId("V24").factory().build(mol, SOLVENT, tess);
                MolecularSurface oracle = SurfaceCatalog.byId("DISTINCT_PACKED_V3").factory().build(mol, SOLVENT, tess);
                SurfaceAccuracy.Report r = SurfaceAccuracy.compare(cand, oracle);
                if (r.areaTotalRelErr() > maxTotal) { maxTotal = r.areaTotalRelErr(); worst = s.name(); }
                maxAtom = Math.max(maxAtom, r.areaAtomMaxRelErr());
                maxSym = Math.max(maxSym, r.pointSetSymDiffFrac());
            }
            boolean pass = maxTotal <= 1e-4 && maxAtom <= 0.02 && maxSym <= 1e-3;
            System.out.printf("%n=== Phase 6 V24 float-build tolerance gate @ tess %d (vs V3, corpus) ===%n", tess);
            System.out.printf("  total area rel err (max): %.3e  (limit 1e-4)%n", maxTotal);
            System.out.printf("  per-atom rel err (max):   %.3e  (limit 2e-2)  worst struct=%s%n", maxAtom, worst);
            System.out.printf("  point-set symdiff frac:   %.3e  (limit 1e-3)%n", maxSym);
            System.out.printf("  GATE: %s%n", pass ? "PASS" : "FAIL");
        }
        System.out.println();
    }
}
