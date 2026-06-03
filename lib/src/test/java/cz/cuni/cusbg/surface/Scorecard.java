package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Opt-in accuracy + intrinsic-quality scorecard over the whole {@link SurfaceCatalog}: ./gradlew scorecard
 *
 * <p>Prints one table — it asserts nothing (the per-variant invariants live in {@code SurfaceQualityTest}
 * and {@code FloatNumericalSurfaceTest}); this is the human-readable Pareto view that pairs with the JMH
 * speed numbers. Each variant is evaluated by its {@link SurfaceCatalog.Fidelity} tier, so a future
 * {@code SAMPLING} surface reports area agreement + quality but is never compared point-for-point.
 */
@Tag("scorecard")
class Scorecard {

    private static final double SOLVENT = 1.4;
    private static final int TESS = 2;   // p2rank's operating point; small enough for the O(points²) quality scan
    private static final TestStructures.Structure STRUCTURE = TestStructures.Structure.UBIQUITIN;

    @Test
    void report() {
        IAtomContainer mol = STRUCTURE.load();
        MolecularSurface oracle = new FasterNumericalSurface(mol, SOLVENT, TESS);   // exact area reference

        System.out.printf("%n=== Surface scorecard: %s, solvent=%.1f, tess=%d ===%n",
                STRUCTURE.pdbId, SOLVENT, TESS);
        System.out.printf("(areaRelErr vs FasterNumericalSurface; quality metrics are oracle-free)%n%n");
        System.out.printf("%-30s %-10s %12s %9s %9s %9s %9s%n",
                "variant", "fidelity", "areaRelErr", "dupRatio", "tooClose", "evenCV", "min/mean");
        System.out.println("-".repeat(96));

        for (SurfaceCatalog.Variant v : SurfaceCatalog.ALL) {
            MolecularSurface s;
            try {
                s = v.factory().build(mol, SOLVENT, TESS);
            } catch (RuntimeException e) {
                System.out.printf("%-30s %-10s  (skipped: %s)%n", v.label(), v.fidelity(), e);
                continue;
            }

            String areaErr = switch (v.fidelity()) {
                case REFERENCE -> "  baseline";
                // BIT_EXACT / TOLERANCE / SAMPLING all report total-area agreement vs the exact oracle;
                // only the (here unused) point-set metric would distinguish them, and sampling has none.
                default -> String.format("%.2e", SurfaceAccuracy.compare(s, oracle).areaTotalRelErr());
            };

            SurfaceQuality.Report q = SurfaceQuality.measure(s);
            System.out.printf("%-30s %-10s %12s %9.3f %9.3f %9.3f %9.3f%n",
                    v.label(), v.fidelity(), areaErr,
                    q.duplicateRatio(), q.tooCloseRatio(), q.evennessCV(), q.minMeanNNRatio());
        }
        System.out.println();
    }
}
