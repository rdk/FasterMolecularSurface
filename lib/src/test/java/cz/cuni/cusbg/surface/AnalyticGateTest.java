package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Phase 4 (Lead B1 — analytic per-atom SASA) cheap gate. Analytic SASA is the converged limit of the
 * dot-sampled area as the tessellation refines, so a finely-tessellated surface is a faithful proxy for
 * the "analytic truth". This measures how far the sampled tess-2 and tess-3 areas are from that limit
 * (total and worst per-atom relative error), over the corpus — bounding what an analytic method could
 * gain on accuracy. Area-only (no point-set keying, so it scales to fine tessellations). Load-immune.
 *
 * <p>Run: {@code ./gradlew scorecard --tests '*AnalyticGateTest'}
 */
@Tag("scorecard")
class AnalyticGateTest {

    private static final double SOLVENT = 1.4;
    private static final int REF = 5;   // reference tessellation (≈analytic; 2562 directions)

    @Test
    void gate() {
        int[] levels = {2, 3, 4};
        // per-level accumulators
        double[] sumAbsTotalErr = new double[levels.length];
        double[] sumRefTotal = new double[levels.length];
        double[] maxAtomRel = new double[levels.length];
        double[] sumAtomRel = new double[levels.length];
        long[] atomCount = new long[levels.length];

        long t0 = System.currentTimeMillis();
        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            IAtomContainer mol = s.load();
            double[] refAreas = new FasterNumericalSurface(mol, SOLVENT, REF).getAllSurfaceAreas();
            double refTotal = 0; for (double a : refAreas) refTotal += a;

            for (int li = 0; li < levels.length; li++) {
                double[] ca = new FasterNumericalSurface(mol, SOLVENT, levels[li]).getAllSurfaceAreas();
                double ct = 0; for (double a : ca) ct += a;
                sumAbsTotalErr[li] += Math.abs(ct - refTotal);
                sumRefTotal[li] += refTotal;
                int n = Math.min(ca.length, refAreas.length);
                for (int i = 0; i < n; i++) {
                    double denom = Math.max(1.0, refAreas[i]);
                    double rel = Math.abs(ca[i] - refAreas[i]) / denom;
                    maxAtomRel[li] = Math.max(maxAtomRel[li], rel);
                    sumAtomRel[li] += rel;
                    atomCount[li]++;
                }
            }
        }
        long elapsed = System.currentTimeMillis() - t0;

        System.out.printf("%n=== Phase 4 B1 analytic gate: sampled tess vs tess-%d (≈analytic) over corpus ===%n", REF);
        System.out.printf("p2rank tolerance: total ≤ 1e-4 rel, per-atom ≤ 2%%   (elapsed %d ms)%n", elapsed);
        for (int li = 0; li < levels.length; li++) {
            double totalRel = sumRefTotal[li] == 0 ? 0 : sumAbsTotalErr[li] / sumRefTotal[li];
            double meanAtomRel = atomCount[li] == 0 ? 0 : sumAtomRel[li] / atomCount[li];
            System.out.printf("  tess %d: total area rel err %.3e | per-atom rel err mean %.3e max %.3e%n",
                    levels[li], totalRel, meanAtomRel, maxAtomRel[li]);
        }
        System.out.println();
    }
}
