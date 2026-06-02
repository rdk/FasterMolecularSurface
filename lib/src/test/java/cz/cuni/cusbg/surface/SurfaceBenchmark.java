package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.Arrays;
import java.util.function.DoubleSupplier;

/**
 * Opt-in performance harness for comparing surface implementations.
 *
 * <p>Not part of the normal test run (tagged {@code benchmark}, excluded by the {@code test}
 * task). Run with:
 * <pre>./gradlew benchmark</pre>
 *
 * <p>It reports median wall-clock time per structure for {@link FasterNumericalSurface} against
 * CDK's reference {@code NumericalSurface}, giving a speedup factor. When a new optimized variant
 * is added, copy the {@code FasterNumericalSurface} row to benchmark it on the same corpus and
 * parameters. This is a coarse timing harness (median of a few runs after warm-up), not JMH — it
 * is meant for relative comparison between variants on the same machine, not absolute numbers.
 */
@Tag("benchmark")
class SurfaceBenchmark {

    private static final int WARMUP = 2;
    private static final int MEASURE = 3;

    /** Accumulates results so the JIT cannot eliminate the timed work. */
    private static volatile double blackhole;

    /**
     * Times CDK's reference {@code NumericalSurface} and the optimized variants
     * ({@link FasterNumericalSurface}, {@link SoaNumericalSurface}, {@link GridSoaNumericalSurface})
     * on the same corpus, at both the library-default tessellation and p2rank's lower one. Add a
     * column when a new variant lands.
     */
    private static final int[] TESS_LEVELS = {4, 3, 2};   // 4 = library default, 2 = p2rank's operating point

    @Test
    void benchmarkVariants() {
        System.out.printf("vectorized scan active (jdk.incubator.vector available): %s%n",
                VectorizedSymmetricHintedGridSoaNumericalSurface.isVectorized());

        // per-level detail tables, collecting the mean speed-vs-CDK of each generator
        double[][] meanVsCdk = new double[TESS_LEVELS.length][]; // [level] -> {faster, soa, grid}
        for (int i = 0; i < TESS_LEVELS.length; i++) {
            meanVsCdk[i] = runVariantTable(1.4, TESS_LEVELS[i]);
        }

        // summary: all four generators, mean speedup vs CDK (over the CDK-runnable structures)
        System.out.printf("%n=== Summary: mean speedup vs CDK NumericalSurface (corpus) ===%n");
        System.out.printf("%-26s", "generator");
        for (int t : TESS_LEVELS) System.out.printf(" %9s", "tess" + t);
        System.out.println();
        System.out.println("------------------------------------------------------");
        printSummaryRow("CDK NumericalSurface", meanVsCdk, -1);
        printSummaryRow("FasterNumericalSurface", meanVsCdk, 0);
        printSummaryRow("SoaNumericalSurface", meanVsCdk, 1);
        printSummaryRow("GridSoaNumericalSurface", meanVsCdk, 2);
        printSummaryRow("OrderedGridSoaNumericalSurface", meanVsCdk, 3);
        printSummaryRow("HintedGridSoaNumericalSurface", meanVsCdk, 4);
        printSummaryRow("SymmetricHintedGridSoaNumericalSurface", meanVsCdk, 5);
        printSummaryRow("LowAllocSymmetricHintedGridSoaNumericalSurface", meanVsCdk, 6);
        printSummaryRow("VectorizedSymmetricHintedGridSoaNumericalSurface", meanVsCdk, 7);
        System.out.println();
    }

    private static void printSummaryRow(String name, double[][] meanVsCdk, int gen) {
        System.out.printf("%-26s", name);
        for (double[] level : meanVsCdk) {
            double v = gen < 0 ? 1.0 : level[gen];   // CDK is the 1.00x baseline
            System.out.printf(" %8.2fx", v);
        }
        System.out.println();
    }

    /** Prints the per-structure table at one tessellation level; returns mean {faster,soa,grid,ord,hint,sym,symLA,vec} speedup vs CDK. */
    private double[] runVariantTable(double solvent, int tess) {
        System.out.printf("%n=== CDK vs variants @ solventRadius=%.1f, tessLevel=%d ===%n", solvent, tess);
        System.out.printf("%-12s %7s %8s %8s %8s %8s %8s %8s %8s %8s %8s %9s %9s%n",
                "structure", "atoms", "cdk(ms)", "faster", "soa", "grid", "ord", "hint", "sym", "symLA", "vec", "sym/cdk", "vec/cdk");
        System.out.println("------------------------------------------------------------------------------------------------------------------------------");
        double sumFst = 0, sumSoa = 0, sumGrid = 0, sumOrd = 0, sumHint = 0, sumSym = 0, sumSymLA = 0, sumVec = 0; int cdkCount = 0;
        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            IAtomContainer mol = s.load();

            // CDK's reference NumericalSurface NPEs on cobalt (3CI3); skip its column there.
            Double cdk = null;
            if (s != TestStructures.Structure.COBALAMIN) {
                cdk = medianMillis(() -> {
                    NumericalSurface ns = new NumericalSurface(mol, solvent, tess);
                    ns.calculateSurface();
                    return ns.getTotalSurfaceArea();
                });
            }
            double faster = medianMillis(() -> new FasterNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double soa    = medianMillis(() -> new SoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double grid   = medianMillis(() -> new GridSoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double ord    = medianMillis(() -> new OrderedGridSoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double hint   = medianMillis(() -> new HintedGridSoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double sym    = medianMillis(() -> new SymmetricHintedGridSoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double symLA  = medianMillis(() -> new LowAllocSymmetricHintedGridSoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());
            double vec    = medianMillis(() -> new VectorizedSymmetricHintedGridSoaNumericalSurface(mol, solvent, tess).getTotalSurfaceArea());

            String cdkStr = cdk == null ? "n/a(Co)" : String.format("%.1f", cdk);
            String sy = "-", ve = "-";
            if (cdk != null) {
                sy = String.format("%.2fx", cdk / sym);
                ve = String.format("%.2fx", cdk / vec);
                sumFst += cdk / faster; sumSoa += cdk / soa; sumGrid += cdk / grid; sumOrd += cdk / ord;
                sumHint += cdk / hint; sumSym += cdk / sym; sumSymLA += cdk / symLA; sumVec += cdk / vec; cdkCount++;
            }
            System.out.printf("%-12s %7d %8s %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %9s %9s%n",
                    s.pdbId, s.atomCount, cdkStr, faster, soa, grid, ord, hint, sym, symLA, vec, sy, ve);
        }
        System.out.println();
        int c = Math.max(1, cdkCount);
        return new double[]{ sumFst / c, sumSoa / c, sumGrid / c, sumOrd / c, sumHint / c, sumSym / c, sumSymLA / c, sumVec / c };
    }

    /** Median wall-clock time in milliseconds over MEASURE runs after WARMUP warm-up runs. */
    private static double medianMillis(DoubleSupplier work) {
        for (int i = 0; i < WARMUP; i++) blackhole += work.getAsDouble();
        double[] times = new double[MEASURE];
        for (int i = 0; i < MEASURE; i++) {
            long t0 = System.nanoTime();
            blackhole += work.getAsDouble();
            long t1 = System.nanoTime();
            times[i] = (t1 - t0) / 1_000_000.0;
        }
        Arrays.sort(times);
        return times[times.length / 2];
    }
}
