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
     * ({@link FasterNumericalSurface}, {@link DevSurfaceV1Soa}, {@link DevSurfaceV2Grid})
     * on the same corpus, at both the library-default tessellation and p2rank's lower one. Add a
     * column when a new variant lands.
     */
    private static final int[] TESS_LEVELS = {4, 3, 2};   // 4 = library default, 2 = p2rank's operating point

    @Test
    void benchmarkVariants() {
        System.out.printf("vectorized scan active (jdk.incubator.vector available): %s%n",
                DevSurfaceV7Simd.isVectorized());

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
        printSummaryRow("DevSurfaceV1Soa", meanVsCdk, 1);
        printSummaryRow("DevSurfaceV2Grid", meanVsCdk, 2);
        printSummaryRow("DevSurfaceV3Sorted", meanVsCdk, 3);
        printSummaryRow("DevSurfaceV4Hinted", meanVsCdk, 4);
        printSummaryRow("DevSurfaceV5Symmetric", meanVsCdk, 5);
        printSummaryRow("DevSurfaceV6LowAlloc", meanVsCdk, 6);
        printSummaryRow("DevSurfaceV7Simd", meanVsCdk, 7);
        printSummaryRow("DevSurfaceV8Pruned", meanVsCdk, 8);
        printSummaryRow("DevSurfaceV9Dedup", meanVsCdk, 9);
        printSummaryRow("DevSurfaceV10CachedMap", meanVsCdk, 10);
        printSummaryRow("DevSurfaceV11CachedTess", meanVsCdk, 11);
        printSummaryRow("DevSurfaceV12Flat", meanVsCdk, 12);
        printSummaryRow("DevSurfaceV13Arena", meanVsCdk, 13);
        printSummaryRow("DevSurfaceV14ArenaFlat", meanVsCdk, 14);
        printSummaryRow("DevSurfaceV15LeanNbr", meanVsCdk, 15);
        printSummaryRow("DevSurfaceV16DirectNbr", meanVsCdk, 16);
        printSummaryRow("DevSurfaceV17PackedNbr", meanVsCdk, 17);
        printSummaryRow("DevSurfaceV18SortedCoords", meanVsCdk, 18);
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
                "structure", "atoms", "cdk(ms)", "faster", "soa", "grid", "ord", "hint", "sym", "tcd", "flat", "tcd/cdk", "flat/cdk");
        System.out.println("------------------------------------------------------------------------------------------------------------------------------");
        double sumFst = 0, sumSoa = 0, sumGrid = 0, sumOrd = 0, sumHint = 0, sumSym = 0, sumSymLA = 0, sumVec = 0, sumPrn = 0, sumDed = 0, sumGded = 0, sumTcd = 0, sumFlat = 0, sumAtcd = 0, sumAflat = 0, sumLean = 0, sumDirect = 0, sumPacked = 0, sumSorted = 0; int cdkCount = 0;
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
            double soa    = medianMillis(() -> new DevSurfaceV1Soa(mol, solvent, tess).getTotalSurfaceArea());
            double grid   = medianMillis(() -> new DevSurfaceV2Grid(mol, solvent, tess).getTotalSurfaceArea());
            double ord    = medianMillis(() -> new DevSurfaceV3Sorted(mol, solvent, tess).getTotalSurfaceArea());
            double hint   = medianMillis(() -> new DevSurfaceV4Hinted(mol, solvent, tess).getTotalSurfaceArea());
            double sym    = medianMillis(() -> new DevSurfaceV5Symmetric(mol, solvent, tess).getTotalSurfaceArea());
            double symLA  = medianMillis(() -> new DevSurfaceV6LowAlloc(mol, solvent, tess).getTotalSurfaceArea());
            double vec    = medianMillis(() -> new DevSurfaceV7Simd(mol, solvent, tess).getTotalSurfaceArea());
            double prn    = medianMillis(() -> new DevSurfaceV8Pruned(mol, solvent, tess).getTotalSurfaceArea());
            double ded    = medianMillis(() -> new DevSurfaceV9Dedup(mol, solvent, tess).getTotalSurfaceArea());
            double gded   = medianMillis(() -> new DevSurfaceV10CachedMap(mol, solvent, tess).getTotalSurfaceArea());
            double tcd    = medianMillis(() -> new DevSurfaceV11CachedTess(mol, solvent, tess).getTotalSurfaceArea());
            double flat   = medianMillis(() -> new DevSurfaceV12Flat(mol, solvent, tess).getTotalSurfaceArea());
            double atcd   = medianMillis(() -> new DevSurfaceV13Arena(mol, solvent, tess).getTotalSurfaceArea());
            double aflat  = medianMillis(() -> new DevSurfaceV14ArenaFlat(mol, solvent, tess).getTotalSurfaceArea());
            double lean   = medianMillis(() -> new DevSurfaceV15LeanNbr(mol, solvent, tess).getTotalSurfaceArea());
            double direct = medianMillis(() -> new DevSurfaceV16DirectNbr(mol, solvent, tess).getTotalSurfaceArea());
            double packed = medianMillis(() -> new DevSurfaceV17PackedNbr(mol, solvent, tess).getTotalSurfaceArea());
            double sorted = medianMillis(() -> new DevSurfaceV18SortedCoords(mol, solvent, tess).getTotalSurfaceArea());

            String cdkStr = cdk == null ? "n/a(Co)" : String.format("%.1f", cdk);
            String tc = "-", fl = "-";
            if (cdk != null) {
                tc = String.format("%.2fx", cdk / tcd);
                fl = String.format("%.2fx", cdk / flat);
                sumFst += cdk / faster; sumSoa += cdk / soa; sumGrid += cdk / grid; sumOrd += cdk / ord;
                sumHint += cdk / hint; sumSym += cdk / sym; sumSymLA += cdk / symLA; sumVec += cdk / vec;
                sumPrn += cdk / prn; sumDed += cdk / ded; sumGded += cdk / gded; sumTcd += cdk / tcd; sumFlat += cdk / flat;
                sumAtcd += cdk / atcd; sumAflat += cdk / aflat; sumLean += cdk / lean; sumDirect += cdk / direct; sumPacked += cdk / packed; sumSorted += cdk / sorted; cdkCount++;
            }
            System.out.printf("%-12s %7d %8s %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %9s %9s%n",
                    s.pdbId, s.atomCount, cdkStr, faster, soa, grid, ord, hint, sym, tcd, flat, tc, fl);
        }
        System.out.println();
        int c = Math.max(1, cdkCount);
        return new double[]{ sumFst / c, sumSoa / c, sumGrid / c, sumOrd / c, sumHint / c, sumSym / c, sumSymLA / c, sumVec / c, sumPrn / c, sumDed / c, sumGded / c, sumTcd / c, sumFlat / c, sumAtcd / c, sumAflat / c, sumLean / c, sumDirect / c, sumPacked / c, sumSorted / c };
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
