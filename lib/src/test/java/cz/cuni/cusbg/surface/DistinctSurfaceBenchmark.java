package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import javax.vecmath.Point3d;
import java.util.Arrays;
import java.util.function.Supplier;

/**
 * Opt-in benchmark for the distinct-point surfaces against the originals they are derived from:
 * <ul>
 *   <li>{@link DistinctFasterNumericalSurface} vs {@link FasterNumericalSurface}</li>
 *   <li>{@link DistinctPackedNumericalSurface} vs {@link PackedNumericalSurface}</li>
 * </ul>
 *
 * <p>The distinct variants emit ONE point per surviving distinct direction (no ~5.7x coincident
 * duplicates) while keeping the area bit-exact. So we report not just wall-clock time but also the
 * emitted point count, since the whole point of the distinct surfaces is to skip the downstream
 * sparsification the originals require: a fair end-to-end comparison must weigh "compute + sparsify"
 * for the originals against "compute only" for the distinct variants. This harness measures the raw
 * surface-compute time here; see the p2rank {@code analyze surface-strategies} command for the
 * end-to-end (compute + sparsify) picture.
 *
 * <p>Not part of the normal test run (tagged {@code benchmark}). Run with:
 * <pre>./gradlew benchmark --tests '*DistinctSurfaceBenchmark*'</pre>
 * Coarse median-of-N timing after warm-up (not JMH); for relative comparison on one machine.
 */
@Tag("benchmark")
class DistinctSurfaceBenchmark {

    private static final int WARMUP = 3;
    private static final int MEASURE = 5;
    private static final int[] TESS_LEVELS = {4, 3, 2};   // 4 = library default, 2 = p2rank's operating point

    /** Keeps the JIT from eliminating the timed work. */
    private static volatile double blackhole;

    @Test
    void benchmarkDistinctVsOriginal() {
        System.out.printf("vectorized scan active (jdk.incubator.vector available): %s%n%n",
                DevSurfaceV7Simd.isVectorized());

        for (int tess : TESS_LEVELS) {
            runTable("Faster", tess,
                    mol -> new FasterNumericalSurface(mol, 1.4, tess),
                    mol -> new DistinctFasterNumericalSurface(mol, 1.4, tess));
            runTable("Packed", tess,
                    mol -> new PackedNumericalSurface(mol, 1.4, tess),
                    mol -> new DistinctPackedNumericalSurface(mol, 1.4, tess));
        }
    }

    private void runTable(String family, int tess,
                          java.util.function.Function<IAtomContainer, MolecularSurface> original,
                          java.util.function.Function<IAtomContainer, MolecularSurface> distinct) {
        System.out.printf("=== %s: original vs Distinct%s @ tess=%d ===%n", family, family, tess);
        System.out.printf("%-12s %7s | %9s %10s %10s | %9s %10s %10s | %8s %8s %10s%n",
                "structure", "atoms",
                "orig(ms)", "orig pts", "orig area",
                "dist(ms)", "dist pts", "dist area",
                "time x", "pts x", "area==");
        System.out.println("-------------------------------------------------------------------------------------------------------------------------------");

        double sumTimeRatio = 0, sumPtsRatio = 0;
        int n = 0;
        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            IAtomContainer mol = s.load();

            double origMs = medianMillis(() -> original.apply(mol).getTotalSurfaceArea());
            double distMs = medianMillis(() -> distinct.apply(mol).getTotalSurfaceArea());

            MolecularSurface o = original.apply(mol);
            MolecularSurface d = distinct.apply(mol);
            int origPts = countPoints(o);
            int distPts = countPoints(d);
            double origArea = o.getTotalSurfaceArea();
            double distArea = d.getTotalSurfaceArea();
            boolean areaEq = origArea == distArea;   // distinct surfaces are bit-exact in area

            double timeRatio = origMs / distMs;          // >1 = distinct faster
            double ptsRatio = (double) origPts / distPts; // ~5.7 = dedup factor
            sumTimeRatio += timeRatio; sumPtsRatio += ptsRatio; n++;

            System.out.printf("%-12s %7d | %9.1f %10d %10.1f | %9.1f %10d %10.1f | %7.2fx %7.2fx %10s%n",
                    s.pdbId, s.atomCount,
                    origMs, origPts, origArea,
                    distMs, distPts, distArea,
                    timeRatio, ptsRatio, areaEq ? "yes" : "NO!");
        }
        System.out.printf("%-12s %7s | %9s %10s %10s | %9s %10s %10s | %7.2fx %7.2fx%n%n",
                "MEAN", "", "", "", "", "", "", "",
                sumTimeRatio / n, sumPtsRatio / n);
    }

    private static int countPoints(MolecularSurface s) {
        Point3d[] pts = s.getAllSurfacePoints();
        return pts.length;
    }

    /** Median wall-clock time in milliseconds over MEASURE runs after WARMUP warm-up runs. */
    private static double medianMillis(Supplier<Double> work) {
        for (int i = 0; i < WARMUP; i++) blackhole += work.get();
        double[] times = new double[MEASURE];
        for (int i = 0; i < MEASURE; i++) {
            long t0 = System.nanoTime();
            blackhole += work.get();
            long t1 = System.nanoTime();
            times[i] = (t1 - t0) / 1_000_000.0;
        }
        Arrays.sort(times);
        return times[times.length / 2];
    }
}
