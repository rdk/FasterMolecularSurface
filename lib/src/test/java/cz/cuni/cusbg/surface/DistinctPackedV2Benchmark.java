package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import javax.vecmath.Point3d;
import java.util.Arrays;
import java.util.function.Function;
import java.util.function.Supplier;

/**
 * Opt-in head-to-head benchmark: {@link DistinctPackedNumericalSurfaceV2} (vectorized weighted dedup
 * scan + right-sized distinct store) vs {@link DistinctPackedNumericalSurface} (scalar weighted dedup
 * scan + full-multiplicity-sized store). Both emit the same distinct surface with bit-exact areas, so
 * point counts should match (ratio ~1.0) and areas must be equal; the figure of merit is wall-clock.
 *
 * <p>Not part of the normal test run (tagged {@code benchmark}). Run with:
 * <pre>./gradlew benchmark --tests '*DistinctPackedV2Benchmark*'</pre>
 * Coarse median-of-N timing after warm-up (not JMH); for relative comparison on one machine.
 */
@Tag("benchmark")
class DistinctPackedV2Benchmark {

    private static final int WARMUP = 3;
    private static final int MEASURE = 5;
    private static final int[] TESS_LEVELS = {4, 3, 2};   // 4 = library default, 2 = p2rank's operating point

    /** Keeps the JIT from eliminating the timed work. */
    private static volatile double blackhole;

    @Test
    void benchmarkV2VsV1() {
        System.out.printf("V1 vectorized: %s   V2 vectorized: %s%n%n",
                PackedNumericalSurface.isVectorized(), DistinctPackedNumericalSurfaceV2.isVectorized());

        for (int tess : TESS_LEVELS) {
            runTable(tess,
                    mol -> new DistinctPackedNumericalSurface(mol, 1.4, tess),
                    mol -> new DistinctPackedNumericalSurfaceV2(mol, 1.4, tess));
        }
    }

    private void runTable(int tess,
                          Function<IAtomContainer, MolecularSurface> v1,
                          Function<IAtomContainer, MolecularSurface> v2) {
        System.out.printf("=== DistinctPacked V1 vs V2 @ tess=%d ===%n", tess);
        System.out.printf("%-12s %7s | %9s %10s %10s | %9s %10s %10s | %8s %8s %10s%n",
                "structure", "atoms",
                "v1(ms)", "v1 pts", "v1 area",
                "v2(ms)", "v2 pts", "v2 area",
                "speedup", "pts x", "area==");
        System.out.println("-------------------------------------------------------------------------------------------------------------------------------");

        double sumSpeedup = 0;
        int n = 0;
        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            IAtomContainer mol = s.load();

            double v1Ms = medianMillis(() -> v1.apply(mol).getTotalSurfaceArea());
            double v2Ms = medianMillis(() -> v2.apply(mol).getTotalSurfaceArea());

            MolecularSurface a = v1.apply(mol);
            MolecularSurface b = v2.apply(mol);
            int v1Pts = countPoints(a);
            int v2Pts = countPoints(b);
            double v1Area = a.getTotalSurfaceArea();
            double v2Area = b.getTotalSurfaceArea();
            boolean areaEq = v1Area == v2Area;   // both distinct surfaces are bit-exact in area

            double speedup = v1Ms / v2Ms;            // >1 = V2 faster
            double ptsRatio = (double) v1Pts / v2Pts; // expect ~1.0 (same point set)
            sumSpeedup += speedup; n++;

            System.out.printf("%-12s %7d | %9.1f %10d %10.1f | %9.1f %10d %10.1f | %7.2fx %7.2fx %10s%n",
                    s.pdbId, s.atomCount,
                    v1Ms, v1Pts, v1Area,
                    v2Ms, v2Pts, v2Area,
                    speedup, ptsRatio, areaEq ? "yes" : "NO!");
        }
        System.out.printf("%-12s %7s | %9s %10s %10s | %9s %10s %10s | %7.2fx%n%n",
                "MEAN", "", "", "", "", "", "", "", sumSpeedup / n);
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
