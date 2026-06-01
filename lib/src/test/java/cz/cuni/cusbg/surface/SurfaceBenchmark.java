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

    @Test
    void benchmarkFasterVsCdk() {
        System.out.printf("%n%-12s %10s %14s %14s %9s%n",
                "structure", "atoms", "faster (ms)", "cdk (ms)", "speedup");
        System.out.println("-----------------------------------------------------------------");

        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            IAtomContainer mol = s.load();

            double faster = medianMillis(() ->
                    new FasterNumericalSurface(mol).getTotalSurfaceArea());

            // CDK's reference NPEs on cobalt (3CI3); skip the comparison column there.
            Double cdk = null;
            if (s != TestStructures.Structure.COBALAMIN) {
                cdk = medianMillis(() -> {
                    NumericalSurface ns = new NumericalSurface(mol, 1.4, 4);
                    ns.calculateSurface();
                    return ns.getTotalSurfaceArea();
                });
            }

            String cdkStr = cdk == null ? "n/a (Co)" : String.format("%.1f", cdk);
            String speedup = cdk == null ? "-" : String.format("%.2fx", cdk / faster);
            System.out.printf("%-12s %10d %14.1f %14s %9s%n",
                    s.pdbId, s.atomCount, faster, cdkStr, speedup);
        }
        System.out.println();
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
