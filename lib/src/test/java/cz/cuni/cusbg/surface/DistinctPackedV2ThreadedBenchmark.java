package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.function.Function;

/**
 * Opt-in multi-threaded throughput benchmark: {@link DistinctPackedNumericalSurfaceV2} vs
 * {@link DistinctPackedNumericalSurface}, building many surfaces concurrently on a fixed thread pool
 * (default 16 threads; override with {@code -Dbench.threads=N}).
 *
 * <p>Where the single-thread harness ({@link DistinctPackedV2Benchmark}) isolates compute time, this
 * one measures aggregate throughput under contention - which is where V2's process-wide shared
 * {@link DirectionMapping} and right-sized {@link DistinctFlatSurfacePointStoreV2} (less per-build
 * allocation -> less GC/allocator contention) are expected to matter more than in the single-thread
 * numbers.
 *
 * <p>Each structure is parsed once and the read-only {@link IAtomContainer} is shared across builds
 * (surface construction only reads it). A round submits every structure {@code REPS} times; we time
 * how long the pool takes to drain a round, median of {@code MEASURE} rounds after {@code WARMUP}.
 *
 * <p>Not part of the normal test run (tagged {@code benchmark}). Run with:
 * <pre>./gradlew benchmark --tests '*DistinctPackedV2ThreadedBenchmark*'
 * ./gradlew benchmark --tests '*DistinctPackedV2ThreadedBenchmark*' -Dbench.threads=8</pre>
 */
@Tag("benchmark")
class DistinctPackedV2ThreadedBenchmark {

    private static final int THREADS = Integer.getInteger("bench.threads", 16);
    private static final int REPS = 80;     // copies of each structure per timed round
    private static final int WARMUP = 2;
    private static final int MEASURE = 3;
    private static final int[] TESS_LEVELS = {4, 3, 2};

    private static final DoubleAdder blackhole = new DoubleAdder();

    @Test
    void benchmarkV2VsV1Threaded() throws Exception {
        System.out.printf("threads=%d  cpus=%d  V1 vectorized=%s  V2 vectorized=%s%n%n",
                THREADS, Runtime.getRuntime().availableProcessors(),
                PackedNumericalSurface.isVectorized(), DistinctPackedNumericalSurfaceV2.isVectorized());

        // Parse once; share the read-only containers across all concurrent builds.
        List<IAtomContainer> mols = new ArrayList<>();
        for (TestStructures.Structure s : TestStructures.Structure.values()) mols.add(s.load());

        ExecutorService pool = Executors.newFixedThreadPool(THREADS);
        try {
            for (int tess : TESS_LEVELS) {
                final int t = tess;
                double v1Ms = timeRound(pool, mols, mol -> new DistinctPackedNumericalSurface(mol, 1.4, t));
                double v2Ms = timeRound(pool, mols, mol -> new DistinctPackedNumericalSurfaceV2(mol, 1.4, t));
                int builds = mols.size() * REPS;
                System.out.printf("=== tess=%d : %d builds/round on %d threads ===%n", tess, builds, THREADS);
                System.out.printf("  V1: %8.1f ms/round  (%8.0f surfaces/s)%n", v1Ms, builds / (v1Ms / 1000.0));
                System.out.printf("  V2: %8.1f ms/round  (%8.0f surfaces/s)%n", v2Ms, builds / (v2Ms / 1000.0));
                System.out.printf("  speedup (V1/V2): %.2fx%n%n", v1Ms / v2Ms);
            }
        } finally {
            pool.shutdown();
        }
    }

    /** Median wall-clock (ms) to build every structure REPS times across the pool, after warmup. */
    private static double timeRound(ExecutorService pool, List<IAtomContainer> mols,
                                    Function<IAtomContainer, MolecularSurface> build) throws Exception {
        List<Callable<Void>> tasks = new ArrayList<>(mols.size() * REPS);
        for (int r = 0; r < REPS; r++) {
            for (IAtomContainer mol : mols) {
                tasks.add(() -> { blackhole.add(build.apply(mol).getTotalSurfaceArea()); return null; });
            }
        }
        for (int w = 0; w < WARMUP; w++) drain(pool, tasks);
        double[] times = new double[MEASURE];
        for (int i = 0; i < MEASURE; i++) {
            long t0 = System.nanoTime();
            drain(pool, tasks);
            times[i] = (System.nanoTime() - t0) / 1_000_000.0;
        }
        Arrays.sort(times);
        return times[times.length / 2];
    }

    private static void drain(ExecutorService pool, List<Callable<Void>> tasks) throws Exception {
        for (Future<Void> f : pool.invokeAll(tasks)) f.get();   // propagate any task failure
    }
}
