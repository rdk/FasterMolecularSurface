package cz.cuni.cusbg.surface;

import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.infra.Blackhole;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.concurrent.TimeUnit;

/**
 * The statistically-rigorous surface benchmark: one parameterized JMH benchmark replacing the four
 * hand-rolled timing classes. Forks + warmup + measurement iterations give run-to-run variance/CI the
 * legacy harness lacked. Run: {@code ./gradlew jmh}.
 *
 * <ul>
 *   <li>{@code variantId} — which surface (the shared {@link SurfaceCatalog}); defaults to the
 *       production + champion subset, full V1..V19 ladder via {@code -p variantId=V1,V9,V18,...}.</li>
 *   <li>{@code tess} — tessellation level (2 = p2rank's operating point, 4 = library default).</li>
 *   <li>{@code consume} — {@code AREA} (just {@code getTotalSurfaceArea()}) vs {@code POINTS} (drains the
 *       points the way p2rank does: zero-copy {@code surfacePointsXYZ()} when available, else Point3d[]).</li>
 * </ul>
 *
 * <p>Thread scaling is JMH's {@code -t} flag, not a separate class. The corpus is the CDK-safe set
 * (excludes 3CI3) so the CDK baseline and every variant are timed on the identical structures, keeping
 * the speedup ratios clean.
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MILLISECONDS)
public class SurfaceBench {

    @Param({"CDK", "FASTER", "PACKED", "DISTINCT_PACKED_V2", "FLOAT", "V18", "V19"})
    public String variantId;

    @Param({"2", "4"})
    public int tess;

    @Param({"AREA", "POINTS"})
    public String consume;

    private static final double SOLVENT = 1.4;

    private IAtomContainer[] corpus;
    private SurfaceCatalog.Factory factory;
    private boolean area;

    @Setup(Level.Trial)
    public void setup() {
        factory = SurfaceCatalog.byId(variantId).factory();
        area = "AREA".equals(consume);
        corpus = TestStructures.loadAllCdkSafe();
    }

    @Benchmark
    public void buildCorpus(Blackhole bh) {
        for (IAtomContainer mol : corpus) {
            MolecularSurface s = factory.build(mol, SOLVENT, tess);
            if (area) {
                bh.consume(s.getTotalSurfaceArea());   // double overload, no boxing
            } else {
                bh.consume(materialize(s));
            }
        }
    }

    /** Drains the surface points the way p2rank does: zero-copy doubles when available, else Point3d[]. */
    private static Object materialize(MolecularSurface s) {
        if (s instanceof PackedSurfaceAccess p) {
            return p.surfacePointsXYZ();
        }
        return s.getAllSurfacePoints();
    }
}
