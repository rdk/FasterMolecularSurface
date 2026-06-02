package cz.cuni.cusbg.surface;

import org.openscience.cdk.geometry.surface.Tessellate;

import javax.vecmath.Point3d;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The unit-sphere tessellation as flat structure-of-arrays: the per-point direction coordinates
 * ({@link #tx}/{@link #ty}/{@link #tz}), the point count, and the {@code pointDensity} used to
 * normalize areas. It is a pure function of the tessellation level (CDK's {@code "ico"} tessellation
 * is deterministic: identical points in identical order every time), so it can be built once and
 * shared read-only.
 *
 * <p>Use {@link #build(int)} to construct fresh (the engine default), or {@link #cached(int)} to fetch
 * a process-wide shared instance keyed by level (see {@link TessellationProvider}). The arrays are
 * never mutated after construction, so a cached instance is safe to read from any number of threads.
 */
final class Tessellation {

    final double[] tx, ty, tz;
    final int numTess;
    final int pointDensity;

    private static final ConcurrentHashMap<Integer, Tessellation> CACHE = new ConcurrentHashMap<>();

    private Tessellation(double[] tx, double[] ty, double[] tz, int numTess, int pointDensity) {
        this.tx = tx; this.ty = ty; this.tz = tz;
        this.numTess = numTess;
        this.pointDensity = pointDensity;
    }

    /** Build the tessellation for a level from scratch (allocates the CDK {@code Tessellate} graph). */
    static Tessellation build(int tesslevel) {
        Tessellate tess = new Tessellate("ico", tesslevel);
        tess.doTessellate();
        Point3d[] tessPoints = tess.getTessAsPoint3ds();
        int numTess = tessPoints.length;
        double[] tx = new double[numTess], ty = new double[numTess], tz = new double[numTess];
        for (int t = 0; t < numTess; t++) {
            tx[t] = tessPoints[t].x; ty[t] = tessPoints[t].y; tz[t] = tessPoints[t].z;
        }
        return new Tessellation(tx, ty, tz, numTess, tess.getNumberOfTriangles() * 3);
    }

    /** Return the shared immutable tessellation for a level, building it once on first request. */
    static Tessellation cached(int tesslevel) {
        return CACHE.computeIfAbsent(tesslevel, Tessellation::build);
    }
}
