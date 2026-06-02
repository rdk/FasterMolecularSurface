package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Low-allocation {@link SurfacePointStore}: surviving point coordinates are written into a single flat
 * {@code double[]} ({@code xyzxyz...}) with per-atom CSR offsets, instead of one {@code Point3d} object
 * per point plus a per-atom list. {@code Point3d} objects are materialized lazily, only when an object
 * accessor is called; callers that need only areas (or the flat coordinates) never allocate them.
 *
 * <p>This is optimization B: it removes the per-point {@code Point3d} churn (a top allocation source),
 * which mainly helps GC pressure / memory bandwidth at scale. The point values, ordering (atom-major,
 * scan order within an atom), and counts are identical to {@link ListSurfacePointStore}, so the surface
 * output is unchanged.
 */
final class FlatSurfacePointStore implements SurfacePointStore {

    private final int[] atomStart;   // CSR row pointers, length numAtoms+1 (in points, not doubles)
    private double[] coords;         // 3 doubles per point
    private int pointCount;
    private int cur;

    FlatSurfacePointStore(int numAtoms) {
        atomStart = new int[numAtoms + 1];
        coords = new double[Math.max(64, numAtoms) * 3];   // grows on demand
    }

    @Override public void startAtom(int atomIndex) { cur = atomIndex; atomStart[atomIndex] = pointCount; }

    @Override
    public void add(double x, double y, double z) {
        int base = pointCount * 3;
        if (base + 3 > coords.length) coords = Arrays.copyOf(coords, Math.max(coords.length * 2, base + 3));
        coords[base] = x; coords[base + 1] = y; coords[base + 2] = z;
        pointCount++;
    }

    @Override public int finishAtom() { atomStart[cur + 1] = pointCount; return pointCount - atomStart[cur]; }

    @Override
    public Point3d[] allPoints() {
        Point3d[] ret = new Point3d[pointCount];
        for (int p = 0; p < pointCount; p++) ret[p] = new Point3d(coords[3 * p], coords[3 * p + 1], coords[3 * p + 2]);
        return ret;
    }

    @Override
    public List<Point3d> atomPoints(int atomIndex) {
        int s = atomStart[atomIndex], e = atomStart[atomIndex + 1];
        List<Point3d> list = new ArrayList<>(e - s);
        for (int p = s; p < e; p++) list.add(new Point3d(coords[3 * p], coords[3 * p + 1], coords[3 * p + 2]));
        return Collections.unmodifiableList(list);
    }

    @Override public int atomPointCount(int atomIndex) { return atomStart[atomIndex + 1] - atomStart[atomIndex]; }
}
