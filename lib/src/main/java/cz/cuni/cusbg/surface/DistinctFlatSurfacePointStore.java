package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Flat point store that keeps ONE point per distinct surviving direction while preserving the
 * multiplicity-weighted count needed for an exact area. Used by {@link DistinctPackedNumericalSurface} with
 * {@link GlobalDedupWeightedOcclusionScan}: the scan calls {@link #addWeighted} once per distinct
 * direction with its tessellation multiplicity, so this store retains a single point per location
 * (no coincident duplicates) but {@link #finishAtom} returns the summed weight, so the engine's
 * {@code area = 4*pi*r^2 * count / pointDensity} stays bit-exact.
 *
 * <p>It therefore produces the de-duplicated surface directly - the redundant ~5.7x coincident points
 * are never emitted - so no downstream sparsification is needed, while areas match the standard
 * (full-multiplicity) surfaces exactly. Point storage is otherwise like {@link FlatSurfacePointStore}
 * (flat {@code double[]} + per-atom CSR offsets; {@code Point3d} materialized lazily).
 */
final class DistinctFlatSurfacePointStore implements SurfacePointStore {

    private final int[] atomStart;   // CSR row pointers over DISTINCT points, length numAtoms+1
    private double[] coords;         // 3 doubles per distinct point
    private int pointCount;          // distinct points stored
    private int cur;
    private long curWeight;          // multiplicity-weighted surviving count for the current atom (for area)

    DistinctFlatSurfacePointStore(int numAtoms) {
        atomStart = new int[numAtoms + 1];
        coords = new double[Math.max(64, numAtoms) * 3];
    }

    DistinctFlatSurfacePointStore(int numAtoms, int estimatedPoints) {
        atomStart = new int[numAtoms + 1];
        coords = new double[Math.max(192, estimatedPoints * 3)];
    }

    @Override public void startAtom(int atomIndex) { cur = atomIndex; atomStart[atomIndex] = pointCount; curWeight = 0; }

    /** Plain add = one distinct point of weight 1 (kept for contract completeness; the scan uses addWeighted). */
    @Override public void add(double x, double y, double z) { store(x, y, z); curWeight += 1; }

    @Override public void addWeighted(double x, double y, double z, int weight) { store(x, y, z); curWeight += weight; }

    private void store(double x, double y, double z) {
        int base = pointCount * 3;
        if (base + 3 > coords.length) coords = Arrays.copyOf(coords, Math.max(coords.length * 2, base + 3));
        coords[base] = x; coords[base + 1] = y; coords[base + 2] = z;
        pointCount++;
    }

    /** Returns the multiplicity-weighted surviving count (NOT the distinct point count) so the area is exact. */
    @Override public int finishAtom() { atomStart[cur + 1] = pointCount; return (int) curWeight; }

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

    @Override public int totalPoints() { return pointCount; }

    @Override public double[] packedXYZ() { return coords; }
}
