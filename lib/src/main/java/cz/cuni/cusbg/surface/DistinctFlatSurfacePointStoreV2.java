package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * V2 distinct-point flat store: keeps ONE point per distinct surviving direction while preserving the
 * multiplicity-weighted count needed for an exact area, like {@link DistinctFlatSurfacePointStore}, but
 * with two improvements that matter for code built on top of it:
 *
 * <ul>
 *   <li><b>Right-sized buffer.</b> The engine's capacity hint estimates FULL-multiplicity survival
 *       (~5.7x more points than a distinct store keeps), so the original store over-allocated its
 *       {@code coords} buffer ~6-8x - ironic for the one surface designed for a small footprint. This
 *       store scales the hint down by the tessellation's expected dedup ratio (see
 *       {@link #DISTINCT_ESTIMATE_DIVISOR}), so it actually allocates a small buffer. It still grows by
 *       doubling if the hint is exceeded, but the divisor keeps headroom so that growth-copy path is
 *       essentially never taken.</li>
 *   <li><b>Unambiguous counts.</b> the distinct point count and the area weight are tracked in two
 *       clearly named fields. {@link #finishAtom()} returns the multiplicity WEIGHT (the engine's area
 *       input); {@link #atomPointCount}/{@link #totalPoints} return the DISTINCT point count. The
 *       original store overloaded a single "count" notion across the two, which was a documented
 *       footgun for anything reading its accessors.</li>
 * </ul>
 *
 * <p><b>Pairing invariant:</b> this store is correct only with a weighted dedup scan
 * ({@link WeightedDedupOcclusionScan} / {@link Vectorized256WeightedDedupOcclusionScan}), which calls
 * {@link #addWeighted} once per distinct surviving direction with that direction's multiplicity.
 * {@link DistinctPackedNumericalSurfaceV2} is the single assembly that wires this pairing, so callers
 * cannot mismatch a scan and store. Point storage is otherwise like {@link FlatSurfacePointStore}
 * (flat {@code double[]} + per-atom CSR offsets; {@code Point3d} materialized lazily).
 */
final class DistinctFlatSurfacePointStoreV2 implements SurfacePointStore {

    /**
     * Divisor applied to the engine's full-multiplicity capacity hint to size the distinct buffer. The
     * icosahedral dedup ratio is ~5.7x and the engine's hint already carries headroom, so the true
     * distinct need is roughly 1/8 of the hint; dividing by 6 (less than the full ratio) leaves ~40%
     * headroom so the doubling-grow path is effectively never taken, while still cutting the bulk of
     * the original ~6-8x over-allocation.
     */
    private static final int DISTINCT_ESTIMATE_DIVISOR = 6;

    private final int[] atomStart;     // CSR row pointers over DISTINCT points, length numAtoms+1
    private double[] coords;           // 3 doubles per distinct point
    private int pointCount;            // DISTINCT points stored so far (across all atoms)
    private int curAtom;
    private long curAtomWeight;        // multiplicity-weighted surviving count for the current atom (area input)

    DistinctFlatSurfacePointStoreV2(int numAtoms, int estimatedPoints) {
        atomStart = new int[numAtoms + 1];
        int distinctEstimate = estimatedPoints / DISTINCT_ESTIMATE_DIVISOR;
        coords = new double[Math.max(192, distinctEstimate * 3)];
    }

    @Override public void startAtom(int atomIndex) { curAtom = atomIndex; atomStart[atomIndex] = pointCount; curAtomWeight = 0; }

    /**
     * Plain add = one distinct point of weight 1. Present for contract completeness; the weighted dedup
     * scans use {@link #addWeighted}.
     */
    @Override public void add(double x, double y, double z) { store(x, y, z); curAtomWeight += 1; }

    @Override public void addWeighted(double x, double y, double z, int weight) { store(x, y, z); curAtomWeight += weight; }

    private void store(double x, double y, double z) {
        int base = pointCount * 3;
        if (base + 3 > coords.length) coords = Arrays.copyOf(coords, Math.max(coords.length * 2, base + 3));
        coords[base] = x; coords[base + 1] = y; coords[base + 2] = z;
        pointCount++;
    }

    /**
     * Finalize the current atom. Returns the multiplicity-WEIGHTED surviving count (NOT the distinct
     * point count): the engine computes {@code area = 4*pi*r^2 * count / pointDensity}, and the weight
     * is what keeps that area bit-for-bit identical to the full-multiplicity surfaces. Distinct point
     * counts are served by {@link #atomPointCount}/{@link #totalPoints}.
     */
    @Override public int finishAtom() { atomStart[curAtom + 1] = pointCount; return (int) curAtomWeight; }

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
