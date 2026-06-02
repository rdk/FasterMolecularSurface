package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Default {@link SurfacePointStore}: one {@code Point3d} object per surviving point, stored as a
 * per-atom {@code List<Point3d>}. Reproduces the engine's original storage exactly (a reused builder
 * list copied into a fresh per-atom list), so variants using it are byte-identical to before.
 */
final class ListSurfacePointStore implements SurfacePointStore {

    private final List<Point3d>[] perAtom;
    private final List<Point3d> reused = new ArrayList<>(512);
    private int cur;

    @SuppressWarnings("unchecked")
    ListSurfacePointStore(int numAtoms) {
        perAtom = (List<Point3d>[]) new List[numAtoms];
    }

    /** The point-count hint is irrelevant to the list store (no monolithic buffer to pre-size); ignored. */
    ListSurfacePointStore(int numAtoms, int estimatedPoints) {
        this(numAtoms);
    }

    @Override public void startAtom(int atomIndex) { cur = atomIndex; reused.clear(); }
    @Override public void add(double x, double y, double z) { reused.add(new Point3d(x, y, z)); }
    @Override public int finishAtom() { perAtom[cur] = new ArrayList<>(reused); return perAtom[cur].size(); }

    @Override
    public Point3d[] allPoints() {
        int npt = 0;
        for (List<Point3d> s : perAtom) npt += s.size();
        Point3d[] ret = new Point3d[npt];
        int j = 0;
        for (List<Point3d> pts : perAtom)
            for (Point3d p : pts) ret[j++] = p;
        return ret;
    }

    @Override public List<Point3d> atomPoints(int atomIndex) { return Collections.unmodifiableList(perAtom[atomIndex]); }
    @Override public int atomPointCount(int atomIndex) { return perAtom[atomIndex].size(); }

    @Override
    public int totalPoints() {
        int npt = 0;
        for (List<Point3d> s : perAtom) npt += s.size();
        return npt;
    }

    /** Fallback (not zero-copy): assembles a fresh {@code double[]} from the per-point {@code Point3d}s. */
    @Override
    public double[] packedXYZ() {
        double[] xyz = new double[3 * totalPoints()];
        int j = 0;
        for (List<Point3d> pts : perAtom)
            for (Point3d p : pts) { xyz[j++] = p.x; xyz[j++] = p.y; xyz[j++] = p.z; }
        return xyz;
    }
}
