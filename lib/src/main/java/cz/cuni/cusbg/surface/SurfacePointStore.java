package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;
import java.util.List;

/**
 * Accumulates the surviving surface points as the engine scans atoms in order, and serves the
 * point-valued accessors of {@link MolecularSurface}. It is the {@link SurfacePointSink} the scan
 * writes to; the engine calls {@link #startAtom}/{@link #finishAtom} around each atom's scan.
 *
 * <p>Two implementations: {@link ListSurfacePointStore} (the default, one {@code Point3d} object per
 * point, matching the original behavior) and {@link FlatSurfacePointStore} (coordinates in a flat
 * {@code double[]}, {@code Point3d} materialized lazily only when an accessor is called).
 */
interface SurfacePointStore extends SurfacePointSink {

    /** Begin accumulating points for atom {@code atomIndex} (atoms are processed in increasing order). */
    void startAtom(int atomIndex);

    /** Finalize the current atom; returns how many points it contributed. */
    int finishAtom();

    /** All surface points, atom-major in scan order. */
    Point3d[] allPoints();

    /** Atom {@code atomIndex}'s surface points as an unmodifiable list. */
    List<Point3d> atomPoints(int atomIndex);

    /** Number of surface points for atom {@code atomIndex} (cheap; no materialization). */
    int atomPointCount(int atomIndex);

    /** Total number of surface points across all atoms. */
    int totalPoints();

    /**
     * Surface point coordinates packed {@code xyzxyz...} in atom-major scan order, for the zero-copy
     * {@link PackedSurfaceAccess} path. The returned array may be longer than {@code 3*totalPoints()}
     * (a backing buffer exposed by reference); only {@code [0, 3*totalPoints())} is valid. A flat store
     * returns its internal array with no copy; a {@code Point3d}-backed store assembles one.
     */
    double[] packedXYZ();
}
