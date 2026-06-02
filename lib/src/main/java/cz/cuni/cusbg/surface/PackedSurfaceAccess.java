package cz.cuni.cusbg.surface;

/**
 * Optional capability: zero-copy access to the surface points as a packed primitive array, so a caller
 * that builds its own point representation (e.g. p2rank's {@code Point}, three primitive doubles) never
 * has to go through the {@code Point3d[]} of {@link MolecularSurface#getAllSurfacePoints()} - which, for
 * a store that already holds coordinates in a flat {@code double[]}, would force materializing one
 * {@code Point3d} per point (plus the array) only to read three doubles back out and discard it.
 *
 * <p>A surface backed by {@link FlatSurfacePointStore} serves this by returning its internal coordinate
 * array <em>by reference</em> (no allocation, no copy); a {@link ListSurfacePointStore}-backed surface
 * falls back to assembling the array once. Callers detect support with {@code instanceof}.
 */
public interface PackedSurfaceAccess {

    /**
     * Surface point coordinates packed as {@code x0,y0,z0, x1,y1,z1, ...}, in the same atom-major scan
     * order as {@link MolecularSurface#getAllSurfacePoints()}.
     *
     * <p><b>The returned array may be longer than needed</b> (a reused/grown buffer exposed by
     * reference for zero copy): valid data is exactly {@code [0, 3 * surfacePointCount())}. Treat it as
     * read-only and do not retain it beyond the surface's lifetime.
     */
    double[] surfacePointsXYZ();

    /** Number of surface points; valid coordinate data is {@code surfacePointsXYZ()[0 .. 3*count)}. */
    int surfacePointCount();
}
