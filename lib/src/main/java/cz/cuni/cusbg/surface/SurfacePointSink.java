package cz.cuni.cusbg.surface;

/**
 * Primitive sink the occlusion scan writes surviving surface points to, one {@code (x,y,z)} at a time.
 *
 * <p>Replaces passing a {@code List<Point3d>} to the scan: by handing the scan a primitive sink, an
 * implementation can store coordinates in a flat {@code double[]} (see {@link FlatSurfacePointStore})
 * instead of allocating a {@code Point3d} per point. The default store ({@link ListSurfacePointStore})
 * still builds {@code Point3d} objects, so existing variants are unchanged.
 *
 * <p>Coordinates are the final surface-point positions ({@code totalRadius * direction + atomCenter}),
 * already transformed by the scan.
 */
interface SurfacePointSink {
    void add(double x, double y, double z);
}
