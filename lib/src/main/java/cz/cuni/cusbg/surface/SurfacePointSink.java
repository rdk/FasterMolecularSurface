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

    /**
     * Add a surface point that stands for {@code weight} coincident tessellation points (the icosahedral
     * tessellation emits each direction with multiplicity >=5, so a surviving direction maps to that many
     * exact-coincident points). The default emits {@code weight} copies via {@link #add}, so a normal
     * store reproduces the full multiplicity (byte-identical to emitting each point). A deduplicating
     * store (see {@link DistinctFlatSurfacePointStore}) overrides this to keep one point while counting
     * {@code weight} toward the area, yielding a distinct-point surface with a bit-exact area.
     */
    default void addWeighted(double x, double y, double z, int weight) {
        for (int i = 0; i < weight; i++) add(x, y, z);
    }
}
