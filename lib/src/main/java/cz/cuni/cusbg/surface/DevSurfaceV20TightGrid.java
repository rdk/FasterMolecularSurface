package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * {@link DistinctPackedNumericalSurfaceV2} (the recommended production surface) with its neighbor build
 * swapped for {@link TightCellGridNeighborList} (backlog idea A2): a finer cell grid (cells half the
 * cutoff, ±2 stencil) that visits a smaller candidate volume in the distance pass. Everything else — the
 * weighted dedup SIMD scan, cached tessellation/VdW, distinct right-sized flat store, copy-free CSR
 * access — is V2's, so the only change is inside the neighbor build, which profiling puts at ~47% of CPU
 * at 16 threads / tess 2 on V2 (p2rank's regime). Based on V2 rather than a full-multiplicity dev rung so
 * the comparison is apples-to-apples with the surface that actually ships and so the build is a larger
 * (more visible) share of runtime.
 *
 * <p>Output is <b>bit-for-bit identical</b> to {@link DistinctPackedNumericalSurfaceV2}: the tight grid
 * yields the exact same neighbor set, and the rest of the pipeline is unchanged. Scalar fallback when
 * {@code jdk.incubator.vector} is absent; {@link #isVectorized()} reports the active path.
 */
public class DevSurfaceV20TightGrid extends DevSurfaceV1Soa {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new Vectorized256WeightedDedupOcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the vectorized weighted dedup scan is in use; false on a JVM without the Vector API (scalar fallback). */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public DevSurfaceV20TightGrid(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV20TightGrid(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new TightCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new Vectorized256WeightedDedupOcclusionScan(tesslevel)
                                 : new WeightedDedupOcclusionScan(tesslevel),
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStoreV2::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access
    }
}
