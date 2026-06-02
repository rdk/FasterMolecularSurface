package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 10: {@link GlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface} that also
 * caches the tessellation arrays process-wide.
 *
 * <p>Even with the dedup mapping cached, the engine still rebuilt the tessellation itself
 * ({@code Tessellate.doTessellate()} plus 240-3840 {@code Point3d} and three {@code double[]}) on
 * every surface. That is a pure function of the tessellation level, so this variant passes
 * {@link TessellationProvider#CACHED} to reuse one shared immutable tessellation across all builds,
 * removing that per-build allocation. It stacks on the full optimized engine: pruned neighbor source,
 * last-occluder hint, and the process-cached 256-bit dedup scan.
 *
 * <p>Output is bit-for-bit identical to {@link FasterNumericalSurface} (the cached tessellation holds
 * the same values the engine would build). Scalar fallback when {@code jdk.incubator.vector} is absent.
 * {@link #isVectorized()} reports the active path.
 */
public class TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface extends SoaNumericalSurface {

    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            new GlobalDedupVectorized256OcclusionScan(4);
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the dedup SIMD scan is in use; false if this JVM lacks the Vector API (scalar fallback). */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED);
    }
}
