package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 8: {@link PrunedVectorizedSymmetricHintedGridSoaNumericalSurface} plus tessellation
 * direction deduplication. The full optimized stack: SoA engine, symmetric neighbor precompute, the
 * per-pair occlusion cutoff ({@link PrunedSymmetricCellGridNeighborList}), the last-occluder hint, and
 * a 256-bit SIMD scan that evaluates each <em>distinct</em> tessellation direction once
 * ({@link DedupVectorized256OcclusionScan}) instead of every (≈5.7x duplicated) point.
 *
 * <p>All three layers are bit-for-bit output-preserving, so the result is identical to
 * {@link FasterNumericalSurface}: the cutoff only drops neighbors that can never bury; the dedup scan
 * re-expands surviving directions into the original point order and multiplicity.
 *
 * <p>The dedup scan memoizes the per-build direction mapping, so unlike the other variants' stateless
 * shared scans it must be created fresh per surface (done here). When the {@code jdk.incubator.vector}
 * module is absent the SIMD scan is unavailable and this variant falls back to the scalar
 * {@link OcclusionScan#LAST_OCCLUDER_FIRST} (the pruning still applies, but the dedup does not).
 * {@link #isVectorized()} reports whether the dedup SIMD path is active.
 */
public class DedupVectorizedSymmetricHintedGridSoaNumericalSurface extends SoaNumericalSurface {

    /** Whether the dedup SIMD scan is usable on this JVM (probed once; the scan itself is per-build). */
    private static final boolean VECTOR_AVAILABLE = probeVector();

    private static boolean probeVector() {
        try {
            // touching the class runs its static init, which resolves the incubator Vector API
            new DedupVectorized256OcclusionScan();
            return true;
        } catch (Throwable notAvailable) {
            return false;
        }
    }

    /** True if the dedup SIMD scan is in use; false if this JVM lacks the Vector API (scalar fallback). */
    public static boolean isVectorized() {
        return VECTOR_AVAILABLE;
    }

    public DedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DedupVectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        // pruned neighbor factory captures solvent (no instance state); the dedup scan is created fresh
        // per surface because it memoizes per-build state. The new-scan branch is only taken when the
        // Vector API is available, so it never loads the incubator class on a JVM without it.
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new DedupVectorized256OcclusionScan() : OcclusionScan.LAST_OCCLUDER_FIRST);
    }
}
