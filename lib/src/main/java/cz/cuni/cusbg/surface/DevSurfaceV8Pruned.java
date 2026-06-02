package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 7: {@link DevSurfaceV7Simd} plus the per-pair
 * occlusion cutoff, and a 256-bit SIMD scan.
 *
 * <p>Two changes over the (unmodified) vectorized variant, both bit-for-bit output-preserving:
 * <ol>
 *   <li><b>Pruned neighbors.</b> The neighbor index is {@link PrunedSymmetricCellGridNeighborList},
 *       which keeps only neighbors whose expanded spheres overlap ({@code d < R_i + R_j}). Neighbors
 *       farther than that can never bury any tessellation point, so dropping them does not change the
 *       result; it shrinks the per-atom neighbor list the occlusion scan walks (the dominant cost).</li>
 *   <li><b>256-bit scan.</b> Once the neighbor lists are short, 512-bit SIMD bursts (the platform
 *       preferred width on AVX-512) are too brief to amortize AVX-512's frequency-license downclock and
 *       actually run slower; this variant pins the scan to {@link Vectorized256OcclusionScan}.</li>
 * </ol>
 *
 * <p>Same scalar fallback as the vectorized variant: if the {@code jdk.incubator.vector} module is
 * absent the SIMD scan is unavailable and the scalar {@link OcclusionScan#LAST_OCCLUDER_FIRST} is used
 * (the pruning still applies). Output is bit-for-bit identical to {@link FasterNumericalSurface} on
 * either path. {@link #isVectorized()} reports the active scan path.
 */
public class DevSurfaceV8Pruned extends DevSurfaceV1Soa {

    private static final OcclusionScan SCAN = chooseScan();

    private static OcclusionScan chooseScan() {
        try {
            return new Vectorized256OcclusionScan();
        } catch (Throwable notAvailable) {
            return OcclusionScan.LAST_OCCLUDER_FIRST;
        }
    }

    /** True if the SIMD scan is in use; false if this JVM lacks the Vector API and the scalar scan is used. */
    public static boolean isVectorized() {
        return SCAN != OcclusionScan.LAST_OCCLUDER_FIRST;
    }

    public DevSurfaceV8Pruned(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV8Pruned(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        // The neighbor factory captures the solvent radius so the pruned grid can compute per-atom
        // expanded radii; it captures no instance state, so it is safe to pass to super().
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new PrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius),
                NeighborOrdering.NONE, SCAN);
    }
}
