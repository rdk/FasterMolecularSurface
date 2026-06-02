package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 6: {@link SymmetricHintedGridSoaNumericalSurface} with the occlusion scan - the
 * ~70% hot spot in that variant's profile - vectorized via the {@code jdk.incubator.vector} API. SoA
 * layout, symmetric neighbor precompute, and the last-occluder-first hint are all unchanged; only the
 * inner per-neighbor buried test runs a lane-width at a time (see {@link VectorizedOcclusionScan}).
 *
 * <p><b>Fallback for older / unconfigured JVMs.</b> The Vector API is an incubator module
 * ({@code jdk.incubator.vector}) that must be present and added to the module graph
 * ({@code --add-modules jdk.incubator.vector}). When it is not - an older JDK, or a run without that
 * flag - selecting the vectorized scan throws (a {@code NoClassDefFoundError} /
 * {@code ExceptionInInitializerError} from {@link VectorizedOcclusionScan}'s initializer), which is
 * caught once at class-load time and this variant transparently uses the scalar
 * {@link OcclusionScan#LAST_OCCLUDER_FIRST} instead. In that case it is byte-for-byte the same
 * computation as {@link SymmetricHintedGridSoaNumericalSurface}. {@link #isVectorized()} reports which
 * path is active.
 *
 * <p>Output is bit-for-bit identical to {@link SymmetricHintedGridSoaNumericalSurface} /
 * {@link FasterNumericalSurface} on either path (the vectorized scan reproduces the scalar dot product
 * lane-for-lane, no FMA).
 */
public class VectorizedSymmetricHintedGridSoaNumericalSurface extends SoaNumericalSurface {

    private static final NeighborSourceFactory SYMMETRIC_GRID =
            (atoms, ax, ay, az, radius) -> new SymmetricCellGridNeighborList(ax, ay, az, radius);

    /** Vectorized scan when the Vector API is available on this JVM, else the scalar hint scan. */
    private static final OcclusionScan SCAN = chooseScan();

    private static OcclusionScan chooseScan() {
        try {
            // touching VectorizedOcclusionScan runs its static initializer, which resolves the
            // incubator Vector API; if the module is absent that throws here and we fall back.
            return new VectorizedOcclusionScan();
        } catch (Throwable notAvailable) {
            return OcclusionScan.LAST_OCCLUDER_FIRST;
        }
    }

    /** True if the SIMD scan is in use; false if this JVM lacks the Vector API and the scalar scan is used. */
    public static boolean isVectorized() {
        // compare against the scalar fallback (not the vector class literal, which could itself fail to load)
        return SCAN != OcclusionScan.LAST_OCCLUDER_FIRST;
    }

    public VectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public VectorizedSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel, SYMMETRIC_GRID, NeighborOrdering.NONE, SCAN);
    }
}
