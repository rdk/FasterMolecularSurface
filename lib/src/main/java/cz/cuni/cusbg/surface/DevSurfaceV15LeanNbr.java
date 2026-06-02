package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step on top of {@link DevSurfaceV11CachedTess}, targeting the neighbor-list build that
 * profiling identified as the new bottleneck once the occlusion scan was vectorized and cached. At
 * tess 2 the neighbor source construction is ~25% of single-thread time and ~44% at 16 threads (the
 * hottest method there), and its {@code int[]} allocation is 51-61% of all allocation - so it, not the
 * SAS math, now bounds throughput.
 *
 * <p>Two changes, both bit-for-bit output-preserving:
 * <ul>
 *   <li><b>#1 Lean two-pass neighbor build.</b> Swaps {@link PrunedSymmetricCellGridNeighborList} for
 *       {@link LowAllocPrunedSymmetricCellGridNeighborList}: instead of buffering kept edges in two
 *       growable {@code IntArrayList(n*8)} and then assembling the CSR, it counts degrees in one stencil
 *       pass and fills the CSR directly in a second, eliminating both edge buffers (the dominant
 *       allocation) and the {@code IntArrayList.add} churn (15% of single-thread CPU). It trades that
 *       for a repeated distance pass - cheap arithmetic against the bandwidth-bound allocation it
 *       removes, which is the right trade in the 16-thread regime.</li>
 *   <li><b>#3 Cached VdW radii, engine side too.</b> The van der Waals radius is a pure function of
 *       element symbol but was resolved per atom via a {@code toLowerCase} + periodic-table lookup
 *       <em>twice</em> per build (once in the engine's {@link #init()} for the expanded radii, once in
 *       the neighbor source for its cutoff). Both now go through {@link VdwRadiusCache}, so each distinct
 *       symbol is resolved once per process. (V12/V14 already cached the neighbor-source call; this also
 *       caches the engine call.)</li>
 * </ul>
 *
 * <p>Keeps the rest of the V11 stack unchanged - process-cached 256-bit dedup scan, last-occluder hint,
 * cached tessellation, and the same {@link ListSurfacePointStore} as V11 (so this isolates the
 * neighbor-build + VdW changes; the flat-store change is V12's axis). Result is identical to
 * {@link FasterNumericalSurface}. Scalar fallback when {@code jdk.incubator.vector} is absent;
 * {@link #isVectorized()} reports the active path.
 *
 * <p><b>Measured outcome (GraalVM 25, median of 5-7 reps; kept as a documented negative result, not a
 * champion):</b> change #1 does <em>not</em> pay off, because the neighbor build is CPU-bound on its
 * distance test, not allocation-bound, so doubling the distance pass costs more than the {@code int[]}
 * it removes - the same tradeoff {@link LowAllocSymmetricCellGridNeighborList} (V6) showed for the
 * unpruned source, here with a heavier per-pair cutoff:
 * <ul>
 *   <li><b>tess 2</b> (p2rank's operating point): ~0.86x V11 single-thread, ~0.88x V11 at 16 threads -
 *       a 12-14% <em>regression</em>, since the neighbor build is the largest single cost at tess 2.</li>
 *   <li><b>tess 4</b>: ~0.94x V11 single-thread; ~1.00x V11 at 16 threads (break-even, where the
 *       occlusion scan dominates and GC pressure is high enough that the allocation saving offsets the
 *       extra compute).</li>
 * </ul>
 * No measured configuration beats V11. The right lever for the neighbor build is instead to remove the
 * redundant per-query CSR-&gt;{@code IntArrayList} copy in {@code getNeighborsInto} (no added compute),
 * not to trade allocation for a repeated distance pass. The engine-side cached-VdW plumbing (#3) added
 * here is sound and reusable by a future variant; its benefit is small and swamped by the #1 regression.
 * See {@code docs/performance-lessons.md}.
 */
public class DevSurfaceV15LeanNbr extends DevSurfaceV1Soa {

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

    public DevSurfaceV15LeanNbr(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV15LeanNbr(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new LowAllocPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                VECTOR_AVAILABLE ? new GlobalDedupVectorized256OcclusionScan(tesslevel) : OcclusionScan.LAST_OCCLUDER_FIRST,
                TessellationProvider.CACHED,
                ListSurfacePointStore::new,
                false,
                VdwRadiusCache::get);
    }
}
