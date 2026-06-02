package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 5: {@link DevSurfaceV4Hinted} with the per-atom cell-grid query
 * replaced by the one-shot symmetric precompute of {@link SymmetricCellGridNeighborList}.
 *
 * <p>Profiling the hinted variant showed that, once the occlusion scan and the neighbor sort were
 * dealt with, the per-atom neighbor query ({@link CellGridNeighborList}) had become the second-largest
 * cost (~26% at tess 2), because it evaluates every pair's distance twice (once from each endpoint).
 * This variant keeps everything else identical (SoA engine, flat cell grid, last-occluder-first scan,
 * no neighbor sort) and only swaps the neighbor index for the symmetric one, which computes each
 * pair's distance once and serves queries as a CSR copy.
 *
 * <p>The neighbor <em>set</em> per atom is unchanged, so the output is bit-for-bit identical to
 * {@link DevSurfaceV4Hinted} / {@link FasterNumericalSurface}.
 */
public class DevSurfaceV5Symmetric extends DevSurfaceV1Soa {

    private static final NeighborSourceFactory SYMMETRIC_GRID =
            (atoms, ax, ay, az, radius) -> new SymmetricCellGridNeighborList(ax, ay, az, radius);

    public DevSurfaceV5Symmetric(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DevSurfaceV5Symmetric(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel, SYMMETRIC_GRID, NeighborOrdering.NONE, OcclusionScan.LAST_OCCLUDER_FIRST);
    }
}
