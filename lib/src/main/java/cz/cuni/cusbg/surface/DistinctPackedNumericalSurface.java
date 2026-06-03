package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Surface that fixes the tessellation's point redundancy at the source: it emits ONE point per distinct
 * surviving direction instead of the ~5.7x exact-coincident duplicates the icosahedral tessellation
 * otherwise produces (which {@link FasterNumericalSurface} / {@link PackedNumericalSurface} emit, and
 * which p2rank then sparsifies away).
 *
 * <p>It uses the {@link PackedNumericalSurface} compute pipeline (cell-sorted pruned neighbour build,
 * cached tessellation, cached VdW, copy-free CSR access) but pairs the
 * {@link GlobalDedupWeightedOcclusionScan} (emits one weighted point per surviving distinct direction)
 * with the {@link DistinctFlatSurfacePointStore} (keeps one point per location, counts the multiplicity
 * toward the area). Result: ~5.7x fewer surface points and no need for downstream sparsification, while
 * the <em>areas are bit-for-bit identical</em> to {@code FasterNumericalSurface} (the weight preserves
 * the multiplicity-based count).
 *
 * <p>NOTE: unlike the other variants this is deliberately NOT point-set-identical to CDK
 * {@code NumericalSurface} - it omits the coincident duplicates. The surviving point set equals what
 * sparsification at the exact-coincidence distance produces; areas are exact.
 */
public class DistinctPackedNumericalSurface extends DevSurfaceV1Soa {

    public DistinctPackedNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public DistinctPackedNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel,
                (atoms, ax, ay, az, radius) -> new CoordSortedPrunedSymmetricCellGridNeighborList(atoms, ax, ay, az, radius, solventRadius, VdwRadiusCache::get),
                NeighborOrdering.NONE,
                new GlobalDedupWeightedOcclusionScan(tesslevel),   // scalar weighted dedup; emits distinct directions
                TessellationProvider.CACHED,
                DistinctFlatSurfacePointStore::new,
                false,
                VdwRadiusCache::get,
                true);   // directNeighbors: copy-free CSR access
    }
}
