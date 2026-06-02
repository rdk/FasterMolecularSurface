package cz.cuni.cusbg.surface;

/**
 * Strategy for reordering an atom's per-neighbor scratch (the parallel
 * {@code diffX/diffY/diffZ/thresh} arrays, first {@code numNeighbors} entries) before the occlusion
 * loop in {@link SoaNumericalSurface}.
 *
 * <p>Reordering is purely a performance lever and must be output-preserving: a tessellation point is
 * buried iff <em>any</em> neighbor buries it, so the order in which neighbors are scanned cannot change
 * which points survive, their coordinates, or the computed areas — only how early
 * {@code collectPoints}' {@code break} fires. Implementations must permute all four arrays in lockstep.
 *
 * <p>Supplied to {@link SoaNumericalSurface} through its constructor (not via an overridable method),
 * so the strategy is fixed before any computation runs and there is no virtual dispatch on a
 * partially-constructed subclass. For that reason implementations must be stateless (read only their
 * array arguments), which the provided implementations — {@link #NONE} and the threshold sort in
 * {@link OrderedGridSoaNumericalSurface} — are.
 */
@FunctionalInterface
interface NeighborOrdering {

    /** Reorder entries {@code [0, numNeighbors)} of the four parallel arrays together. */
    void order(double[] diffX, double[] diffY, double[] diffZ, double[] thresh, int numNeighbors);

    /** Identity ordering: leaves the neighbor scratch untouched. */
    NeighborOrdering NONE = (diffX, diffY, diffZ, thresh, numNeighbors) -> { };
}
