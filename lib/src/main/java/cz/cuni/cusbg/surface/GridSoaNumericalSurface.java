package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 2: the SoA surface ({@link SoaNumericalSurface}) with the hash-box
 * {@link NeighborList} replaced by the flat {@link CellGridNeighborList}. The neighbor sets are
 * identical (the grid is a drop-in for the same cutoff), so the surface output matches the golden
 * baseline, {@link SoaNumericalSurface}, and {@link FasterNumericalSurface} bit-for-bit; only the
 * neighbor-index cost differs.
 *
 * <p>The grid is selected by passing a {@link NeighborSourceFactory} to the superclass constructor
 * (not by overriding a method), so there is no virtual call during construction.
 */
public class GridSoaNumericalSurface extends SoaNumericalSurface {

    private static final NeighborSourceFactory GRID =
            (atoms, ax, ay, az, radius) -> new CellGridNeighborList(ax, ay, az, radius);

    public GridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public GridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel, GRID, NeighborOrdering.NONE, OcclusionScan.STANDARD);
    }

    /**
     * For subclasses that keep the flat cell grid but supply a neighbor ordering (e.g.
     * {@link OrderedGridSoaNumericalSurface}). The ordering is passed up to the engine constructor
     * rather than installed via an overridable method, so it is fixed before any computation runs.
     */
    protected GridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel,
                                      NeighborOrdering ordering) {
        super(atomContainer, solventRadius, tesslevel, GRID, ordering, OcclusionScan.STANDARD);
    }

    /**
     * For subclasses that keep the flat cell grid but supply a custom occlusion scan (e.g.
     * {@link HintedGridSoaNumericalSurface}), optionally combined with a neighbor ordering. Both
     * strategies are passed up to the engine constructor rather than installed via overridable
     * methods, so they are fixed before any computation runs.
     */
    protected GridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel,
                                      NeighborOrdering ordering, OcclusionScan scan) {
        super(atomContainer, solventRadius, tesslevel, GRID, ordering, scan);
    }
}
