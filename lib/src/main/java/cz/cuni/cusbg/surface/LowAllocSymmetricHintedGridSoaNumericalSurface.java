package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Low-allocation counterpart of {@link SymmetricHintedGridSoaNumericalSurface}: identical engine (SoA,
 * symmetric neighbor precompute, last-occluder-first scan, no sort) but with the neighbor index swapped
 * for {@link LowAllocSymmetricCellGridNeighborList}, which builds its CSR adjacency in two half-stencil
 * passes instead of buffering the discovered edges, dropping the {@code int[]} churn that profiling of
 * {@link SymmetricHintedGridSoaNumericalSurface} flagged.
 *
 * <p><b>WARNING - this variant is intentionally slower.</b> It is a GC-vs-latency tradeoff, not a
 * strict improvement. Measured on 4HHB at tess 2 (GraalVM 25) it allocates ~63% less per build but runs
 * ~2-7% slower than {@link SymmetricHintedGridSoaNumericalSurface} (the gap widens at lower
 * tessellation, where the neighbor build is a larger share of total time), because it repeats the
 * distance pass rather than store the edges. Choose it only when allocation / GC pressure / throughput
 * at scale matters more than single-run latency. For lowest latency, use
 * {@link SymmetricHintedGridSoaNumericalSurface}.
 *
 * <p>Output is bit-for-bit identical to {@link SymmetricHintedGridSoaNumericalSurface} /
 * {@link FasterNumericalSurface}.
 */
public class LowAllocSymmetricHintedGridSoaNumericalSurface extends SoaNumericalSurface {

    private static final NeighborSourceFactory LOW_ALLOC_SYMMETRIC_GRID =
            (atoms, ax, ay, az, radius) -> new LowAllocSymmetricCellGridNeighborList(ax, ay, az, radius);

    public LowAllocSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public LowAllocSymmetricHintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel, LOW_ALLOC_SYMMETRIC_GRID, NeighborOrdering.NONE, OcclusionScan.LAST_OCCLUDER_FIRST);
    }
}
