package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

import javax.vecmath.Point3d;

/**
 * Optimization step 4: {@link GridSoaNumericalSurface} plus a <em>last-occluder-first hint</em> in
 * the occlusion scan, and <em>no</em> neighbor sort.
 *
 * <p>Motivation: profiling {@link OrderedGridSoaNumericalSurface} at low tessellation showed the
 * per-atom {@code thresh} sort had become the dominant cost — it front-loads the largest occluders
 * once per atom, but pays an O(m log m) sort to do it. This variant chases the same early-exit win
 * without sorting: the tessellation points of one atom are visited in a fixed order and are spatially
 * coherent, so the neighbor that buried the previous point very often buries the next one too. The
 * scan therefore remembers the last burying neighbor and tests it first; on a hit it exits after a
 * single comparison, on a miss it falls back to the full in-order scan. The hint is one {@code int}
 * local per atom and costs nothing to maintain.
 *
 * <p>Output is identical to {@link GridSoaNumericalSurface} / {@link FasterNumericalSurface}
 * bit-for-bit: a point is buried iff ANY neighbor buries it, so testing one neighbor out of order
 * first cannot change which points survive, their coordinates, or the areas - only how early the scan
 * stops. (On a hint miss every neighbor, including the hinted one, is still tested.)
 *
 * <p>The hint is a {@link OcclusionScan} whose only state is a local of its {@code collect} method,
 * reset per atom, so it is stateless across atoms and safe to fix before construction (it is supplied
 * through the constructor, not an overridable method).
 */
public class HintedGridSoaNumericalSurface extends GridSoaNumericalSurface {

    /**
     * Stateless occlusion scan with a last-occluder-first hint. Identical predicate to
     * {@link OcclusionScan#STANDARD}; only the per-point test order differs.
     */
    private static final OcclusionScan LAST_OCCLUDER_FIRST =
            (tx, ty, tz, numTess, numNeighbors, diffX, diffY, diffZ, thresh,
             totalRadius, atomX, atomY, atomZ, points) -> {
        int last = -1;   // neighbor that buried the previous point, or -1 (per-atom local; reset here)
        for (int t = 0; t < numTess; t++) {
            double px = tx[t], py = ty[t], pz = tz[t];

            // try the cached last occluder first: consecutive tessellation points are spatially close,
            // so the same neighbor frequently buries both
            if (last >= 0 && diffX[last] * px + diffY[last] * py + diffZ[last] * pz > thresh[last]) {
                continue;   // buried, last unchanged
            }

            boolean buried = false;
            for (int k = 0; k < numNeighbors; k++) {
                if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) {
                    last = k;   // remember this occluder for the next point
                    buried = true;
                    break;
                }
            }
            if (!buried) {
                points.add(new Point3d(totalRadius * px + atomX, totalRadius * py + atomY, totalRadius * pz + atomZ));
            }
        }
    };

    public HintedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public HintedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel, NeighborOrdering.NONE, LAST_OCCLUDER_FIRST);
    }
}
