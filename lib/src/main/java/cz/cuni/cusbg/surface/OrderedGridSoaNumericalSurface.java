package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Optimization step 3: {@link GridSoaNumericalSurface} plus neighbor ordering.
 *
 * <p>Before the occlusion loop, the per-atom neighbors are sorted by ascending {@code thresh}. The
 * occlusion test buries a tessellation point as soon as some neighbor satisfies
 * {@code diff.p > thresh}; a neighbor with a smaller {@code thresh} covers a larger spherical cap
 * (occludes more points), so visiting those first makes {@code collectPoints}' early-exit
 * {@code break} fire sooner, shortening the average inner loop.
 *
 * <p>Output is identical to {@link GridSoaNumericalSurface} / {@link FasterNumericalSurface}
 * bit-for-bit: a point is buried iff ANY neighbor buries it, so reordering the neighbor scan cannot
 * change which points survive, their coordinates, or the areas - only how early the scan stops.
 *
 * <p>The sort is in-place over the four parallel scratch arrays and reads no instance state, so it
 * is safe to invoke from the superclass constructor (via the {@code orderNeighbors} hook).
 */
public class OrderedGridSoaNumericalSurface extends GridSoaNumericalSurface {

    private static final int INSERTION_CUTOFF = 16;

    public OrderedGridSoaNumericalSurface(IAtomContainer atomContainer) {
        super(atomContainer);
    }

    public OrderedGridSoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        super(atomContainer, solventRadius, tesslevel);
    }

    @Override
    protected void orderNeighbors(double[] diffX, double[] diffY, double[] diffZ, double[] thresh, int numNeighbors) {
        qsort(diffX, diffY, diffZ, thresh, 0, numNeighbors - 1);
    }

    /** In-place quicksort of the four parallel arrays, keyed ascending by {@code th}. Stateless. */
    private static void qsort(double[] dX, double[] dY, double[] dZ, double[] th, int lo, int hi) {
        while (lo < hi) {
            if (hi - lo < INSERTION_CUTOFF) {
                insertionSort(dX, dY, dZ, th, lo, hi);
                return;
            }
            double pivot = th[lo + ((hi - lo) >> 1)];   // middle element (robust on sorted input)
            int i = lo, j = hi;
            while (i <= j) {
                while (th[i] < pivot) i++;
                while (th[j] > pivot) j--;
                if (i <= j) {
                    swap(dX, dY, dZ, th, i, j);
                    i++; j--;
                }
            }
            // recurse into the smaller partition, loop on the larger (bounded stack depth)
            if (j - lo < hi - i) {
                qsort(dX, dY, dZ, th, lo, j);
                lo = i;
            } else {
                qsort(dX, dY, dZ, th, i, hi);
                hi = j;
            }
        }
    }

    private static void insertionSort(double[] dX, double[] dY, double[] dZ, double[] th, int lo, int hi) {
        for (int i = lo + 1; i <= hi; i++) {
            double t = th[i], x = dX[i], y = dY[i], z = dZ[i];
            int j = i - 1;
            while (j >= lo && th[j] > t) {
                th[j + 1] = th[j]; dX[j + 1] = dX[j]; dY[j + 1] = dY[j]; dZ[j + 1] = dZ[j];
                j--;
            }
            th[j + 1] = t; dX[j + 1] = x; dY[j + 1] = y; dZ[j + 1] = z;
        }
    }

    private static void swap(double[] dX, double[] dY, double[] dZ, double[] th, int i, int j) {
        double t = th[i]; th[i] = th[j]; th[j] = t;
        double x = dX[i]; dX[i] = dX[j]; dX[j] = x;
        double y = dY[i]; dY[i] = dY[j]; dY[j] = y;
        double z = dZ[i]; dZ[i] = dZ[j]; dZ[j] = z;
    }
}
