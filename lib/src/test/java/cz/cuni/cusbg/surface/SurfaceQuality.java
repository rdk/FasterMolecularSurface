package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtom;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Intrinsic, oracle-free quality metrics on an emitted surface point set. Needs NO reference surface,
 * so it works for ANY family — tessellation, distinct, or the planned density samplers — and is the
 * primary quality story for sampling surfaces that have no bit-exact / point-set oracle.
 *
 * <p>Metrics are computed per atom over that atom's own points (via {@link MolecularSurface#getAtomSurfaceMap()})
 * and aggregated. The nearest-neighbour search is O(points²) per atom; call it on small structures /
 * low tessellation (it is a correctness/quality measurement, not a hot path).
 */
final class SurfaceQuality {

    private SurfaceQuality() {}

    /** Two points within this Euclidean distance count as coincident duplicates. */
    private static final double DUP_EPS = 1e-9;
    /** A point is "too close" if its nearest distinct same-atom neighbour is below this fraction of the median. */
    private static final double TOO_CLOSE_FRACTION = 0.5;

    /**
     * @param totalPoints     total emitted points
     * @param duplicateRatio  fraction coincident (within EPS) with another point on the same atom
     *                        (~0 for distinct/sampling surfaces; high for CDK/Packed full multiplicity)
     * @param tooCloseRatio   fraction whose nearest distinct neighbour is below half the global median spacing
     * @param evennessCV      mean over atoms of (stddev/mean of nearest-distinct-neighbour distances); lower = more even
     * @param minMeanNNRatio  mean over atoms of (min NN / mean NN); 1.0 = perfectly even, →0 = clumped
     */
    record Report(long totalPoints, double duplicateRatio, double tooCloseRatio,
                  double evennessCV, double minMeanNNRatio) {}

    static Report measure(MolecularSurface s) {
        Map<IAtom, List<Point3d>> map = s.getAtomSurfaceMap();
        long total = 0, dup = 0, tooClose = 0;
        double cvSum = 0, minMeanSum = 0;
        int atomsScored = 0;
        List<double[]> perAtomNN = new ArrayList<>();

        for (List<Point3d> pts : map.values()) {
            int n = pts.size();
            total += n;
            if (n < 2) continue;
            double[] nn = new double[n];   // nearest *distinct* neighbour distance per point
            for (int i = 0; i < n; i++) {
                Point3d pi = pts.get(i);
                double best = Double.POSITIVE_INFINITY;
                boolean coincident = false;
                for (int j = 0; j < n; j++) {
                    if (i == j) continue;
                    double d = pi.distance(pts.get(j));
                    if (d <= DUP_EPS) { coincident = true; continue; }
                    if (d < best) best = d;
                }
                if (coincident) dup++;
                nn[i] = best;
            }
            perAtomNN.add(nn);

            double[] finite = Arrays.stream(nn).filter(Double::isFinite).toArray();
            if (finite.length >= 2) {
                double mean = Arrays.stream(finite).average().orElse(0);
                if (mean > 0) {
                    double var = 0;
                    for (double d : finite) var += (d - mean) * (d - mean);
                    double sd = Math.sqrt(var / finite.length);
                    double min = Arrays.stream(finite).min().orElse(0);
                    cvSum += sd / mean;
                    minMeanSum += min / mean;
                    atomsScored++;
                }
            }
        }

        // global median nearest-neighbour distance drives the "too close" threshold
        double[] allNN = perAtomNN.stream().flatMapToDouble(Arrays::stream).filter(Double::isFinite).sorted().toArray();
        double median = allNN.length == 0 ? 0 : allNN[allNN.length / 2];
        double thresh = TOO_CLOSE_FRACTION * median;
        for (double[] nn : perAtomNN)
            for (double d : nn)
                if (Double.isFinite(d) && d < thresh) tooClose++;

        return new Report(
                total,
                total == 0 ? 0 : (double) dup / total,
                total == 0 ? 0 : (double) tooClose / total,
                atomsScored == 0 ? 0 : cvSum / atomsScored,
                atomsScored == 0 ? 0 : minMeanSum / atomsScored);
    }
}
