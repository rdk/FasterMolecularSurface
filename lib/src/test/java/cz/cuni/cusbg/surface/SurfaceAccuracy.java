package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;
import java.util.HashSet;
import java.util.Set;

/**
 * Tolerance comparison of a candidate surface against an exact oracle, generalizing the bespoke checks
 * in {@code FloatNumericalSurfaceTest} so any non-bit-exact variant can be scored uniformly.
 *
 * <p>For the {@code SAMPLING} fidelity tier (alternative unit-sphere sampling) only
 * {@link Report#areaTotalRelErr()} is meaningful — the point set is intentionally different, so the
 * caller ignores {@link Report#pointSetSymDiffFrac()} there (the framework never requires a sampling
 * surface to match a reference point-for-point).
 */
final class SurfaceAccuracy {

    private SurfaceAccuracy() {}

    /**
     * @param areaTotalRelErr     |candidate − oracle| / oracle for the total area
     * @param areaAtomMaxRelErr   worst per-atom relative area error
     * @param pointSetSymDiffFrac symmetric difference of the distinct point sets, as a fraction of the
     *                            oracle's point count (meaningful only when both sample the same set)
     */
    record Report(double areaTotalRelErr, double areaAtomMaxRelErr, double pointSetSymDiffFrac) {}

    static Report compare(MolecularSurface candidate, MolecularSurface oracle) {
        double ct = candidate.getTotalSurfaceArea(), ot = oracle.getTotalSurfaceArea();
        double areaTotalRelErr = ot == 0 ? 0 : Math.abs(ct - ot) / ot;

        double[] ca = candidate.getAllSurfaceAreas(), oa = oracle.getAllSurfaceAreas();
        double maxRel = 0;
        int n = Math.min(ca.length, oa.length);
        for (int i = 0; i < n; i++) {
            maxRel = Math.max(maxRel, Math.abs(ca[i] - oa[i]) / Math.max(1.0, oa[i]));
        }

        Set<String> cs = keys(candidate.getAllSurfacePoints());
        Set<String> os = keys(oracle.getAllSurfacePoints());
        int diff = 0;
        for (String k : cs) if (!os.contains(k)) diff++;
        for (String k : os) if (!cs.contains(k)) diff++;
        double symDiff = os.isEmpty() ? 0 : (double) diff / os.size();

        return new Report(areaTotalRelErr, maxRel, symDiff);
    }

    private static Set<String> keys(Point3d[] pts) {
        Set<String> set = new HashSet<>(pts.length * 2);
        for (Point3d p : pts) set.add(p.x + "," + p.y + "," + p.z);
        return set;
    }
}
