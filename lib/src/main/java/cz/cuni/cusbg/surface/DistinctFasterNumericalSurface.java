package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.surface.Tessellate;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Distinct-point variant of {@link FasterNumericalSurface}.
 *
 * <p>Same exact pipeline as {@link FasterNumericalSurface} (icosahedral {@code Tessellate}, the hash-box
 * {@link NeighborList}, the scalar occlusion test, and the same van der Waals radii via
 * {@link FasterNumericalSurface#getVdwRadius}), but it emits ONE point per surviving <em>distinct</em>
 * direction instead of the ~5.7x exact-coincident duplicates the tessellation produces (it shares
 * triangle vertices, so each direction appears with multiplicity >=5). The per-atom area is computed
 * from the multiplicity-weighted surviving count, so it is <b>bit-for-bit identical</b> to
 * {@link FasterNumericalSurface}; only the redundant coincident points are removed, so no downstream
 * sparsification is needed.
 *
 * <p>This is the {@code FasterNumericalSurface}-based counterpart of {@link DistinctPackedNumericalSurface}
 * (which builds on the faster {@link PackedNumericalSurface} engine pipeline). The two produce the
 * identical surviving point set and identical areas; they differ only in how fast they compute it.
 *
 * <p>NOTE: deliberately NOT point-set-identical to CDK {@code NumericalSurface} - it omits the
 * coincident duplicates; areas are exact.
 */
public class DistinctFasterNumericalSurface implements MolecularSurface {

    private double solventRadius = 1.4;
    private int tesslevel = 4;
    private final IAtom[] atoms;
    private List<Point3d>[] surfPoints;
    private double[] areas;

    public DistinctFasterNumericalSurface(IAtomContainer atomContainer) {
        this.atoms = AtomContainerManipulator.getAtomArray(atomContainer);
        init();
    }

    public DistinctFasterNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        this.solventRadius = solventRadius;
        this.tesslevel = tesslevel;
        this.atoms = AtomContainerManipulator.getAtomArray(atomContainer);
        init();
    }

    @SuppressWarnings("unchecked")
    private void init() {
        for (IAtom atom : atoms) {
            if (atom.getPoint3d() == null)
                throw new IllegalArgumentException("One or more atoms had no 3D coordinate set");
        }

        double maxRadius = 0;
        for (IAtom atom : atoms) {
            double r = FasterNumericalSurface.getVdwRadius(atom) + solventRadius;
            if (r > maxRadius) maxRadius = r;
        }

        Tessellate tess = new Tessellate("ico", tesslevel);
        tess.doTessellate();
        Point3d[] tessPoints = tess.getTessAsPoint3ds();
        int pointDensity = tess.getNumberOfTriangles() * 3;   // = tessPoints.length (full multiplicity)

        // collapse the tessellation into distinct directions + multiplicity (exact-coincident dedup);
        // all copies of a direction are bit-equal, so an exact-value key groups them.
        int numTess = tessPoints.length;
        Map<String, Integer> index = new HashMap<>(numTess * 2);
        double[] ddx = new double[numTess], ddy = new double[numTess], ddz = new double[numTess];
        int[] mult = new int[numTess];
        int numDir = 0;
        for (Point3d p : tessPoints) {
            String key = p.x + "," + p.y + "," + p.z;
            Integer idx = index.get(key);
            if (idx == null) {
                idx = numDir;
                index.put(key, idx);
                ddx[numDir] = p.x; ddy[numDir] = p.y; ddz[numDir] = p.z;
                numDir++;
            }
            mult[idx]++;
        }

        NeighborList nbrlist = new NeighborList(atoms, maxRadius + solventRadius);

        this.surfPoints = (List<Point3d>[]) new List[atoms.length];
        this.areas = new double[atoms.length];

        IntArrayList neighbors = new IntArrayList(32);
        double[] diffX = new double[512], diffY = new double[512], diffZ = new double[512], thresh = new double[512];

        for (int i = 0; i < atoms.length; i++) {
            IAtom atom = atoms[i];
            double totalRadius = FasterNumericalSurface.getVdwRadius(atom) + solventRadius;
            double totalRadius2 = totalRadius * totalRadius;
            double twiceTotalRadius = 2 * totalRadius;
            Point3d atomPoint = atom.getPoint3d();

            neighbors.clear();
            nbrlist.getNeighborsInto(i, neighbors);
            int numNeighbors = neighbors.size();
            if (diffX.length < numNeighbors) {
                diffX = new double[numNeighbors]; diffY = new double[numNeighbors];
                diffZ = new double[numNeighbors]; thresh = new double[numNeighbors];
            }
            int[] nb = neighbors.buffer;
            for (int k = 0; k < numNeighbors; k++) {
                IAtom neighborAtom = atoms[nb[k]];
                Point3d np = neighborAtom.getPoint3d();
                double x12 = np.x - atomPoint.x, y12 = np.y - atomPoint.y, z12 = np.z - atomPoint.z;
                double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                double tmp = FasterNumericalSurface.getVdwRadius(neighborAtom) + solventRadius;
                tmp = tmp * tmp;
                diffX[k] = x12; diffY[k] = y12; diffZ[k] = z12;
                thresh[k] = (d2 + totalRadius2 - tmp) / twiceTotalRadius;
            }

            // one point per surviving distinct direction; area weighted by the direction multiplicity
            List<Point3d> points = new ArrayList<>();
            long survivingFull = 0;
            for (int d = 0; d < numDir; d++) {
                double px = ddx[d], py = ddy[d], pz = ddz[d];
                boolean buried = false;
                for (int k = 0; k < numNeighbors; k++) {
                    if (diffX[k] * px + diffY[k] * py + diffZ[k] * pz > thresh[k]) { buried = true; break; }
                }
                if (!buried) {
                    points.add(new Point3d(totalRadius * px + atomPoint.x, totalRadius * py + atomPoint.y, totalRadius * pz + atomPoint.z));
                    survivingFull += mult[d];
                }
            }
            this.surfPoints[i] = points;
            this.areas[i] = 4 * Math.PI * totalRadius2 * survivingFull / pointDensity;
        }
    }

    // --- MolecularSurface API (same semantics as FasterNumericalSurface) ---

    @Override
    public Point3d[] getAllSurfacePoints() {
        int npt = 0;
        for (List<Point3d> s : surfPoints) npt += s.size();
        Point3d[] ret = new Point3d[npt];
        int j = 0;
        for (List<Point3d> points : surfPoints)
            for (Point3d p : points) ret[j++] = p;
        return ret;
    }

    @Override
    public Map<IAtom, List<Point3d>> getAtomSurfaceMap() {
        Map<IAtom, List<Point3d>> map = new HashMap<>();
        for (int i = 0; i < surfPoints.length; i++)
            if (!surfPoints[i].isEmpty()) map.put(atoms[i], Collections.unmodifiableList(surfPoints[i]));
        return map;
    }

    @Override
    public Point3d[] getSurfacePoints(int atomIdx) throws CDKException {
        if (atomIdx < 0 || atomIdx >= surfPoints.length) throw new CDKException("Atom index was out of bounds");
        return surfPoints[atomIdx].toArray(new Point3d[0]);
    }

    @Override
    public double getSurfaceArea(int atomIdx) throws CDKException {
        if (atomIdx < 0 || atomIdx >= surfPoints.length) throw new CDKException("Atom index was out of bounds");
        return areas[atomIdx];
    }

    @Override public double[] getAllSurfaceAreas() { return areas; }

    @Override
    public double getTotalSurfaceArea() {
        double ta = 0.0;
        for (double a : areas) ta += a;
        return ta;
    }
}
