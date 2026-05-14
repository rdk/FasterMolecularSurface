/*
 * Copyright (C) 2004-2007  The Chemistry Development Kit (CDK) project
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.geometry.surface.Tessellate;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.vecmath.Point3d;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.*;

/**
 * A class representing the solvent accessible surface area surface of a molecule.
 *
 * <p>This class is based on the Python implementation of the DCLM method
 * ({@cdk.cite EIS95}) by Peter McCluskey, which is a non-analytical method to generate a set of points
 * representing the solvent accessible surface area of a molecule.
 *
 * <p>The neighbor list is a simplified version of that
 * described in {@cdk.cite EIS95} and as a result, the surface areas of the atoms may not be exact
 * (compared to analytical calculations). The tessellation is slightly different from
 * that described by McCluskey and uses recursive subdivision starting from an icosahedral
 * representation.
 *
 * <p>The default solvent radius used is 1.4A and setting this to 0 will give the
 * Van der Waals surface. The accuracy can be increased by increasing the tessellation
 * level, though the default of 4 is a good balance between accuracy and speed.
 *
 * @author      Rajarshi Guha
 * @cdk.created 2005-05-08
 * @cdk.module  qsarmolecular
 * @cdk.githash
 * @cdk.bug     1846421
 */
public class FasterNumericalSurface {

    private static ILoggingTool logger         = LoggingToolFactory.createLoggingTool(NumericalSurface.class);
    private static final Logger log = LoggerFactory.getLogger(FasterNumericalSurface.class);

    /**
     * Fallback VdW radii indexed by element symbol.
     *
     * CDK's {@link org.openscience.cdk.config.Elements} enum hard-codes the VdW radius as {@code null}
     * for several biologically common elements (Co, Ni, Cu, Rh, Os, Ir, several lanthanides),
     * even though CDK ships the values in {@code radii-vdw.txt}. Calling
     * {@link PeriodicTable#getVdwRadius(String)} on those returns {@code null} and any
     * downstream {@code rw + solventRadius} arithmetic throws NPE.
     *
     * This table is populated from a copy of CDK's own {@code radii-vdw.txt} bundled as a
     * resource, so the fallback values match CDK's data file exactly.
     */
    private static final Map<String, Double> VDW_FALLBACK_BY_SYMBOL = loadVdwFallbackTable();

    /** Last-resort default when neither CDK enum nor radii-vdw.txt has a value. Matches the
     * comment "Missing values are represented by 2 as a mean value" at the top of radii-vdw.txt. */
    private static final double VDW_DEFAULT = 2.0;

    private static Map<String, Double> loadVdwFallbackTable() {
        Map<String, Double> map = new HashMap<>();
        try (InputStream in = FasterNumericalSurface.class.getResourceAsStream("radii-vdw.txt")) {
            if (in == null) {
                log.warn("Bundled radii-vdw.txt not found on classpath; falling back to {} A default for unknown elements", VDW_DEFAULT);
                return map;
            }
            try (BufferedReader r = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
                String line;
                while ((line = r.readLine()) != null) {
                    String trimmed = line.trim();
                    if (trimmed.isEmpty() || trimmed.startsWith("#")) continue;
                    String[] parts = trimmed.split("\\s+");
                    if (parts.length < 2) continue;
                    try {
                        int z = Integer.parseInt(parts[0]);
                        double rw = Double.parseDouble(parts[1]);
                        String sym = PeriodicTable.getSymbol(z);
                        if (sym != null && !sym.isEmpty()) {
                            map.put(sym, rw);
                        }
                    } catch (NumberFormatException ignored) {
                        // skip malformed lines
                    }
                }
            }
        } catch (IOException e) {
            log.warn("Failed to load bundled radii-vdw.txt: {}", e.toString());
        }
        return map;
    }

    /**
     * Get the VdW radius for an atom with graceful fallback.
     *
     * Order of preference: CDK's {@link PeriodicTable#getVdwRadius} (enum-backed),
     * then the bundled {@code radii-vdw.txt} table, then a 2.0 A default.
     */
    private static double getVdwRadius(IAtom atom) {
        String symbol = atom.getSymbol();
        Double r = PeriodicTable.getVdwRadius(symbol);
        if (r != null) return r;
        Double fallback = VDW_FALLBACK_BY_SYMBOL.get(symbol);
        if (fallback != null) {
            log.debug("CDK enum has no VdW radius for {} - using {} A from radii-vdw.txt", symbol, fallback);
            return fallback;
        }
        log.warn("No VdW radius for element {} - using {} A default", symbol, VDW_DEFAULT);
        return VDW_DEFAULT;
    }

    double          solventRadius  = 1.4;
    int             tesslevel      = 4;
    IAtom[]         atoms;
    List<Point3d>[] surfPoints;
    double[]        areas;
    double[]        volumes;

    /**
     * Constructor to initialize the surface calculation with default values.
     *
     * This constructor use the Van der Waals radii as defined in <i>org/openscience/cdk/config/data/jmol_atomtypes.txt</i>
     * of the source distribution. Also uses a tesselation level of 4 and solvent radius of 1.4A.
     *
     * @param atomContainer The {@link IAtomContainer} for which the surface is to be calculated
     */
    public FasterNumericalSurface(IAtomContainer atomContainer) {
        this.atoms = AtomContainerManipulator.getAtomArray(atomContainer);
        init();
    }

    /**
     * Constructor to initialize the surface calculation with user specified values.
     *
     * This constructor use the Van der Waals radii as defined in <i>org/openscience/cdk/config/data/jmol_atomtypes.txt</i>
     * of the source distribution
     *
     * @param atomContainer The {@link IAtomContainer} for which the surface is to be calculated
     * @param solventRadius The radius of a solvent molecule that is used to extend
     * the radius of each atom. Setting to 0 gives the Van der Waals surface
     * @param tesslevel The number of levels that the subdivision algorithm for tessellation
     * should use
     */
    public FasterNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        this.solventRadius = solventRadius;
        this.atoms = AtomContainerManipulator.getAtomArray(atomContainer);
        this.tesslevel = tesslevel;
        init();
    }

    /**
     * Evaluate the surface.
     *
     * This method generates the points on the accessible surface area of each atom
     * as well as calculating the surface area of each atom.
     * @deprecated
     */
    @Deprecated
    public void calculateSurface() {
        // NO-OP
    }

    /**
     * Initialize the surface, generating the points on the accessible surface
     * area of each atom as well as calculating the surface area of each atom.
     */
    private void init() {

        // invariants
        for (IAtom atom : atoms) {
            if (atom.getPoint3d() == null)
                throw new IllegalArgumentException("One or more atoms had no 3D coordinate set");
        }

        // get r_f and geometric center
        Point3d cp = new Point3d(0, 0, 0);
        double maxRadius = 0;
        for (IAtom atom : atoms) {
            double vdwr = getVdwRadius(atom);
            if (vdwr + solventRadius > maxRadius)
                maxRadius = vdwr + solventRadius;

            cp.x = cp.x + atom.getPoint3d().x;
            cp.y = cp.y + atom.getPoint3d().y;
            cp.z = cp.z + atom.getPoint3d().z;
        }
        cp.x = cp.x / atoms.length;
        cp.y = cp.y / atoms.length;
        cp.z = cp.z / atoms.length;

        // do the tesselation
        Tessellate tess = new Tessellate("ico", tesslevel);
        tess.doTessellate();
        logger.info("Got tesselation, number of triangles = " + tess.getNumberOfTriangles());
        Point3d[] tessPoints = tess.getTessAsPoint3ds();

        // get neighbor list
        NeighborList nbrlist = new NeighborList(atoms, maxRadius + solventRadius);
        logger.info("Got neighbor list");

        /*
         * for (int i = 0; i < atoms.length; i++) { int[] nlist =
         * nbrlist.getNeighbors(i); logger.debug("Atom "+i+": "); for (int j =
         * 0; j < nlist.length; j++) logger.debug(j+" "); logger.debug(""); }
         */

        // loop over atoms and get surface points
        this.surfPoints = (List<Point3d>[]) new List[atoms.length];
        this.areas = new double[atoms.length];
        this.volumes = new double[atoms.length];

        List<Point3d[]> reusedPointList = new ArrayList<>(512);
        PointDiffData[] reusedDataBuffer = createEmptyDiffData(512);  // required length is the max number of neighbors of individual atoms
        IntArrayList neighborBuffer = new IntArrayList(32);
        int pointDensity = tess.getNumberOfTriangles() * 3;
        
        for (int i = 0; i < atoms.length; i++) {
            List<Point3d[]> points = atomicSurfacePoints(nbrlist, i, atoms[i], tessPoints, reusedPointList, reusedDataBuffer, neighborBuffer);
            translatePoints(i, points, pointDensity, atoms[i], cp);
        }
        logger.info("Obtained points, areas and volumes");

    }

    /**
     * Get an array of all the points on the molecular surface.
     *
     * This returns an array of Point3d objects representing all the points
     * on the molecular surface
     *
     * @return  An array of Point3d objects
     */
    public Point3d[] getAllSurfacePoints() {
        int npt = 0;
        for (List<Point3d> surfPoint : this.surfPoints)
            npt += surfPoint.size();
        Point3d[] ret = new Point3d[npt];
        int j = 0;
        for (List<Point3d> points : this.surfPoints) {
            for (Point3d p : points) {
                ret[j++] = p;
            }
        }
        return (ret);
    }

    /**
     * Get the map from atom to surface points. If an atom does not appear in
     * the map it is buried. Atoms may share surface points with other atoms.
     *
     * @return surface atoms and associated points on the surface
     */
    public Map<IAtom, List<Point3d>> getAtomSurfaceMap() {
        Map<IAtom,List<Point3d>> map = new HashMap<>();
        for (int i = 0; i < this.surfPoints.length; i++) {
            if (!this.surfPoints[i].isEmpty())
                map.put(this.atoms[i], Collections.unmodifiableList(this.surfPoints[i]));
        }
        return map;
    }

    /**
     * Get an array of the points on the accessible surface of a specific atom.
     *
     * @param atomIdx The index of the atom. Ranges from 0 to n-1, where n is the
     * number of atoms in the AtomContainer that the surface was calculated for
     * @return  An array of Point3d objects
     * @throws CDKException if the atom index is outside the range of allowable indices
     */
    public Point3d[] getSurfacePoints(int atomIdx) throws CDKException {
        if (atomIdx >= this.surfPoints.length)
            throw new CDKException("Atom index was out of bounds");
        return this.surfPoints[atomIdx].toArray(new Point3d[0]);
    }

    /**
     * Get the surface area for the specified atom.
     *
     * @param atomIdx The index of the atom. Ranges from 0 to n-1, where n is the
     * number of atoms in the AtomContainer that the surface was calculated for
     * @return A double representing the accessible surface area of the atom
     * @throws CDKException if the atom index is outside the range of allowable indices
     */
    public double getSurfaceArea(int atomIdx) throws CDKException {
        if (atomIdx >= this.surfPoints.length)
            throw new CDKException("Atom index was out of bounds");
        return this.areas[atomIdx];
    }

    /**
     * Get an array containing the accessible surface area for each atom.
     *
     * @return An array of double giving the surface areas of all the atoms
     */
    public double[] getAllSurfaceAreas() {
        return this.areas;
    }

    /**
     * Get the total surface area for the AtomContainer.
     *
     * @return A double containing the total surface area of the AtomContainer for
     * which the surface was calculated for
     */
    public double getTotalSurfaceArea() {
        double ta = 0.0;
        for (double area : this.areas)
            ta += area;
        return ta;
    }

    private void translatePoints(int atmIdx, List<Point3d[]> points, int pointDensity, IAtom atom, Point3d cp) {
        double totalRadius = getVdwRadius(atom) + solventRadius;

        double area = 4 * Math.PI * (totalRadius * totalRadius) * points.size() / pointDensity;

        double sumx = 0.0;
        double sumy = 0.0;
        double sumz = 0.0;
        for (Point3d[] point : points) {
            Point3d p = point[1];
            sumx += p.x;
            sumy += p.y;
            sumz += p.z;
        }
        double vconst = 4.0 / 3.0 * Math.PI / (double) pointDensity;
        double dotp1 = (atom.getPoint3d().x - cp.x) * sumx + (atom.getPoint3d().y - cp.y) * sumy
                + (atom.getPoint3d().z - cp.z) * sumz;
        double volume = vconst * (totalRadius * totalRadius) * dotp1 + (totalRadius * totalRadius * totalRadius)
                * points.size();

        this.areas[atmIdx] = area;
        this.volumes[atmIdx] = volume;

        List<Point3d> surfPoints = new ArrayList<>(points.size());
        for (Point3d[] point : points)
            surfPoints.add(point[0]);
        this.surfPoints[atmIdx] = surfPoints;
    }

    private static class PointDiffData {
        double diffX;
        double diffY;
        double diffZ;
        double thresh;
    }

    PointDiffData[] createEmptyDiffData(int size) {
        PointDiffData[] res = new PointDiffData[size];
        for (int i = 0; i != size; i++) {
            res[i] = new PointDiffData();
        }
        return res;
    }

    private List<Point3d[]> atomicSurfacePoints(NeighborList nbrlist, int currAtomIdx, IAtom atom, Point3d[] tessPoints,
                                                List<Point3d[]> reusedPointList, PointDiffData[] reusedDataBuffer, IntArrayList neighbors) {

        double totalRadius = getVdwRadius(atom) + solventRadius;
        double totalRadius2 = totalRadius * totalRadius;
        double twiceTotalRadius = 2 * totalRadius;

        Point3d atomPoint = atom.getPoint3d();

        neighbors.clear();  // reuse buffer
        nbrlist.getNeighborsInto(currAtomIdx, neighbors);
        int numNeighbors = neighbors.size();

        PointDiffData[] data = reusedDataBuffer;
        if (data.length < numNeighbors) {
            log.info("Must allocate array with n={}", numNeighbors);
            data = createEmptyDiffData(numNeighbors);
        }

        int[] neighBuff = neighbors.buffer;
        for (int i = 0; i != numNeighbors; i++) {
            IAtom neighborAtom = atoms[neighBuff[i]];
            Point3d neighborPoint = neighborAtom.getPoint3d();

            double x12 = neighborPoint.x - atomPoint.x;
            double y12 = neighborPoint.y - atomPoint.y;
            double z12 = neighborPoint.z - atomPoint.z;

            double d2 = x12 * x12 + y12 * y12 + z12 * z12;
            double tmp = getVdwRadius(neighborAtom) + solventRadius;
            tmp = tmp * tmp;
            double thresh = (d2 + totalRadius2 - tmp) / twiceTotalRadius;

            PointDiffData data_i = data[i];
            data_i.diffX = x12;
            data_i.diffY = y12;
            data_i.diffZ = z12;
            data_i.thresh = thresh;
        }

        //Point3d[] tessPoints = tess.getTessAsPoint3ds();

        List<Point3d[]> points = reusedPointList;
        points.clear();

        collectPoints(tessPoints, numNeighbors, data, totalRadius, atomPoint, points);

        // the first column contains the transformed points
        // and the second column contains the points from the
        // original unit tesselation
        return points;
    }

    /**
     * takes most of the time
     */
    private static void collectPoints(Point3d[] tessPoints, int numNeighbors, PointDiffData[] data, double totalRadius, Point3d atomPoint, List<Point3d[]> points) {
        boolean buried;
        double px, py, pz;

        for (Point3d pt : tessPoints) {
            px = pt.x;
            py = pt.y;
            pz = pt.z;
            buried = false;

            for (int i = 0; i != numNeighbors; i++) {
                PointDiffData datum = data[i];
                if (datum.diffX * px + datum.diffY * py + datum.diffZ * pz > datum.thresh) {
                    buried = true;
                    break;
                }
            }
            if (!buried) {
                Point3d[] tmp = new Point3d[2];
                tmp[0] = new Point3d(totalRadius * px + atomPoint.x, totalRadius * py
                    + atomPoint.y, totalRadius * pz + atomPoint.z);
                tmp[1] = pt;
                points.add(tmp);
            }
        }
    }

}
