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

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;

import javax.vecmath.Point3d;
import java.util.List;
import java.util.Map;

/**
 * Common contract for numerical solvent-accessible-surface implementations.
 *
 * <p>This interface captures the public API shared by {@link FasterNumericalSurface}
 * and any future, differently-optimized variant. It mirrors the surface-querying
 * methods of CDK's {@code org.openscience.cdk.geometry.surface.NumericalSurface}, so the
 * same battery of contract / equivalence / golden tests can be run against every variant
 * simply by supplying a factory that produces a {@code MolecularSurface}.
 *
 * <p>Implementations are expected to compute the surface eagerly (typically in their
 * constructor); all accessor methods below return already-computed results.
 */
public interface MolecularSurface {

    /**
     * Get an array of all the points on the molecular surface.
     *
     * @return an array of {@link Point3d} objects representing all points on the surface
     */
    Point3d[] getAllSurfacePoints();

    /**
     * Get the map from atom to surface points. If an atom does not appear in the map it is
     * buried. Atoms may share surface points with other atoms.
     *
     * @return surface atoms and their associated points on the surface
     */
    Map<IAtom, List<Point3d>> getAtomSurfaceMap();

    /**
     * Get an array of the points on the accessible surface of a specific atom.
     *
     * @param atomIdx the index of the atom (0 .. n-1, n = number of atoms)
     * @return an array of {@link Point3d} objects
     * @throws CDKException if the atom index is outside the range of allowable indices
     */
    Point3d[] getSurfacePoints(int atomIdx) throws CDKException;

    /**
     * Get the accessible surface area for the specified atom.
     *
     * @param atomIdx the index of the atom (0 .. n-1, n = number of atoms)
     * @return the accessible surface area of the atom
     * @throws CDKException if the atom index is outside the range of allowable indices
     */
    double getSurfaceArea(int atomIdx) throws CDKException;

    /**
     * Get an array containing the accessible surface area for each atom.
     *
     * @return an array giving the surface areas of all atoms
     */
    double[] getAllSurfaceAreas();

    /**
     * Get the total surface area for the molecule.
     *
     * @return the total accessible surface area
     */
    double getTotalSurfaceArea();
}
