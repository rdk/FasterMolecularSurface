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
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import javax.vecmath.Point3d;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.IntFunction;

/**
 * Structure-of-arrays (SoA) variant of {@link FasterNumericalSurface}.
 *
 * <p>Algorithmically identical to {@link FasterNumericalSurface} (same tessellation, same neighbor
 * sets, same occlusion test and arithmetic, so it reproduces the golden baseline bit-for-bit), but
 * it removes per-element overhead from the two hot loops:
 *
 * <ul>
 *   <li>Atom coordinates and per-atom expanded radii ({@code vdw + solvent}) and their squares are
 *       extracted once into flat {@code double[]}, so the neighbor pre-pass no longer calls
 *       {@code IAtom.getPoint3d()} / {@code getVdwRadius()} per neighbor.</li>
 *   <li>The reused per-neighbor scratch is parallel {@code double[]} arrays (diffX/diffY/diffZ/thresh)
 *       instead of an array of {@code PointDiffData} objects, so the innermost occlusion loop strides
 *       contiguous memory instead of chasing a pointer per neighbor.</li>
 *   <li>Tessellation points are held as flat {@code double[]} (tx/ty/tz).</li>
 * </ul>
 *
 * <p>This is optimization step 1 (SoA); the {@link NeighborList} itself is unchanged. The neighbor
 * index is supplied through a {@link NeighborSourceFactory} constructor argument, so a subclass (e.g.
 * {@link GridSoaNumericalSurface}) can swap in a different index without overriding a method during
 * construction. All scratch derived from the atoms is local to {@link #init()} and not retained.
 */
public class SoaNumericalSurface implements MolecularSurface {

    /** Default neighbor index: the hash-box {@link NeighborList} (identical sets to FasterNumericalSurface). */
    private static final NeighborSourceFactory DEFAULT_NEIGHBORS =
            (atoms, ax, ay, az, radius) -> {
                NeighborList nl = new NeighborList(atoms, radius);
                return nl::getNeighborsInto;
            };

    private final double  solventRadius;
    private final int     tesslevel;
    private final IAtom[] atoms;
    private final NeighborSourceFactory neighborFactory;
    private final NeighborOrdering ordering;
    private final OcclusionScan scan;
    private final TessellationProvider tessProvider;
    private final IntFunction<SurfacePointStore> storeFactory;
    private final boolean useArena;

    /** Per-thread reusable transient scratch, used only when {@code useArena} is set. */
    private static final ThreadLocal<EngineScratch> ARENA = ThreadLocal.withInitial(EngineScratch::new);

    private SurfacePointStore store;
    private double[]          areas;

    public SoaNumericalSurface(IAtomContainer atomContainer) {
        this(atomContainer, 1.4, 4);
    }

    public SoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel) {
        this(atomContainer, solventRadius, tesslevel, DEFAULT_NEIGHBORS, NeighborOrdering.NONE, OcclusionScan.STANDARD);
    }

    /**
     * @param neighborFactory supplies the neighbor index over the extracted coordinate arrays.
     * @param ordering        reorders the per-neighbor scratch before the occlusion loop (a no-op for
     *        the base class; see {@link OrderedGridSoaNumericalSurface}).
     * @param scan            performs the occlusion test over each atom's tessellation points (the
     *        reference loop for the base class; see {@link HintedGridSoaNumericalSurface}). All three
     *        strategies are passed here (not via overridable methods) so they are fixed before
     *        {@link #init()} runs and cannot observe uninitialized subclass state; all must be
     *        stateless for the same reason.
     */
    protected SoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel,
                                  NeighborSourceFactory neighborFactory, NeighborOrdering ordering, OcclusionScan scan) {
        this(atomContainer, solventRadius, tesslevel, neighborFactory, ordering, scan, TessellationProvider.FRESH);
    }

    /**
     * @param tessProvider supplies the unit-sphere tessellation. {@link TessellationProvider#FRESH}
     *        (the default for the other constructors) rebuilds it per surface, matching the original
     *        behavior; a variant may pass {@link TessellationProvider#CACHED} to reuse a process-wide
     *        shared tessellation. Like the other strategies it is fixed before {@link #init()} runs.
     */
    protected SoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel,
                                  NeighborSourceFactory neighborFactory, NeighborOrdering ordering, OcclusionScan scan,
                                  TessellationProvider tessProvider) {
        this(atomContainer, solventRadius, tesslevel, neighborFactory, ordering, scan, tessProvider, ListSurfacePointStore::new);
    }

    /**
     * @param storeFactory builds the per-surface point store from the atom count. The default
     *        {@link ListSurfacePointStore} keeps one {@code Point3d} per point (original behavior); a
     *        variant may pass {@link FlatSurfacePointStore}{@code ::new} to store coordinates flat and
     *        materialize {@code Point3d} lazily. Fixed before {@link #init()} like the other strategies.
     */
    protected SoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel,
                                  NeighborSourceFactory neighborFactory, NeighborOrdering ordering, OcclusionScan scan,
                                  TessellationProvider tessProvider, IntFunction<SurfacePointStore> storeFactory) {
        this(atomContainer, solventRadius, tesslevel, neighborFactory, ordering, scan, tessProvider, storeFactory, false);
    }

    /**
     * @param useArena when true, the engine's transient per-build scratch (coordinate/radius arrays,
     *        neighbor list, diff/thresh) is drawn from a reusable per-thread {@link EngineScratch} arena
     *        instead of allocated per surface (optimization A). Output-identical; only allocation differs.
     *        Requires a neighbor source that tolerates oversized coordinate arrays (takes the atom count
     *        explicitly) - the pruned source does.
     */
    protected SoaNumericalSurface(IAtomContainer atomContainer, double solventRadius, int tesslevel,
                                  NeighborSourceFactory neighborFactory, NeighborOrdering ordering, OcclusionScan scan,
                                  TessellationProvider tessProvider, IntFunction<SurfacePointStore> storeFactory,
                                  boolean useArena) {
        this.solventRadius = solventRadius;
        this.tesslevel = tesslevel;
        this.atoms = AtomContainerManipulator.getAtomArray(atomContainer);
        this.neighborFactory = neighborFactory;
        this.ordering = ordering;
        this.scan = scan;
        this.tessProvider = tessProvider;
        this.storeFactory = storeFactory;
        this.useArena = useArena;
        init();
    }

    private void init() {
        int n = atoms.length;

        // transient scratch: reused from a per-thread arena when enabled, else allocated per build.
        // The coordinate arrays may then be larger than n; only [0, n) is used and the neighbor source
        // takes n explicitly, so the extra capacity is inert.
        EngineScratch sc = useArena ? ARENA.get() : null;
        double[] ax, ay, az, atomRadius, atomRadius2;
        if (sc != null) {
            sc.ensureAtoms(n);
            ax = sc.ax; ay = sc.ay; az = sc.az; atomRadius = sc.atomRadius; atomRadius2 = sc.atomRadius2;
        } else {
            ax = new double[n]; ay = new double[n]; az = new double[n];
            atomRadius = new double[n]; atomRadius2 = new double[n];
        }
        double maxRadius = 0;
        for (int i = 0; i < n; i++) {
            Point3d p = atoms[i].getPoint3d();
            if (p == null) throw new IllegalArgumentException("One or more atoms had no 3D coordinate set");
            ax[i] = p.x; ay[i] = p.y; az[i] = p.z;
            double r = FasterNumericalSurface.getVdwRadius(atoms[i]) + solventRadius;
            atomRadius[i]  = r;
            atomRadius2[i] = r * r;
            if (r > maxRadius) maxRadius = r;
        }

        // tessellation (flat SoA): fresh per build by default; a cached provider shares it across builds
        Tessellation tessel = tessProvider.get(tesslevel);
        double[] tx = tessel.tx, ty = tessel.ty, tz = tessel.tz;
        int numTess = tessel.numTess;
        int pointDensity = tessel.pointDensity;

        NeighborSource neighbors = neighborFactory.create(atoms, ax, ay, az, maxRadius + solventRadius);

        this.store = storeFactory.apply(n);
        this.areas = new double[n];

        // reused buffers (from the arena when enabled)
        IntArrayList nbr = sc != null ? sc.nbr : new IntArrayList(32);
        double[] diffX, diffY, diffZ, thresh;
        if (sc != null) { diffX = sc.diffX; diffY = sc.diffY; diffZ = sc.diffZ; thresh = sc.thresh; }
        else { diffX = new double[512]; diffY = new double[512]; diffZ = new double[512]; thresh = new double[512]; }

        for (int i = 0; i < n; i++) {
            double totalRadius      = atomRadius[i];
            double totalRadius2     = atomRadius2[i];
            double twiceTotalRadius = 2 * totalRadius;
            double atomX = ax[i], atomY = ay[i], atomZ = az[i];

            nbr.clear();
            neighbors.getNeighborsInto(i, nbr);
            int numNeighbors = nbr.size();

            if (diffX.length < numNeighbors) {
                diffX = new double[numNeighbors]; diffY = new double[numNeighbors];
                diffZ = new double[numNeighbors]; thresh = new double[numNeighbors];
                if (sc != null) { sc.diffX = diffX; sc.diffY = diffY; sc.diffZ = diffZ; sc.thresh = thresh; }
            }

            int[] nb = nbr.buffer;
            for (int k = 0; k < numNeighbors; k++) {
                int j = nb[k];
                double x12 = ax[j] - atomX;
                double y12 = ay[j] - atomY;
                double z12 = az[j] - atomZ;
                double d2 = x12 * x12 + y12 * y12 + z12 * z12;
                // identical arithmetic to FasterNumericalSurface: tmp = (vdw(j)+solvent)^2 = atomRadius2[j]
                diffX[k] = x12;
                diffY[k] = y12;
                diffZ[k] = z12;
                thresh[k] = (d2 + totalRadius2 - atomRadius2[j]) / twiceTotalRadius;
            }

            // reorder the per-neighbor scratch so collectPoints' early-exit can fire sooner. A no-op
            // by default; output-preserving regardless (a point is buried iff ANY neighbor buries it,
            // so the order of the OR does not matter). See NeighborOrdering.
            ordering.order(diffX, diffY, diffZ, thresh, numNeighbors);

            store.startAtom(i);
            scan.collect(tx, ty, tz, numTess, numNeighbors, diffX, diffY, diffZ, thresh,
                    totalRadius, atomX, atomY, atomZ, store);
            int count = store.finishAtom();

            this.areas[i] = 4 * Math.PI * (totalRadius * totalRadius) * count / pointDensity;
        }
    }

    // --- MolecularSurface API (identical semantics to FasterNumericalSurface) ---

    @Override
    public Point3d[] getAllSurfacePoints() {
        return store.allPoints();
    }

    @Override
    public Map<IAtom, List<Point3d>> getAtomSurfaceMap() {
        Map<IAtom, List<Point3d>> map = new HashMap<>();
        for (int i = 0; i < atoms.length; i++)
            if (store.atomPointCount(i) > 0)
                map.put(atoms[i], store.atomPoints(i));
        return map;
    }

    @Override
    public Point3d[] getSurfacePoints(int atomIdx) throws CDKException {
        if (atomIdx < 0 || atomIdx >= atoms.length)
            throw new CDKException("Atom index was out of bounds");
        return store.atomPoints(atomIdx).toArray(new Point3d[0]);
    }

    @Override
    public double getSurfaceArea(int atomIdx) throws CDKException {
        if (atomIdx < 0 || atomIdx >= atoms.length)
            throw new CDKException("Atom index was out of bounds");
        return areas[atomIdx];
    }

    @Override
    public double[] getAllSurfaceAreas() {
        return areas;
    }

    @Override
    public double getTotalSurfaceArea() {
        double ta = 0.0;
        for (double a : areas) ta += a;
        return ta;
    }
}
