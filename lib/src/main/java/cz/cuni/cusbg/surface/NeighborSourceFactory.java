package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtom;

/**
 * Creates a {@link NeighborSource} over the already-extracted coordinate arrays.
 *
 * <p>Supplied to {@link SoaNumericalSurface} through its constructor (rather than via an
 * overridable method called during construction), so the spatial-index choice is fixed before any
 * computation and a subclass cannot accidentally observe its own uninitialized state.
 */
@FunctionalInterface
interface NeighborSourceFactory {
    NeighborSource create(IAtom[] atoms, double[] ax, double[] ay, double[] az, double radius);
}
