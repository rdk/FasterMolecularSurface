package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV18SortedCoords}. Confirms that
 * reading candidate coordinates from a cell-sorted interleaved copy (rather than gathering at the atom
 * index) reproduces the reference result. Runs on whichever scan path the test JVM selected.
 */
class DevSurfaceV18SortedCoordsContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV18SortedCoords(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV18SortedCoords(mol, solventRadius, tessLevel);
    }
}
