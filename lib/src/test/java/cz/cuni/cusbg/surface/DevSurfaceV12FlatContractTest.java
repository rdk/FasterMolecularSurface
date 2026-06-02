package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for
 * {@link DevSurfaceV12Flat}. Exercises the flat
 * {@code double[]} point store (and its lazy {@code Point3d} materialization in the accessors) and the
 * cached VdW lookup, confirming both reproduce the reference result. Runs on whichever scan path the
 * test JVM selected.
 */
class DevSurfaceV12FlatContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV12Flat(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV12Flat(mol, solventRadius, tessLevel);
    }
}
