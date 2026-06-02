package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV19FlatStore}. Confirms that V18
 * with the flat point store (lazily materializing {@code Point3d}) reproduces the reference result.
 */
class DevSurfaceV19FlatStoreContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV19FlatStore(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV19FlatStore(mol, solventRadius, tessLevel);
    }
}
