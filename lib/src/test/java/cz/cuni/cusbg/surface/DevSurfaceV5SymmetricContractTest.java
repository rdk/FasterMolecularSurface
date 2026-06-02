package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV5Symmetric}
 * (SoA + symmetric-precompute cell grid + last-occluder-first scan, no sort). The exact point/area
 * golden values confirm that computing each neighbor pair once does not change the result.
 */
class DevSurfaceV5SymmetricContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV5Symmetric(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV5Symmetric(mol, solventRadius, tessLevel);
    }
}
