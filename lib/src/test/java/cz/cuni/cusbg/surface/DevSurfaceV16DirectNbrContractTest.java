package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV16DirectNbr}. Confirms that
 * copy-free CSR neighbor access and the engine-side cached VdW lookup reproduce the reference result.
 * Runs on whichever scan path the test JVM selected.
 */
class DevSurfaceV16DirectNbrContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV16DirectNbr(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV16DirectNbr(mol, solventRadius, tessLevel);
    }
}
