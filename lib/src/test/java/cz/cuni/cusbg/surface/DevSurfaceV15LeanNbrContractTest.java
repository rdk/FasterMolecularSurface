package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV15LeanNbr}. Confirms that the
 * two-pass low-allocation pruned neighbor build and the engine-side cached VdW lookup reproduce the
 * reference result. Runs on whichever scan path the test JVM selected.
 */
class DevSurfaceV15LeanNbrContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV15LeanNbr(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV15LeanNbr(mol, solventRadius, tessLevel);
    }
}
