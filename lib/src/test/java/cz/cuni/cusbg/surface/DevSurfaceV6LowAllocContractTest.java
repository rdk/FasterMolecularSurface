package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV6LowAlloc}
 * (the low-allocation two-pass symmetric precompute). The exact point/area golden values confirm that
 * building the CSR in two passes does not change the result.
 */
class DevSurfaceV6LowAllocContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV6LowAlloc(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV6LowAlloc(mol, solventRadius, tessLevel);
    }
}
