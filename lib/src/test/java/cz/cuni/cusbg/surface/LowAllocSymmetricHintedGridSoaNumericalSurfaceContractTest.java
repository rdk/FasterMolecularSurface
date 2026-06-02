package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link LowAllocSymmetricHintedGridSoaNumericalSurface}
 * (the low-allocation two-pass symmetric precompute). The exact point/area golden values confirm that
 * building the CSR in two passes does not change the result.
 */
class LowAllocSymmetricHintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new LowAllocSymmetricHintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new LowAllocSymmetricHintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
