package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link OrderedGridSoaNumericalSurface}
 * (SoA + flat cell grid + neighbor ordering). The exact point/area golden values confirm that the
 * neighbor reordering does not change the result.
 */
class OrderedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new OrderedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new OrderedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
