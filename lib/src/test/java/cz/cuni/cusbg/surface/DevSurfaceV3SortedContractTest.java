package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV3Sorted}
 * (SoA + flat cell grid + neighbor ordering). The exact point/area golden values confirm that the
 * neighbor reordering does not change the result.
 */
class DevSurfaceV3SortedContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV3Sorted(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV3Sorted(mol, solventRadius, tessLevel);
    }
}
