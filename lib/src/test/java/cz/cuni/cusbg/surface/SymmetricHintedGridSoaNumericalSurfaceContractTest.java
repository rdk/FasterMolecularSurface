package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link SymmetricHintedGridSoaNumericalSurface}
 * (SoA + symmetric-precompute cell grid + last-occluder-first scan, no sort). The exact point/area
 * golden values confirm that computing each neighbor pair once does not change the result.
 */
class SymmetricHintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new SymmetricHintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new SymmetricHintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
