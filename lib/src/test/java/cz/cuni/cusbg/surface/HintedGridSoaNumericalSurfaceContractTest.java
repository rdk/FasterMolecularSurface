package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link HintedGridSoaNumericalSurface}
 * (SoA + flat cell grid + last-occluder-first scan hint, no sort). The exact point/area golden
 * values confirm that the scan-order hint does not change the result.
 */
class HintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new HintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new HintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
