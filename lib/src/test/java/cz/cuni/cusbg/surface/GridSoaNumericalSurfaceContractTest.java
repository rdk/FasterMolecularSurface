package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Runs the full {@link MolecularSurface} contract battery (incl. the golden baseline) against
 * {@link GridSoaNumericalSurface} (SoA + flat cell-grid neighbor index).
 */
class GridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new GridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new GridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
