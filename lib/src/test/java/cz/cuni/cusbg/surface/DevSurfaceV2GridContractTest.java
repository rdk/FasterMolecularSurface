package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Runs the full {@link MolecularSurface} contract battery (incl. the golden baseline) against
 * {@link DevSurfaceV2Grid} (SoA + flat cell-grid neighbor index).
 */
class DevSurfaceV2GridContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV2Grid(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV2Grid(mol, solventRadius, tessLevel);
    }
}
