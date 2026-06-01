package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Runs the full {@link MolecularSurface} contract battery against {@link FasterNumericalSurface}.
 *
 * <p>To validate a future variant, create an analogous subclass returning that variant from the
 * two factory methods — it inherits every contract test automatically.
 */
class FasterNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new FasterNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new FasterNumericalSurface(mol, solventRadius, tessLevel);
    }
}
