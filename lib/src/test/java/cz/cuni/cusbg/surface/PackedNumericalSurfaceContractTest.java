package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Runs the full {@link MolecularSurface} contract battery (incl. golden baseline) against
 * {@link PackedNumericalSurface}, the production flat-store + zero-copy delivery implementation.
 */
class PackedNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new PackedNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new PackedNumericalSurface(mol, solventRadius, tessLevel);
    }
}
