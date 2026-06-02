package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link GlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface}.
 * Confirms that caching the distinct-direction mapping process-wide (rather than rebuilding it per
 * surface) reproduces the reference result. Runs on whichever scan path the test JVM selected.
 */
class GlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new GlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new GlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
