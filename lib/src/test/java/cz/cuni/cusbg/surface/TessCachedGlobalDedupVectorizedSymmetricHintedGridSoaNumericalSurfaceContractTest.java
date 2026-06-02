package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for
 * {@link TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface}. Confirms that reusing
 * a process-wide cached tessellation (rather than rebuilding it per surface) reproduces the reference
 * result. Runs on whichever scan path the test JVM selected.
 */
class TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new TessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
