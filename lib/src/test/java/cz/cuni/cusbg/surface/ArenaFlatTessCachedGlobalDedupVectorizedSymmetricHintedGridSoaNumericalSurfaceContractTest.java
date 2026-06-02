package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface}. Building many surfaces of varying
 * sizes through these tests on one thread also exercises the per-thread scratch arena's cross-build
 * reuse, confirming reused buffers never leak stale values. Runs on whichever scan path the test JVM
 * selected.
 */
class ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new ArenaFlatTessCachedGlobalDedupVectorizedSymmetricHintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
