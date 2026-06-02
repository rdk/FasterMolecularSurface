package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link VectorizedSymmetricHintedGridSoaNumericalSurface}.
 * The exact point/area golden values confirm that the SIMD scan reproduces the scalar result. Runs on
 * whichever scan path the test JVM selected (vectorized when {@code jdk.incubator.vector} is on the
 * module graph, scalar fallback otherwise); the contract holds either way.
 */
class VectorizedSymmetricHintedGridSoaNumericalSurfaceContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new VectorizedSymmetricHintedGridSoaNumericalSurface(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new VectorizedSymmetricHintedGridSoaNumericalSurface(mol, solventRadius, tessLevel);
    }
}
