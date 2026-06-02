package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV7Simd}.
 * The exact point/area golden values confirm that the SIMD scan reproduces the scalar result. Runs on
 * whichever scan path the test JVM selected (vectorized when {@code jdk.incubator.vector} is on the
 * module graph, scalar fallback otherwise); the contract holds either way.
 */
class DevSurfaceV7SimdContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV7Simd(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV7Simd(mol, solventRadius, tessLevel);
    }
}
