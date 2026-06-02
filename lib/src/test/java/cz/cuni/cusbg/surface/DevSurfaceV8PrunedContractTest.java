package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV8Pruned}.
 * The exact point/area golden values confirm that the per-pair occlusion cutoff (dropping neighbors
 * that can never bury) and the 256-bit SIMD scan reproduce the reference result. Runs on whichever
 * scan path the test JVM selected (vectorized when {@code jdk.incubator.vector} is on the module graph,
 * scalar fallback otherwise); the contract holds either way.
 */
class DevSurfaceV8PrunedContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV8Pruned(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV8Pruned(mol, solventRadius, tessLevel);
    }
}
