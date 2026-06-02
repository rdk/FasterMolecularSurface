package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV9Dedup}.
 * The exact point/area golden values confirm that scanning each distinct tessellation direction once and
 * re-expanding into the original point order reproduces the reference result. Runs on whichever scan
 * path the test JVM selected (dedup SIMD when {@code jdk.incubator.vector} is on the module graph, scalar
 * fallback otherwise); the contract holds either way.
 */
class DevSurfaceV9DedupContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV9Dedup(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV9Dedup(mol, solventRadius, tessLevel);
    }
}
