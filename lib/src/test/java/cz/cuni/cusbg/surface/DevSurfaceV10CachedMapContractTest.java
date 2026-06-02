package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV10CachedMap}.
 * Confirms that caching the distinct-direction mapping process-wide (rather than rebuilding it per
 * surface) reproduces the reference result. Runs on whichever scan path the test JVM selected.
 */
class DevSurfaceV10CachedMapContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV10CachedMap(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV10CachedMap(mol, solventRadius, tessLevel);
    }
}
