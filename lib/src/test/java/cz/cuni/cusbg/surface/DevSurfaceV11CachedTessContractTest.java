package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for
 * {@link DevSurfaceV11CachedTess}. Confirms that reusing
 * a process-wide cached tessellation (rather than rebuilding it per surface) reproduces the reference
 * result. Runs on whichever scan path the test JVM selected.
 */
class DevSurfaceV11CachedTessContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV11CachedTess(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV11CachedTess(mol, solventRadius, tessLevel);
    }
}
