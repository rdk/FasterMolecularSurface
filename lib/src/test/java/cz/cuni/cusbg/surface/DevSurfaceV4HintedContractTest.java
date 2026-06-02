package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV4Hinted}
 * (SoA + flat cell grid + last-occluder-first scan hint, no sort). The exact point/area golden
 * values confirm that the scan-order hint does not change the result.
 */
class DevSurfaceV4HintedContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV4Hinted(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV4Hinted(mol, solventRadius, tessLevel);
    }
}
