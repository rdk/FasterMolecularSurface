package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Runs the full {@link MolecularSurface} contract battery (incl. the golden baseline) against
 * {@link DevSurfaceV1Soa}. Passing {@code matchesGoldenValues} proves the SoA variant
 * reproduces the same point counts / exposed-atom counts / areas as {@link FasterNumericalSurface}.
 */
class DevSurfaceV1SoaContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV1Soa(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV1Soa(mol, solventRadius, tessLevel);
    }
}
