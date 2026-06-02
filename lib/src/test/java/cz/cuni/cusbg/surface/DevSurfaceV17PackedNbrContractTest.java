package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * Full contract battery (incl. golden baseline) for {@link DevSurfaceV17PackedNbr}. Confirms that
 * packing the kept edges into a single {@code int[]} (rather than two {@code IntArrayList}s) before
 * building the CSR reproduces the reference result. Runs on whichever scan path the test JVM selected.
 */
class DevSurfaceV17PackedNbrContractTest extends AbstractMolecularSurfaceContractTest {

    @Override
    protected MolecularSurface create(IAtomContainer mol) {
        return new DevSurfaceV17PackedNbr(mol);
    }

    @Override
    protected MolecularSurface create(IAtomContainer mol, double solventRadius, int tessLevel) {
        return new DevSurfaceV17PackedNbr(mol, solventRadius, tessLevel);
    }
}
