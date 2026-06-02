package cz.cuni.cusbg.surface;

import javax.vecmath.Point3d;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * Shared assertion: two {@link MolecularSurface} implementations produce bit-for-bit identical output
 * (per-atom areas, total surface-point count, and every surface-point coordinate). Used by the
 * variant equivalence tests to prove a data-layout / neighbor-index / neighbor-ordering change does
 * not alter the result at a given (solvent radius, tessellation level).
 */
final class VariantEquivalence {

    private VariantEquivalence() {}

    static void assertBitForBit(TestStructures.Structure s, double solvent, int tess,
                                MolecularSurface ref, MolecularSurface variant) {
        String cfg = s + " (solvent=" + solvent + ", tess=" + tess + ")";

        double[] ra = ref.getAllSurfaceAreas();
        double[] va = variant.getAllSurfaceAreas();
        assertEquals(ra.length, va.length, () -> cfg + " area array length");
        for (int i = 0; i < ra.length; i++) {
            final int idx = i;
            assertEquals(ra[i], va[i], 0.0, () -> cfg + " per-atom area[" + idx + "]");
        }

        Point3d[] rp = ref.getAllSurfacePoints();
        Point3d[] vp = variant.getAllSurfacePoints();
        assertEquals(rp.length, vp.length, () -> cfg + " total surface-point count");
        for (int i = 0; i < rp.length; i++) {
            final int idx = i;
            assertEquals(rp[i].x, vp[i].x, 0.0, () -> cfg + " point[" + idx + "].x");
            assertEquals(rp[i].y, vp[i].y, 0.0, () -> cfg + " point[" + idx + "].y");
            assertEquals(rp[i].z, vp[i].z, 0.0, () -> cfg + " point[" + idx + "].z");
        }
    }
}
