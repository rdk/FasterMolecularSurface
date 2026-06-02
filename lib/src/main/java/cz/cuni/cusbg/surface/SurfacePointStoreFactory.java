package cz.cuni.cusbg.surface;

/**
 * Builds the per-surface {@link SurfacePointStore}. Given both the atom count (sizing per-atom
 * structures) and an estimate of the total surviving surface points, so a flat store can pre-size its
 * coordinate buffer in one allocation instead of growing it by doubling (the doubling copies were the
 * footprint regression that profiling found). The estimate is a hint only: a store may ignore it
 * (the list store does) and any store must still grow correctly if it is exceeded.
 */
@FunctionalInterface
interface SurfacePointStoreFactory {

    /**
     * @param numAtoms        number of atoms (sizes per-atom structures, e.g. CSR row pointers)
     * @param estimatedPoints rough upper estimate of total surface points across all atoms (capacity
     *                        hint; not a hard bound)
     */
    SurfacePointStore create(int numAtoms, int estimatedPoints);
}
