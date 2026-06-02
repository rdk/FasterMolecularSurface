package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;

/**
 * Flat uniform-cell grid neighbor index (optimization step 2): a {@link CellGrid} queried with a
 * 3x3x3 cell scan per atom.
 *
 * <p>Reproduces {@link NeighborList}'s neighbor sets exactly: cell size = cutoff = {@code 2*radius},
 * a 3x3x3 cell scan, the {@code dist^2 < cutoff^2} test, and self-exclusion. (With cell size equal to
 * the cutoff, any atom within the cutoff differs by at most one cell index per axis, so the 3x3x3
 * neighbourhood is exhaustive.)
 *
 * <p>Where it differs from {@code NeighborList} is the data structure: instead of a
 * {@code HashMap<Key, IntArrayList>} keyed by boxed object keys (object hashing + a {@code new Key}
 * probe per query), it uses the {@link CellGrid} CSR layout over flat coordinate arrays. Lookups are
 * O(1) array indexing, there is no per-query allocation, and {@code IAtom.getPoint3d()} never appears
 * in the hot loop. Each query recomputes the candidates' distances; {@link SymmetricCellGridNeighborList}
 * trades that for a one-shot symmetric precompute when the per-atom query cost dominates.
 */
final class CellGridNeighborList implements NeighborSource {

    private final CellGrid g;

    CellGridNeighborList(double[] ax, double[] ay, double[] az, double radius) {
        this.g = new CellGrid(ax, ay, az, radius);
    }

    @Override
    public void getNeighborsInto(int i, IntArrayList out) {
        double xi = g.ax[i], yi = g.ay[i], zi = g.az[i];
        int bx = g.cx[i], by = g.cy[i], bz = g.cz[i];
        for (int dx = -1; dx <= 1; dx++) {
            int ix = bx + dx; if (ix < g.minX || ix >= g.minX + g.nx) continue;
            for (int dy = -1; dy <= 1; dy++) {
                int iy = by + dy; if (iy < g.minY || iy >= g.minY + g.ny) continue;
                for (int dz = -1; dz <= 1; dz++) {
                    int iz = bz + dz; if (iz < g.minZ || iz >= g.minZ + g.nz) continue;
                    int c = g.cellIndex(ix, iy, iz);
                    for (int p = g.cellStart[c], end = g.cellStart[c + 1]; p < end; p++) {
                        int j = g.cellAtoms[p];
                        if (j == i) continue;
                        double x12 = g.ax[j] - xi, y12 = g.ay[j] - yi, z12 = g.az[j] - zi;
                        if (x12 * x12 + y12 * y12 + z12 * z12 < g.cutoff2) out.add(j);
                    }
                }
            }
        }
    }
}
