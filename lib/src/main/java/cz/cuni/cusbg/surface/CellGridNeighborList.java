package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;

/**
 * Flat uniform-cell grid neighbor index (optimization step 2).
 *
 * <p>Reproduces {@link NeighborList}'s neighbor sets exactly: cell size = cutoff = {@code 2*radius},
 * a 3x3x3 cell scan, the {@code dist^2 < cutoff^2} test, and self-exclusion. (With cell size equal to
 * the cutoff, any atom within the cutoff differs by at most one cell index per axis, so the 3x3x3
 * neighbourhood is exhaustive.)
 *
 * <p>Where it differs from {@code NeighborList} is the data structure: instead of a
 * {@code HashMap<Key, IntArrayList>} keyed by boxed object keys (object hashing + a {@code new Key}
 * probe per query), it uses a CSR cell layout ({@code cellStart[]} + {@code cellAtoms[]}, built with a
 * counting sort) over flat coordinate arrays. Lookups are O(1) array indexing, there is no per-query
 * allocation, and {@code IAtom.getPoint3d()} never appears in the hot loop.
 *
 * <p><b>Precondition:</b> this is a <em>dense</em> grid spanning the input's bounding box, so it suits
 * spatially compact inputs (a molecular structure packs many atoms per cell, so {@code cells << atoms}).
 * For sparse or outlier coordinates (e.g. a stray atom thousands of A away) the cell count explodes; the
 * constructor fails loudly rather than allocating a pathologically large array. The sparse hash-box
 * {@link NeighborList} has no such limit and is the right choice for those inputs.
 */
final class CellGridNeighborList implements NeighborSource {

    private final double[] ax, ay, az;
    private final double cutoff2;
    private final double boxSize;
    private final int minX, minY, minZ, nx, ny, nz;
    private final int[] cellStart;   // prefix sums, length numCells+1
    private final int[] cellAtoms;   // atom indices grouped by cell, length n
    private final int[] cx, cy, cz;  // per-atom integer cell coordinates

    CellGridNeighborList(double[] ax, double[] ay, double[] az, double radius) {
        this.ax = ax; this.ay = ay; this.az = az;
        this.boxSize = 2 * radius;
        this.cutoff2 = boxSize * boxSize;
        int n = ax.length;

        cx = new int[n]; cy = new int[n]; cz = new int[n];
        int mnx = Integer.MAX_VALUE, mny = Integer.MAX_VALUE, mnz = Integer.MAX_VALUE;
        int mxx = Integer.MIN_VALUE, mxy = Integer.MIN_VALUE, mxz = Integer.MIN_VALUE;
        for (int i = 0; i < n; i++) {
            int ix = (int) Math.floor(ax[i] / boxSize);
            int iy = (int) Math.floor(ay[i] / boxSize);
            int iz = (int) Math.floor(az[i] / boxSize);
            cx[i] = ix; cy[i] = iy; cz[i] = iz;
            if (ix < mnx) mnx = ix;  if (ix > mxx) mxx = ix;
            if (iy < mny) mny = iy;  if (iy > mxy) mxy = iy;
            if (iz < mnz) mnz = iz;  if (iz > mxz) mxz = iz;
        }
        if (n == 0) { mnx = mny = mnz = 0; mxx = mxy = mxz = 0; }
        minX = mnx; minY = mny; minZ = mnz;
        nx = mxx - mnx + 1; ny = mxy - mny + 1; nz = mxz - mnz + 1;
        // Guard the dense-grid size. A spatially compact molecular input packs many atoms per cell, so
        // cells << atoms; cells growing far beyond the atom count means sparse/outlier coordinates that
        // a dense grid cannot represent economically. Fail loudly far below an OOM allocation instead of
        // silently allocating hundreds of MB (or overflowing the int cell index). The sparse hash-box
        // NeighborList handles such inputs; this grid intentionally does not.
        long cells = (long) nx * ny * nz;
        long maxCells = Math.max(1L << 20, 64L * n);   // ~1M cell floor, else 64x the atom count
        if (cells > maxCells || cells > Integer.MAX_VALUE - 8) {
            throw new IllegalArgumentException(
                    "coordinate span too sparse for a dense cell grid: " + nx + "x" + ny + "x" + nz
                            + " = " + cells + " cells for " + n + " atoms; use the hash-box NeighborList for sparse/outlier inputs");
        }
        int numCells = (int) cells;

        int[] counts = new int[numCells];
        for (int i = 0; i < n; i++) counts[cellIndex(cx[i], cy[i], cz[i])]++;
        cellStart = new int[numCells + 1];
        for (int c = 0; c < numCells; c++) cellStart[c + 1] = cellStart[c] + counts[c];
        cellAtoms = new int[n];
        int[] cursor = cellStart.clone();
        for (int i = 0; i < n; i++) {
            cellAtoms[cursor[cellIndex(cx[i], cy[i], cz[i])]++] = i;
        }
    }

    private int cellIndex(int ix, int iy, int iz) {
        return (ix - minX) + nx * ((iy - minY) + ny * (iz - minZ));
    }

    @Override
    public void getNeighborsInto(int i, IntArrayList out) {
        double xi = ax[i], yi = ay[i], zi = az[i];
        int bx = cx[i], by = cy[i], bz = cz[i];
        for (int dx = -1; dx <= 1; dx++) {
            int ix = bx + dx; if (ix < minX || ix >= minX + nx) continue;
            for (int dy = -1; dy <= 1; dy++) {
                int iy = by + dy; if (iy < minY || iy >= minY + ny) continue;
                for (int dz = -1; dz <= 1; dz++) {
                    int iz = bz + dz; if (iz < minZ || iz >= minZ + nz) continue;
                    int c = cellIndex(ix, iy, iz);
                    for (int p = cellStart[c], end = cellStart[c + 1]; p < end; p++) {
                        int j = cellAtoms[p];
                        if (j == i) continue;
                        double x12 = ax[j] - xi, y12 = ay[j] - yi, z12 = az[j] - zi;
                        if (x12 * x12 + y12 * y12 + z12 * z12 < cutoff2) out.add(j);
                    }
                }
            }
        }
    }
}
