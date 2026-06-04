package cz.cuni.cusbg.surface;

/**
 * Finer-celled variant of {@link CellGrid}: the cell side is the neighbor cutoff divided by
 * {@code cellsPerCutoff} (vs {@code = cutoff} in {@link CellGrid}). A neighbor build that searches a
 * {@code ±cellsPerCutoff} stencil over these smaller cells covers the same cutoff while sweeping a
 * smaller candidate volume — e.g. a 5×5×5 stencil over half-cutoff cells encloses
 * {@code 125·(c/2)³ = 15.6·c³} vs {@link CellGrid}'s 3×3×3 {@code 27·c³}, ~42% less — so the distance
 * pass visits fewer non-neighbor candidates. Used by {@link TightCellGridNeighborList}; the frozen
 * {@link CellGrid} is left byte-identical so existing variants are unaffected.
 *
 * <p>{@code cellsPerCutoff == 1} reproduces {@link CellGrid}. Same dense-grid precondition and CSR layout.
 */
final class TightCellGrid {

    final double[] ax, ay, az;
    final double boxSize;
    final int minX, minY, minZ, nx, ny, nz;
    final int[] cellStart;   // prefix sums, length numCells+1
    final int[] cellAtoms;   // atom indices grouped by cell, length n
    final int[] cx, cy, cz;  // per-atom integer cell coordinates

    /**
     * @param n              atoms to bin (entries {@code [0, n)}; arrays may be larger, e.g. arena-reused)
     * @param radius         max expanded radius; the cutoff is {@code 2*radius} (as in {@link CellGrid})
     * @param cellsPerCutoff cells per cutoff length; {@code 1} = {@link CellGrid}, {@code 2} = half cells
     */
    TightCellGrid(double[] ax, double[] ay, double[] az, int n, double radius, int cellsPerCutoff) {
        this.ax = ax; this.ay = ay; this.az = az;
        this.boxSize = (2 * radius) / cellsPerCutoff;   // CellGrid uses 2*radius (= cutoff); finer by cellsPerCutoff

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
        // Same dense-grid guard as CellGrid, relaxed by cellsPerCutoff³ (a finer grid legitimately has
        // ~k³ more cells): fail loudly on sparse/outlier coordinates instead of allocating pathologically.
        long cells = (long) nx * ny * nz;
        long k3 = (long) cellsPerCutoff * cellsPerCutoff * cellsPerCutoff;
        long maxCells = Math.max(1L << 20, 64L * k3 * n);
        if (cells > maxCells || cells > Integer.MAX_VALUE - 8) {
            throw new IllegalArgumentException(
                    "coordinate span too sparse for a dense cell grid: " + nx + "x" + ny + "x" + nz
                            + " = " + cells + " cells for " + n + " atoms (cellsPerCutoff=" + cellsPerCutoff
                            + "); use the hash-box NeighborList for sparse/outlier inputs");
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

    int cellIndex(int ix, int iy, int iz) {
        return (ix - minX) + nx * ((iy - minY) + ny * (iz - minZ));
    }
}
