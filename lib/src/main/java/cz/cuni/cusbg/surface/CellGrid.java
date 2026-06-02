package cz.cuni.cusbg.surface;

/**
 * Flat uniform-cell spatial grid over a set of atom coordinates, shared by the cell-grid neighbor
 * sources ({@link CellGridNeighborList}, {@link SymmetricCellGridNeighborList}).
 *
 * <p>Cell size equals the neighbor cutoff ({@code 2*radius}), so any two atoms within the cutoff
 * differ by at most one cell index per axis and a 3x3x3 cell scan around an atom is exhaustive.
 * Atoms are bucketed into a CSR layout ({@link #cellStart} prefix sums + {@link #cellAtoms} grouped
 * indices) built with a counting sort, so cell contents are O(1) array slices with no per-query
 * allocation.
 *
 * <p><b>Precondition:</b> a <em>dense</em> grid spanning the input bounding box; suited to spatially
 * compact inputs (a molecular structure packs many atoms per cell, so {@code cells << atoms}). For
 * sparse or outlier coordinates the cell count explodes and the constructor fails loudly rather than
 * allocating a pathologically large array; the sparse hash-box {@link NeighborList} handles those.
 *
 * <p>Fields are package-private and read directly by the neighbor sources' hot loops.
 */
final class CellGrid {

    final double[] ax, ay, az;
    final double cutoff2;
    final double boxSize;
    final int minX, minY, minZ, nx, ny, nz;
    final int[] cellStart;   // prefix sums, length numCells+1
    final int[] cellAtoms;   // atom indices grouped by cell, length n
    final int[] cx, cy, cz;  // per-atom integer cell coordinates

    CellGrid(double[] ax, double[] ay, double[] az, double radius) {
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

    int cellIndex(int ix, int iy, int iz) {
        return (ix - minX) + nx * ((iy - minY) + ny * (iz - minZ));
    }
}
