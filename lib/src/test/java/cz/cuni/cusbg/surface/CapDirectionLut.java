package cz.cuni.cusbg.surface;

/**
 * Prototype cap→direction lookup table for the Lead #1 bitmask feasibility study (backlog §1b).
 * Given the fixed distinct-direction set of one tessellation level, it precomputes, for each
 * (cube-face grid cell of a cap axis) × (cosθ radius bin), a <em>superset</em> bitmask of the directions
 * that cap could bury. At runtime a neighbor's cap {@code (axis=diff/|diff|, cosθ=thresh/|diff|)} is
 * mapped to a cell + bin and the precomputed mask narrows which directions need an exact dot-test.
 *
 * <p>Soundness (every truly-buried direction is in the returned mask) is what keeps the eventual scan
 * bit-exact; it is verified empirically by {@link BitmaskFeasibilityScan}'s violation counter. This
 * class is a measurement prototype (test scope) — it is deliberately allowed to use {@code sqrt}/{@code acos}
 * so the selectivity it measures is accurate; the sqrt-free runtime formulation is a separate concern
 * once the selectivity is shown to be worth pursuing.
 */
final class CapDirectionLut {

    final int numDir;
    final int words;          // longs per direction bitmask = ceil(numDir/64)
    private final int grid;   // G: cells per face axis = G*G per face, 6 faces
    private final int bins;   // B: cosθ bins over [-1, 1]
    private final long[] masks;   // [cell*bins + bin] flattened, each `words` longs
    private final long[] fullMask;

    /** Build the LUT for the given direction set at resolution {@code grid}×{@code grid} per face, {@code bins} cosθ bins. */
    CapDirectionLut(double[] dx, double[] dy, double[] dz, int numDir, int grid, int bins) {
        this.numDir = numDir;
        this.words = (numDir + 63) >>> 6;
        this.grid = grid;
        this.bins = bins;
        int numCells = 6 * grid * grid;
        this.masks = new long[numCells * bins * words];
        this.fullMask = new long[words];
        for (int d = 0; d < numDir; d++) setBit(fullMask, 0, d);

        // For each cell: geometric center (normalized) + angular radius (max over a fine sub-sampling).
        final int sub = 7;   // sub-samples per axis within a cell patch
        for (int face = 0; face < 6; face++) {
            int axis = face >> 1;            // 0=x,1=y,2=z
            double sign = (face & 1) == 0 ? 1.0 : -1.0;
            for (int iu = 0; iu < grid; iu++) {
                for (int iw = 0; iw < grid; iw++) {
                    double u0 = -1 + 2.0 * iu / grid, u1 = -1 + 2.0 * (iu + 1) / grid;
                    double w0 = -1 + 2.0 * iw / grid, w1 = -1 + 2.0 * (iw + 1) / grid;
                    double cx = 0, cy = 0, cz = 0;
                    // pass 1: normalized mean direction of the patch
                    for (int su = 0; su <= sub; su++) {
                        double u = u0 + (u1 - u0) * su / sub;
                        for (int sw = 0; sw <= sub; sw++) {
                            double w = w0 + (w1 - w0) * sw / sub;
                            double[] v = mapToSphere(axis, sign, u, w);
                            cx += v[0]; cy += v[1]; cz += v[2];
                        }
                    }
                    double cn = Math.sqrt(cx * cx + cy * cy + cz * cz);
                    cx /= cn; cy /= cn; cz /= cn;
                    // pass 2: angular radius = max angle from center to any patch sub-sample
                    double minDot = 1.0;
                    for (int su = 0; su <= sub; su++) {
                        double u = u0 + (u1 - u0) * su / sub;
                        for (int sw = 0; sw <= sub; sw++) {
                            double w = w0 + (w1 - w0) * sw / sub;
                            double[] v = mapToSphere(axis, sign, u, w);
                            double dot = cx * v[0] + cy * v[1] + cz * v[2];
                            if (dot < minDot) minDot = dot;
                        }
                    }
                    // margin for finite sub-sampling: half a sub-step's worth of extra angle
                    double cellRadius = Math.acos(clamp(minDot)) + 2.5 * (Math.PI / 2) / (grid * sub);

                    int cell = (face * grid + iu) * grid + iw;
                    for (int b = 0; b < bins; b++) {
                        // bin b covers cosθ in [edge(b), edge(b+1)); most inclusive cap = smallest cosθ = edge(b)
                        double cosLeft = -1 + 2.0 * b / bins;
                        double thetaMax = Math.acos(clamp(cosLeft));
                        double reach = thetaMax + cellRadius;
                        int base = (cell * bins + b) * words;
                        if (reach >= Math.PI) {
                            System.arraycopy(fullMask, 0, masks, base, words);
                            continue;
                        }
                        double cosReach = Math.cos(reach);
                        for (int d = 0; d < numDir; d++) {
                            if (cx * dx[d] + cy * dy[d] + cz * dz[d] > cosReach) setBit(masks, base, d);
                        }
                    }
                }
            }
        }
    }

    /**
     * Map face-local (u,w) to a unit sphere vector for face {@code axis}/{@code sign}. The off-axis
     * components carry the dominant-axis sign so that {@code maskFor}'s runtime {@code u=offaxis/a}
     * (with {@code a} the SIGNED dominant component) reproduces exactly this {@code (u,w)} — otherwise
     * the negative-dominant faces would be sign-flipped and the mask would miss directions.
     */
    private static double[] mapToSphere(int axis, double sign, double u, double w) {
        double x, y, z;
        switch (axis) {
            case 0 -> { x = sign; y = u * sign; z = w * sign; }
            case 1 -> { x = u * sign; y = sign; z = w * sign; }
            default -> { x = u * sign; y = w * sign; z = sign; }
        }
        double n = Math.sqrt(x * x + y * y + z * z);
        return new double[]{x / n, y / n, z / n};
    }

    /**
     * Return the precomputed superset mask for a neighbor's cap into {@code out} (length {@code words}).
     * {@code diffX/Y/Z} and {@code thresh} are the raw scan inputs (the buried test is {@code diff·p > thresh}).
     */
    void maskFor(double dxN, double dyN, double dzN, double thresh, long[] out) {
        double len2 = dxN * dxN + dyN * dyN + dzN * dzN;
        double len = Math.sqrt(len2);
        double cos = thresh / len;             // cosθ of the cap
        if (cos >= 1.0) { java.util.Arrays.fill(out, 0L); return; }   // empty cap, buries nothing
        if (cos <= -1.0) { System.arraycopy(fullMask, 0, out, 0, words); return; }

        // cube-face cell of the axis (parallel to diff, so magnitude irrelevant)
        double ax = Math.abs(dxN), ay = Math.abs(dyN), az = Math.abs(dzN);
        int axis; double a, bb, cc;
        if (ax >= ay && ax >= az) { axis = 0; a = dxN; bb = dyN; cc = dzN; }
        else if (ay >= az)        { axis = 1; a = dyN; bb = dxN; cc = dzN; }
        else                      { axis = 2; a = dzN; bb = dxN; cc = dyN; }
        int face = axis * 2 + (a < 0 ? 1 : 0);
        double u = bb / a, w = cc / a;          // both in [-1,1]
        int iu = clampIdx((int) ((u + 1) * 0.5 * grid));
        int iw = clampIdx((int) ((w + 1) * 0.5 * grid));
        int cell = (face * grid + iu) * grid + iw;
        int bin = (int) ((cos + 1) * 0.5 * bins);
        if (bin < 0) bin = 0; else if (bin >= bins) bin = bins - 1;
        int base = (cell * bins + bin) * words;
        System.arraycopy(masks, base, out, 0, words);
    }

    private int clampIdx(int i) { return i < 0 ? 0 : (i >= grid ? grid - 1 : i); }
    private static double clamp(double c) { return c < -1 ? -1 : (c > 1 ? 1 : c); }
    private static void setBit(long[] m, int base, int d) { m[base + (d >>> 6)] |= 1L << (d & 63); }
}
