package cz.cuni.cusbg.surface;

/**
 * Copy-free neighbor access for sources that already store every atom's neighbors contiguously in a
 * single CSR array (start offsets + a packed index array). The default {@link NeighborSource} contract
 * appends atom {@code i}'s neighbors into a caller-supplied {@code IntArrayList} - which, for a CSR
 * source, means copying its slice out element by element on every query (profiling showed that copy,
 * {@code IntArrayList.add}, at ~15% of single-thread CPU). A source implementing this interface instead
 * hands the engine its backing array and the slice bounds, so {@link DevSurfaceV1Soa}'s neighbor
 * pre-pass reads {@code adjacency()[neighborStart(i) .. neighborEnd(i))} directly with no per-query
 * copy.
 *
 * <p>Opt-in: the engine uses this path only when explicitly enabled (its {@code directNeighbors} flag),
 * so a variant can keep the copy path for an apples-to-apples baseline even when its source also
 * implements this interface. Output is identical either way (same neighbor set, same order).
 */
interface DirectNeighborSource extends NeighborSource {

    /** Backing array holding all atoms' neighbor indices in CSR order. Not copied; do not mutate. */
    int[] adjacency();

    /** Start offset (inclusive) of atom {@code i}'s neighbors in {@link #adjacency()}. */
    int neighborStart(int i);

    /** End offset (exclusive) of atom {@code i}'s neighbors in {@link #adjacency()}. */
    int neighborEnd(int i);
}
