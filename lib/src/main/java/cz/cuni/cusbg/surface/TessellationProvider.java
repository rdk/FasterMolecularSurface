package cz.cuni.cusbg.surface;

/**
 * Supplies the {@link Tessellation} for a level to {@link SoaNumericalSurface}. Passed to the engine
 * constructor (like the other strategies) so a variant can opt into a shared cached tessellation
 * without changing the default behavior of existing variants.
 *
 * <p>Both implementations return bit-identical data; they differ only in whether the tessellation is
 * rebuilt per surface ({@link #FRESH}) or fetched from a process-wide cache ({@link #CACHED}).
 */
@FunctionalInterface
interface TessellationProvider {

    Tessellation get(int tesslevel);

    /** Build the tessellation fresh for every surface (the engine default; matches the original code). */
    TessellationProvider FRESH = Tessellation::build;

    /** Reuse a process-wide tessellation shared across all builds at a level. */
    TessellationProvider CACHED = Tessellation::cached;
}
