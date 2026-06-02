package cz.cuni.cusbg.surface;

import org.openscience.cdk.interfaces.IAtom;

import java.util.concurrent.ConcurrentHashMap;

/**
 * Process-wide cache of van der Waals radii by element symbol (optimization D).
 *
 * <p>{@link FasterNumericalSurface#getVdwRadius(IAtom)} resolves the radius per atom via a
 * {@code toLowerCase} + periodic-table lookup (and boxes a {@code Double}). The result depends only on
 * the element symbol, so this memoizes it: each distinct symbol is resolved once for the whole process.
 * Returns the exact same value as the underlying lookup, so it is bit-for-bit safe.
 */
final class VdwRadiusCache {

    private static final ConcurrentHashMap<String, Double> CACHE = new ConcurrentHashMap<>();

    private VdwRadiusCache() {}

    static double get(IAtom atom) {
        String symbol = atom.getSymbol();
        Double r = CACHE.get(symbol);
        if (r != null) return r;
        double v = FasterNumericalSurface.getVdwRadius(atom);
        CACHE.putIfAbsent(symbol, v);
        return v;
    }
}
