package cz.cuni.cusbg.surface;

import com.carrotsearch.hppc.IntArrayList;

/**
 * Minimal neighbor-query abstraction shared by the surface variants: append to {@code out} the
 * indices of all atoms within the construction cutoff of atom {@code i} (excluding {@code i}).
 * Lets {@link SoaNumericalSurface} stay agnostic about whether neighbors come from the
 * hash-box {@link NeighborList} or the flat {@link CellGridNeighborList}.
 */
interface NeighborSource {
    void getNeighborsInto(int i, IntArrayList out);
}
