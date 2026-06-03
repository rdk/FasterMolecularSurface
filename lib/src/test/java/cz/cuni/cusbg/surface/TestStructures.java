package cz.cuni.cusbg.surface;

import org.junit.jupiter.params.provider.Arguments;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import java.io.InputStream;
import java.util.List;
import java.util.stream.Stream;

/**
 * Catalogue of representative PDB test structures and helpers to load them.
 *
 * <p>The ten structures span a wide size range (a few hundred to a few thousand atoms)
 * so the contract / equivalence / golden test batteries exercise small, medium and large
 * inputs. They are the shared corpus for every {@link MolecularSurface} variant.
 */
final class TestStructures {

    private TestStructures() {}

    /**
     * Representative structures, ordered by increasing size. Atom counts are the number of
     * atoms in the first {@link IAtomContainer} returned by CDK's PDB reader (what the
     * surface is actually computed over).
     */
    enum Structure {
        CRAMBIN("1CRN", 327),          // crambin, smallest
        BPTI("4PTI", 514),             // bovine pancreatic trypsin inhibitor
        VILLIN("1VII", 596),           // villin headpiece
        UBIQUITIN("1UBQ", 660),        // ubiquitin
        LYSOZYME_HEW("1AKI", 1079),    // hen egg-white lysozyme
        LYSOZYME("2LYZ", 1102),        // lysozyme
        MYOGLOBIN("1MBN", 1260),       // sperm whale myoglobin
        COBALAMIN("3CI3", 1842),       // methylmalonyl-CoA mutase (exercises Co VdW fallback)
        HEMOGLOBIN_HHO("1HHO", 2396),  // oxyhaemoglobin alpha chain
        HEMOGLOBIN_4HHB("4HHB", 4779); // deoxyhaemoglobin, largest

        final String pdbId;
        final int atomCount;

        Structure(String pdbId, int atomCount) {
            this.pdbId = pdbId;
            this.atomCount = atomCount;
        }

        String resourcePath() {
            return "/data/pdb/" + pdbId + ".pdb";
        }

        IAtomContainer load() {
            return TestStructures.load(resourcePath());
        }
    }

    /** Load the first atom container from a PDB resource on the classpath. */
    static IAtomContainer load(String resourcePath) {
        IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
        IChemFile chemFile;
        try (InputStream in = TestStructures.class.getResourceAsStream(resourcePath);
             PDBReader pdbr = new PDBReader(in)) {
            if (in == null) {
                throw new IllegalStateException("Missing test resource: " + resourcePath);
            }
            chemFile = pdbr.read(bldr.newInstance(IChemFile.class));
        } catch (Exception e) {
            throw new IllegalStateException("Failed to read PDB resource " + resourcePath, e);
        }
        List<IAtomContainer> containers = ChemFileManipulator.getAllAtomContainers(chemFile);
        if (containers.isEmpty()) {
            throw new IllegalStateException("No atom containers in " + resourcePath);
        }
        return containers.get(0);
    }

    /** JUnit {@code @MethodSource} provider streaming every structure as a test argument. */
    static Stream<Arguments> all() {
        return Stream.of(Structure.values()).map(Arguments::of);
    }

    /** Load every structure's first atom container (shared corpus for the benchmark/quality harnesses). */
    static IAtomContainer[] loadAll() {
        Structure[] all = Structure.values();
        IAtomContainer[] out = new IAtomContainer[all.length];
        for (int i = 0; i < all.length; i++) out[i] = all[i].load();
        return out;
    }

    /**
     * Load every structure except the ones CDK's reference {@code NumericalSurface} cannot process
     * (3CI3 / cobalt NPEs there). Used as the common benchmark corpus so the CDK baseline and every
     * library variant are timed on the identical structure set, keeping speedup ratios clean.
     */
    static IAtomContainer[] loadAllCdkSafe() {
        List<IAtomContainer> out = new java.util.ArrayList<>();
        for (Structure s : Structure.values()) {
            if (s == Structure.COBALAMIN) continue;
            out.add(s.load());
        }
        return out.toArray(new IAtomContainer[0]);
    }

    /**
     * Parameter combinations for cross-variant equivalence: each {@code (structure, solventRadius,
     * tessLevel)}. Covers the tessellation sweep at the standard 1.4 A solvent radius — including
     * {@code tess=2}, p2rank's operating point — plus the Van der Waals path ({@code solvent=0}).
     */
    static Stream<Arguments> structureConfigs() {
        double[][] configs = {
                {1.4, 2}, {1.4, 3}, {1.4, 4},   // tessellation sweep at the standard solvent radius
                {0.0, 4}                        // Van der Waals surface path
        };
        List<Arguments> out = new java.util.ArrayList<>();
        for (Structure s : Structure.values()) {
            for (double[] c : configs) {
                out.add(Arguments.of(s, c[0], (int) c[1]));
            }
        }
        return out.stream();
    }
}
