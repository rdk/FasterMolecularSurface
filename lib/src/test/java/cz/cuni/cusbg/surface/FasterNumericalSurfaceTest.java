package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemFile;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.PDBReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import javax.vecmath.Point3d;
import java.io.InputStream;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 *
 */
class FasterNumericalSurfaceTest {

    @Test
    public void testCranbinSurface() throws Exception {
        IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
        IChemFile chemFile;
        String path = "/data/pdb/1CRN.pdb";
        try (InputStream in = getClass().getResourceAsStream(path);
             PDBReader pdbr = new PDBReader(in)) {
            chemFile = pdbr.read(bldr.newInstance(IChemFile.class));
        }
        IAtomContainer mol     = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);
        FasterNumericalSurface surface = new FasterNumericalSurface(mol);
        Map<IAtom, List<Point3d>> map = surface.getAtomSurfaceMap();

        assertEquals(222, map.size());
        assertEquals(327, mol.getAtomCount());
    }
    
}