package cz.cuni.cusbg.surface;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import javax.vecmath.Point3d;
import java.util.List;
import java.util.Map;

/**
 * Adapts CDK's reference {@link NumericalSurface} to {@link MolecularSurface} so the benchmark/quality
 * registry can treat the baseline uniformly with the library variants. CDK already exposes every
 * {@code MolecularSurface} method; it just does not declare the interface and needs an explicit
 * {@code calculateSurface()} call, which this adapter performs in the constructor.
 */
final class CdkSurfaceAdapter implements MolecularSurface {

    private final NumericalSurface ns;

    CdkSurfaceAdapter(IAtomContainer mol, double solventRadius, int tessLevel) {
        this.ns = new NumericalSurface(mol, solventRadius, tessLevel);
        ns.calculateSurface();
    }

    @Override public Point3d[] getAllSurfacePoints() { return ns.getAllSurfacePoints(); }
    @Override public Map<IAtom, List<Point3d>> getAtomSurfaceMap() { return ns.getAtomSurfaceMap(); }
    @Override public Point3d[] getSurfacePoints(int atomIdx) throws CDKException { return ns.getSurfacePoints(atomIdx); }
    @Override public double getSurfaceArea(int atomIdx) throws CDKException { return ns.getSurfaceArea(atomIdx); }
    @Override public double[] getAllSurfaceAreas() { return ns.getAllSurfaceAreas(); }
    @Override public double getTotalSurfaceArea() { return ns.getTotalSurfaceArea(); }
}
