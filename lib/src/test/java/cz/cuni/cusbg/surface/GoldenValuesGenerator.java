package cz.cuni.cusbg.surface;

import org.junit.jupiter.api.Test;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assumptions.assumeTrue;

/**
 * Regenerates the golden baseline ({@code /golden/surface-golden.csv}) from
 * {@link FasterNumericalSurface}'s current output.
 *
 * <p>Skipped by default. Run explicitly when the algorithm intentionally changes:
 * <pre>./gradlew test --tests '*GoldenValuesGenerator' -Dgolden.regenerate=true</pre>
 * then re-run the suite and review the diff to {@code surface-golden.csv} before committing.
 */
class GoldenValuesGenerator {

    @Test
    void regenerateGoldenCsv() throws IOException {
        assumeTrue(Boolean.getBoolean("golden.regenerate"),
                "Set -Dgolden.regenerate=true to regenerate the golden baseline");

        List<String> lines = new ArrayList<>();
        lines.add("# pdbId,atomCount,surfaceAtomCount,totalPoints,totalArea");
        lines.add("# golden baseline for FasterNumericalSurface at solventRadius=1.4, tessLevel=4");
        lines.add("# regenerate: ./gradlew test --tests '*GoldenValuesGenerator' -Dgolden.regenerate=true");
        for (TestStructures.Structure s : TestStructures.Structure.values()) {
            IAtomContainer mol = s.load();
            FasterNumericalSurface surf = new FasterNumericalSurface(mol);
            int atomCount = surf.getAllSurfaceAreas().length;
            int surfaceAtomCount = surf.getAtomSurfaceMap().size();
            int totalPoints = surf.getAllSurfacePoints().length;
            double totalArea = surf.getTotalSurfaceArea();
            lines.add(String.format("%s,%d,%d,%d,%.10g",
                    s.pdbId, atomCount, surfaceAtomCount, totalPoints, totalArea));
        }

        // Write into the source resources. Prefer the absolute dir forwarded by the build
        // (-Dgolden.resourcesDir); fall back to the module-relative path when run via Gradle from lib/.
        String resourcesDir = System.getProperty("golden.resourcesDir");
        Path out = (resourcesDir != null ? Paths.get(resourcesDir, "golden", "surface-golden.csv")
                                          : Paths.get("src/test/resources/golden/surface-golden.csv"));
        Files.createDirectories(out.getParent());
        Files.write(out, lines, StandardCharsets.UTF_8);
        System.out.println("Wrote golden baseline to " + out.toAbsolutePath());
    }
}
