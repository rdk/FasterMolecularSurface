package cz.cuni.cusbg.surface;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

/**
 * Loader for the pinned golden surface values.
 *
 * <p>The baseline lives in {@code /golden/surface-golden.csv} and records, per structure at the
 * default parameters (solvent radius 1.4 A, tessellation level 4): atom count, number of exposed
 * (surface) atoms, total number of surface points, and total surface area. These were captured
 * from {@link FasterNumericalSurface} and validated against CDK's reference by
 * {@link CdkEquivalenceTest}. Regenerate with {@link GoldenValuesGenerator}.
 */
final class GoldenValues {

    static final class Row {
        final String pdbId;
        final int atomCount;
        final int surfaceAtomCount;
        final int totalPoints;
        final double totalArea;

        Row(String pdbId, int atomCount, int surfaceAtomCount, int totalPoints, double totalArea) {
            this.pdbId = pdbId;
            this.atomCount = atomCount;
            this.surfaceAtomCount = surfaceAtomCount;
            this.totalPoints = totalPoints;
            this.totalArea = totalArea;
        }
    }

    static final String RESOURCE = "/golden/surface-golden.csv";

    private static final Map<String, Row> ROWS = load();

    private static Map<String, Row> load() {
        Map<String, Row> map = new HashMap<>();
        try (InputStream in = GoldenValues.class.getResourceAsStream(RESOURCE)) {
            if (in == null) {
                throw new IllegalStateException("Missing golden resource " + RESOURCE
                        + " - regenerate with GoldenValuesGenerator");
            }
            try (BufferedReader r = new BufferedReader(new InputStreamReader(in, StandardCharsets.UTF_8))) {
                String line;
                while ((line = r.readLine()) != null) {
                    String t = line.trim();
                    if (t.isEmpty() || t.startsWith("#")) continue;
                    String[] p = t.split(",");
                    if (p.length < 5) continue;
                    map.put(p[0], new Row(p[0],
                            Integer.parseInt(p[1].trim()),
                            Integer.parseInt(p[2].trim()),
                            Integer.parseInt(p[3].trim()),
                            Double.parseDouble(p[4].trim())));
                }
            }
        } catch (Exception e) {
            throw new IllegalStateException("Failed to load " + RESOURCE, e);
        }
        return map;
    }

    static Row get(String pdbId) {
        Row row = ROWS.get(pdbId);
        if (row == null) {
            throw new IllegalStateException("No golden row for " + pdbId
                    + " - regenerate with GoldenValuesGenerator");
        }
        return row;
    }

    private GoldenValues() {}
}
