#!/usr/bin/env python3
"""Render the canonical speedup ladder (markdown) from the JMH CSV produced by ./bench.sh (or ./gradlew jmh).

Usage:  python3 bench-table.py [path/to/results.csv]

Speedup is CDK_score / variant_score (AverageTime mode: lower ms/op is faster). One table per `consume`
mode (AREA vs POINTS). This replaces the legacy hand-maintained ladder table with a single, one-harness,
CI-bearing source of truth - paste it into docs/performance-lessons.md as a dated JMH table (keep the
historical median-of-3 table; do not overwrite it).
"""
import csv
import collections
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "lib/build/results/jmh/results.csv"

try:
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
except FileNotFoundError:
    sys.exit(f"No JMH results at {path}. Run ./bench.sh (or ./gradlew jmh) first.")

if not rows:
    sys.exit(f"{path} is empty.")


def col(row, *names):
    for n in names:
        if n in row:
            return row[n]
    raise KeyError(f"none of {names} in CSV columns {list(row)}")


# (consume) -> variantId -> tess -> (score, error)
data = collections.defaultdict(lambda: collections.defaultdict(dict))
for r in rows:
    consume = col(r, "Param: consume", '"Param: consume"')
    tess = col(r, "Param: tess", '"Param: tess"')
    var = col(r, "Param: variantId", '"Param: variantId"')
    score = float(col(r, "Score", '"Score"'))
    raw_err = col(r, "Score Error (99.9%)", '"Score Error (99.9%)"')
    err = float(raw_err) if raw_err not in ("", "NaN") else float("nan")
    data[consume][var][tess] = (score, err)

tess_levels = sorted({t for c in data.values() for v in c.values() for t in v}, key=int)

for consume in sorted(data):
    print(f"\n### consume = {consume}  (avg ms/build over the CDK-safe corpus; x CDK = speedup)\n")
    header = "| variant | " + " | ".join(f"tess {t} (ms +/- CI)" for t in tess_levels) \
             + " | " + " | ".join(f"xCDK t{t}" for t in tess_levels) + " |"
    print(header)
    print("|" + "---|" * (1 + 2 * len(tess_levels)))
    cdk = data[consume].get("CDK", {})
    for var in data[consume]:
        cells, ratios = [], []
        for t in tess_levels:
            sc = data[consume][var].get(t)
            if sc is None:
                cells.append("-")
                ratios.append("-")
                continue
            cells.append(f"{sc[0]:.2f} +/- {sc[1]:.2f}" if sc[1] == sc[1] else f"{sc[0]:.2f}")
            c = cdk.get(t)
            ratios.append(f"{c[0] / sc[0]:.2f}x" if c else "-")
        print(f"| {var} | " + " | ".join(cells) + " | " + " | ".join(ratios) + " |")
print()
