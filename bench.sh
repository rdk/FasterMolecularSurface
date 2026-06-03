#!/usr/bin/env bash
#
# Reproducible JMH benchmark run for the surface kernel.
#
# Pins the CPU frequency governor + disables turbo (needs root) so sub-1.1x deltas are stable on this
# turbo/downclock-bound box, stamps the environment next to the results, then runs the JMH harness to
# CSV. Degrades gracefully without root (warns, records the un-pinned state, still runs).
#
# Scope/params live in lib/build.gradle's `jmh { }` block and SurfaceBench's @Param defaults (widen the
# variantId list for the full V1..V19 ladder). Profiling: uncomment `profilers` in that block.
#
# After the run, render the canonical ladder with:  python3 bench-table.py
#
set -euo pipefail
cd "$(dirname "$0")"

RESULTS_DIR="lib/build/results/jmh"
mkdir -p "$RESULTS_DIR"
SHA="$(git rev-parse --short HEAD 2>/dev/null || echo nogit)"
STAMP="$RESULTS_DIR/env-$SHA.txt"

pin() {
  if [ "$(id -u)" -ne 0 ]; then
    echo "[bench] not root - skipping governor/turbo pin (numbers will be noisier; sub-1.1x deltas unreliable)"
    return
  fi
  for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do echo performance > "$g" 2>/dev/null || true; done
  echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || true
  echo "[bench] governor=performance, turbo disabled"
}
pin

{
  echo "git:       $(git rev-parse HEAD 2>/dev/null || echo n/a)"
  echo "date_utc:  $(date -u +%FT%TZ)"
  echo "cpu:       $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2- | sed 's/^ //' || echo n/a)"
  echo "cores:     $(nproc 2>/dev/null || echo n/a)"
  echo "governor:  $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo n/a)"
  echo "no_turbo:  $(cat /sys/devices/system/cpu/intel_pstate/no_turbo 2>/dev/null || echo n/a)"
  echo "java:      $(java -version 2>&1 | head -1)"
} | tee "$STAMP"

echo "[bench] running ./gradlew jmh (full @Param matrix; narrow it in lib/build.gradle to go faster)"
./gradlew jmh "$@"

echo "[bench] results: $RESULTS_DIR/results.csv    env stamp: $STAMP"
echo "[bench] render the ladder with:  python3 bench-table.py"
