#!/usr/bin/env bash
#
# Reproducible JMH benchmark run for the surface kernel.
#
# Run as your NORMAL user (not under sudo): the JMH run stays unprivileged; only the CPU-pinning step
# escalates via sudo for the two root-only sysfs writes (governor + turbo). Pinning makes sub-1.1x deltas
# stable on this turbo/downclock-bound box; without sudo the harness still runs, just noisier.
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

# Escalate ONLY for the two root-only sysfs writes; the rest of the script (Gradle, the JVM) stays
# unprivileged so it leaves no root-owned files behind.
pin() {
  if [ "$(id -u)" -eq 0 ]; then
    echo "[bench] note: running as root; prefer plain ./bench.sh so the JVM stays unprivileged"
  elif ! command -v sudo >/dev/null 2>&1 || ! sudo -v 2>/dev/null; then
    echo "[bench] no usable sudo - skipping governor/turbo pin (numbers will be noisier; sub-1.1x deltas unreliable)"
    return
  fi
  local sudo=""; [ "$(id -u)" -ne 0 ] && sudo="sudo"
  for g in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
    echo performance | $sudo tee "$g" >/dev/null 2>&1 || true
  done
  echo 1 | $sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo >/dev/null 2>&1 || true
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
