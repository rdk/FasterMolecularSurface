#!/usr/bin/env bash
#
# Run the COMPLETE unit-test suite, including the experimental DevSurfaceV* ladder contract tests and
# the cross-variant equivalence harness that `./gradlew test` (and ./unit-test.sh) exclude by default.
# Runs in parallel (-PtestForks=auto -> CPU cores / 2); each fork can use up to the test heap, so this
# trades memory for speed. Use this before a release or after touching the engine / a DevSurface rung.
#
# Extra arguments are forwarded to Gradle, e.g.:
#   ./unit-test-all.sh -PtestForks=4          # cap the fork count
#   ./unit-test-all.sh --tests '*DevSurfaceV18*'
#
set -euo pipefail
cd "$(dirname "$0")"
exec ./gradlew test -PallTests -PtestForks=auto "$@"
