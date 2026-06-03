#!/usr/bin/env bash
#
# Run the DEFAULT unit-test suite (the experimental DevSurfaceV* ladder contract tests and the
# cross-variant equivalence harness are excluded - see lib/build.gradle), single-threaded.
#
# Extra arguments are forwarded to Gradle, e.g.:
#   ./unit-test.sh --tests '*DistinctPackedNumericalSurfaceV2Test'
#
set -euo pipefail
cd "$(dirname "$0")"
exec ./gradlew test -PtestForks=1 "$@"
