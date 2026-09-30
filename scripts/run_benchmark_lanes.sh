#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

BUILD_DIR="${1:-$DEFAULT_FIXTURE_BUILD}"
BENCH_ROOT="${2:-$DEFAULT_BENCHMARK_ROOT}"
THREADS="${THREADS:-12}"

TARGET_LIBRARY="$BUILD_DIR/OpenSwathTimSim.target_only.transitions.tsv"
ENTRAPMENT_LIBRARY="$BUILD_DIR/OpenSwathTimSim.transitions.tsv"
TARGET_OUT="$BENCH_ROOT/opendia_target_only"
ENTRAPMENT_OUT="$BENCH_ROOT/opendia_entrapment"
RESULTS_OUT="$BENCH_ROOT/results"

[[ -s "$TARGET_LIBRARY" ]] || { echo "ERROR: missing target-only library: $TARGET_LIBRARY" >&2; exit 1; }
[[ -s "$ENTRAPMENT_LIBRARY" ]] || { echo "ERROR: missing active entrapment library: $ENTRAPMENT_LIBRARY" >&2; exit 1; }

printf '=== OpenDIA target-only ID/quant lane ===\n'
THREADS="$THREADS" "$ROOT/scripts/run_opendia.sh" \
  "$BUILD_DIR" \
  "$TARGET_OUT" \
  "$TARGET_LIBRARY"

printf '\n=== OpenDIA external-null FDR lane ===\n'
THREADS="$THREADS" "$ROOT/scripts/run_opendia.sh" \
  "$BUILD_DIR" \
  "$ENTRAPMENT_OUT" \
  "$ENTRAPMENT_LIBRARY"

printf '\n=== Benchmarking isolated lanes ===\n'
"$ROOT/scripts/benchmark_fixture.sh" \
  "$BUILD_DIR" \
  "$TARGET_OUT" \
  "$RESULTS_OUT" \
  "$ENTRAPMENT_OUT"

printf '\nTarget-only OpenDIA: %s\n' "$TARGET_OUT"
printf 'Entrapment OpenDIA:  %s\n' "$ENTRAPMENT_OUT"
printf 'Benchmark results:  %s\n' "$RESULTS_OUT"
