#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

BUILD_DIR="${1:?Usage: $0 STUDY_DIR BENCH_DIR}"
BENCH_ROOT="${2:?Usage: $0 STUDY_DIR BENCH_DIR}"
THREADS="${THREADS:-12}"
TARGET_LIBRARY="$BUILD_DIR/OpenSwathTimSim.target_only.transitions.tsv"
ENTRAPMENT_LIBRARY="$BUILD_DIR/OpenSwathTimSim.transitions.tsv"
TARGET_OUT="$BENCH_ROOT/opendia_target_only"
ENTRAPMENT_OUT="$BENCH_ROOT/opendia_entrapment"
RESULTS_OUT="$BENCH_ROOT/results"

[[ -s "$BUILD_DIR/fixture_manifest.json" ]] || { echo "ERROR: missing fixture manifest" >&2; exit 1; }
[[ -s "$TARGET_LIBRARY" ]] || { echo "ERROR: missing target-only library: $TARGET_LIBRARY" >&2; exit 1; }
[[ -s "$ENTRAPMENT_LIBRARY" ]] || { echo "ERROR: missing entrapment library: $ENTRAPMENT_LIBRARY" >&2; exit 1; }

printf '=== OpenDIA target-only study lane ===\n'
THREADS="$THREADS" "$ROOT/scripts/run_opendia.sh" "$BUILD_DIR" "$TARGET_OUT" "$TARGET_LIBRARY"

printf '\n=== OpenDIA independent-entrapment study lane ===\n'
THREADS="$THREADS" "$ROOT/scripts/run_opendia.sh" "$BUILD_DIR" "$ENTRAPMENT_OUT" "$ENTRAPMENT_LIBRARY"

printf '
=== Benchmarking completed OpenDIA study lanes ===
'
"$ROOT/scripts/benchmark_study_lanes.sh" "$BUILD_DIR" "$BENCH_ROOT"

printf '
Target-only OpenDIA: %s
' "$TARGET_OUT"
printf 'Entrapment OpenDIA:  %s
' "$ENTRAPMENT_OUT"
printf 'Study benchmark:     %s
' "$RESULTS_OUT/study"
printf 'FDR benchmark:       %s
' "$RESULTS_OUT/entrapment"
