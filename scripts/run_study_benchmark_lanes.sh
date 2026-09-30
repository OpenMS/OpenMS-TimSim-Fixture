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

rm -rf "$RESULTS_OUT"
mkdir -p "$RESULTS_OUT"

printf '\n=== Multi-run identification / quantification benchmark ===\n'
python "$ROOT/tools/benchmark_study.py" \
  --build-dir "$BUILD_DIR" \
  --opendia-dir "$TARGET_OUT" \
  --out-dir "$RESULTS_OUT/study"

printf '\n=== External-null FDR benchmark ===\n'
python "$ROOT/tools/benchmark_entrapment.py" \
  --build-dir "$BUILD_DIR" \
  --opendia-dir "$ENTRAPMENT_OUT" \
  --out-dir "$RESULTS_OUT/entrapment"

printf '\nTarget-only OpenDIA: %s\n' "$TARGET_OUT"
printf 'Entrapment OpenDIA:  %s\n' "$ENTRAPMENT_OUT"
printf 'Study benchmark:     %s\n' "$RESULTS_OUT/study"
printf 'FDR benchmark:       %s\n' "$RESULTS_OUT/entrapment"
