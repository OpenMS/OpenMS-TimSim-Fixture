#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

BUILD_DIR="${1:?Usage: $0 STUDY_DIR BENCH_DIR}"
BENCH_ROOT="${2:?Usage: $0 STUDY_DIR BENCH_DIR}"
TARGET_OUT="$BENCH_ROOT/opendia_target_only"
ENTRAPMENT_OUT="$BENCH_ROOT/opendia_entrapment"
RESULTS_OUT="$BENCH_ROOT/results"
PYTHON="$FIXTURE_VENV/bin/python"

[[ -x "$PYTHON" ]] || { echo "ERROR: run ./scripts/setup.sh first; missing $PYTHON" >&2; exit 1; }
[[ -s "$BUILD_DIR/fixture_manifest.json" ]] || { echo "ERROR: missing fixture manifest: $BUILD_DIR/fixture_manifest.json" >&2; exit 1; }
[[ -s "$TARGET_OUT/OpenDIA.results.tsv" ]] || { echo "ERROR: target-only OpenDIA results missing: $TARGET_OUT/OpenDIA.results.tsv" >&2; exit 1; }
[[ -s "$ENTRAPMENT_OUT/OpenDIA.results.tsv" ]] || { echo "ERROR: entrapment OpenDIA results missing: $ENTRAPMENT_OUT/OpenDIA.results.tsv" >&2; exit 1; }
[[ -s "$ENTRAPMENT_OUT/OpenDIA_intermediates/workflow.oswpq" ]] || { echo "ERROR: entrapment workflow archive missing: $ENTRAPMENT_OUT/OpenDIA_intermediates/workflow.oswpq" >&2; exit 1; }

rm -rf "$RESULTS_OUT"
mkdir -p "$RESULTS_OUT"

printf '=== Multi-run identification / quantification benchmark ===\n'
"$PYTHON" "$ROOT/tools/benchmark_study.py" \
  --build-dir "$BUILD_DIR" \
  --opendia-dir "$TARGET_OUT" \
  --out-dir "$RESULTS_OUT/study"

printf '\n=== Quantification diagnostics ===\n'
"$PYTHON" "$ROOT/tools/diagnose_study_quantification.py" \
  --study-results-dir "$RESULTS_OUT/study"

printf '\n=== External-null FDR benchmark ===\n'
"$PYTHON" "$ROOT/tools/benchmark_entrapment.py" \
  --build-dir "$BUILD_DIR" \
  --opendia-dir "$ENTRAPMENT_OUT" \
  --out-dir "$RESULTS_OUT/entrapment"

printf '\nStudy benchmark: %s\n' "$RESULTS_OUT/study"
printf 'FDR benchmark:   %s\n' "$RESULTS_OUT/entrapment"
