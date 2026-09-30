#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

BUILD_DIR="${1:-$DEFAULT_FIXTURE_BUILD}"
TARGET_OPENDIA_DIR="${2:-$DEFAULT_BENCHMARK_ROOT/opendia_target_only}"
OUT_ROOT="${3:-$DEFAULT_BENCHMARK_ROOT/results}"
ENTRAPMENT_OPENDIA_DIR="${4:-$TARGET_OPENDIA_DIR}"
PYTHON="$FIXTURE_VENV/bin/python"

require_external_output_path "$OUT_ROOT" "benchmark output directory"
[[ -x "$PYTHON" ]] || { echo "ERROR: run ./scripts/setup.sh first" >&2; exit 1; }
[[ -s "$BUILD_DIR/fixture_manifest.json" ]] || { echo "ERROR: fixture not found: $BUILD_DIR" >&2; exit 1; }
mkdir -p "$OUT_ROOT"

IDENT_DIR="$OUT_ROOT/identification_quantification"
CORRECTNESS_DIR="$OUT_ROOT/correctness"
ENTRAPMENT_DIR="$OUT_ROOT/entrapment"

"$PYTHON" "$ROOT/tools/benchmark_opendia.py" \
  --build-dir "$BUILD_DIR" \
  --opendia-dir "$TARGET_OPENDIA_DIR" \
  --out-dir "$IDENT_DIR" \
  --synthetic-db-root "$BUILD_DIR" \
  --realized-coordinate "${REALIZED_COORDINATE:-apex}"

correctness_args=(
  --build-dir "$BUILD_DIR"
  --benchmark-dir "$IDENT_DIR"
  --out-dir "$CORRECTNESS_DIR"
)
[[ "${REQUIRE_PERFECT_RECOVERY:-false}" == "true" ]] && correctness_args+=(--require-perfect-recovery)
"$PYTHON" "$ROOT/tools/benchmark_correctness.py" "${correctness_args[@]}"

if [[ -s "$BUILD_DIR/OpenSwathTimSim.entrapment_truth.tsv" ]]; then
  entrapment_args=(
    --build-dir "$BUILD_DIR"
    --opendia-dir "$ENTRAPMENT_OPENDIA_DIR"
    --out-dir "$ENTRAPMENT_DIR"
    --q-threshold "${ENTRAPMENT_Q_THRESHOLD:-0.01}"
  )
  [[ -n "${OPENDIA_WORKFLOW_OSWPQ:-}" ]] && entrapment_args+=(--workflow-oswpq "$OPENDIA_WORKFLOW_OSWPQ")
  [[ -n "${OPENDIA_RESULTS_TSV:-}" ]] && entrapment_args+=(--results-tsv "$OPENDIA_RESULTS_TSV")
  "$PYTHON" "$ROOT/tools/benchmark_entrapment.py" "${entrapment_args[@]}"
else
  echo "No entrapment truth found; skipping entrapment benchmark."
fi

printf '\nBenchmark results: %s\n' "$OUT_ROOT"
