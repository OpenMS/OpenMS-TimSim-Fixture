#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

STUDY_DIR="${1:?Usage: $0 STUDY_DIR BENCH_DIR [SWEEP_OUT_DIR]}"
BENCH_ROOT="${2:?Usage: $0 STUDY_DIR BENCH_DIR [SWEEP_OUT_DIR]}"
SWEEP_OUT="${3:-$BENCH_ROOT/results/raw_signal_oracle_sweep}"
PYTHON="$FIXTURE_VENV/bin/python"

[[ -x "$PYTHON" ]] || { echo "ERROR: run ./scripts/setup.sh first; missing $PYTHON" >&2; exit 1; }
[[ -s "$STUDY_DIR/OpenSwathTimSim.realized_truth.tsv" ]] || { echo "ERROR: missing realized truth" >&2; exit 1; }
[[ -s "$STUDY_DIR/OpenSwathTimSim.target_only.transitions.tsv" ]] || { echo "ERROR: missing target-only transition library" >&2; exit 1; }
[[ -s "$STUDY_DIR/OpenSwathTimSim_blueprint/synthetic_data.db" ]] || { echo "ERROR: missing TimSim blueprint database" >&2; exit 1; }
[[ -s "$BENCH_ROOT/results/study/precursor_measurements.tsv" ]] || {
  echo "ERROR: missing completed target-only study benchmark" >&2
  echo "Run scripts/benchmark_study_lanes.sh first; TimSim and OpenDIA do not need to be rerun." >&2
  exit 1
}
require_external_output_path "$SWEEP_OUT" "raw-oracle sweep output"

rm -rf "$SWEEP_OUT"
mkdir -p "$SWEEP_OUT"

"$PYTHON" "$ROOT/tools/sweep_raw_signal_oracle.py" \
  --study-dir "$STUDY_DIR" \
  --study-results-dir "$BENCH_ROOT/results/study" \
  --out-dir "$SWEEP_OUT" \
  --rt-windows "${RAW_ORACLE_SWEEP_RT_WINDOWS:-1,2,4,6}" \
  --im-windows "${RAW_ORACLE_SWEEP_IM_WINDOWS:-0.01,0.02,0.03}" \
  --fragment-ppm "${RAW_ORACLE_SWEEP_FRAGMENT_PPM:-10,15,25}" \
  --transition-subsets "${RAW_ORACLE_SWEEP_TRANSITION_SUBSETS:-8,6,4,3}"

printf '\nRaw-oracle interference sweep: %s\n' "$SWEEP_OUT"
