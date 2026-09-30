#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

STUDY_DIR="${1:?Usage: $0 STUDY_DIR BENCH_DIR [ORACLE_OUT_DIR]}"
BENCH_ROOT="${2:?Usage: $0 STUDY_DIR BENCH_DIR [ORACLE_OUT_DIR]}"
ORACLE_OUT="${3:-$BENCH_ROOT/results/raw_signal_oracle}"
PYTHON="$FIXTURE_VENV/bin/python"

[[ -x "$PYTHON" ]] || { echo "ERROR: run ./scripts/setup.sh first; missing $PYTHON" >&2; exit 1; }
[[ -s "$STUDY_DIR/OpenSwathTimSim.realized_truth.tsv" ]] || { echo "ERROR: missing realized truth" >&2; exit 1; }
[[ -s "$STUDY_DIR/OpenSwathTimSim.target_only.transitions.tsv" ]] || { echo "ERROR: missing target-only transition library" >&2; exit 1; }
[[ -s "$BENCH_ROOT/results/study/precursor_measurements.tsv" ]] || {
  echo "ERROR: missing completed target-only study benchmark: $BENCH_ROOT/results/study/precursor_measurements.tsv" >&2
  echo "Run scripts/benchmark_study_lanes.sh first; OpenDIA itself does not need to be rerun." >&2
  exit 1
}
require_external_output_path "$ORACLE_OUT" "raw-signal oracle output"

if [[ "${RAW_ORACLE_PROBE_ONLY:-0}" == "1" ]]; then
  "$PYTHON" "$ROOT/tools/extract_raw_signal_oracle.py" \
    --study-dir "$STUDY_DIR" \
    --out-dir "$ORACLE_OUT" \
    --probe-reader
  exit 0
fi

rm -rf "$ORACLE_OUT"
mkdir -p "$ORACLE_OUT"

"$PYTHON" "$ROOT/tools/extract_raw_signal_oracle.py" \
  --study-dir "$STUDY_DIR" \
  --out-dir "$ORACLE_OUT" \
  --rt-half-window "${RAW_ORACLE_RT_HALF_WINDOW:-6.0}" \
  --im-half-window "${RAW_ORACLE_IM_HALF_WINDOW:-0.03}" \
  --fragment-ppm "${RAW_ORACLE_FRAGMENT_PPM:-25.0}"

"$PYTHON" "$ROOT/tools/benchmark_raw_signal_oracle.py" \
  --oracle-dir "$ORACLE_OUT" \
  --study-results-dir "$BENCH_ROOT/results/study" \
  --out-dir "$ORACLE_OUT"

printf '\nRaw-signal oracle: %s\n' "$ORACLE_OUT"
