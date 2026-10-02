#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=common_plan.sh
source "$ROOT/scripts/scale/common_plan.sh"
REFERENCE_ROOT="${1:?usage: plan_single_cell_calibration.sh SCALE_REFERENCE_ROOT STUDY_ROOT}"
STUDY_ROOT="${2:?missing study root}"
CONTROL_RUNS="${CONTROL_RUNS:-24}"
TREATMENT_RUNS="${TREATMENT_RUNS:-24}"
INPUT_LEVELS="${INPUT_LEVELS:--6,-8,-10,-12}"

BLUEPRINT_DB="$REFERENCE_ROOT/OpenSwathTimSim_blueprint/synthetic_data.db"
FASTA="$REFERENCE_ROOT/generated_inputs/synthetic_proteome.fasta"
SELECTED="$REFERENCE_ROOT/OpenSwathTimSim.high_signal_selection.tsv"
for path in "$BLUEPRINT_DB" "$FASTA" "$SELECTED"; do [[ -s "$path" ]] || { echo "ERROR: missing scale reference input: $path" >&2; exit 1; }; done

rm -rf "$STUDY_ROOT"
mkdir -p "$STUDY_ROOT"
stage_scale_reference "$REFERENCE_ROOT" "$STUDY_ROOT"

python "$ROOT/tools/plan_materialized_study.py" \
  --blueprint-db "$BLUEPRINT_DB" \
  --output-root "$STUDY_ROOT" \
  --selected-precursors "$STUDY_ROOT/OpenSwathTimSim.high_signal_selection.tsv" \
  --fasta "$FASTA" \
  --control-runs "$CONTROL_RUNS" \
  --treatment-runs "$TREATMENT_RUNS" \
  --profile single_cell \
  --input-log2-offset-levels="$INPUT_LEVELS" \
  --event-sampling poisson \
  --allow-zero-events \
  --run-log2-sd 0.10 \
  --cell-size-log2-sd 0.35 \
  --protein-log2-sd 0.20 \
  --peptide-log2-sd 0.15

printf '\nSingle-cell calibration plan ready: %s\nInput offsets: %s\n' "$STUDY_ROOT" "$INPUT_LEVELS"
printf 'Choose a production input offset only after calibration/OpenDIA depth is measured.\n'
