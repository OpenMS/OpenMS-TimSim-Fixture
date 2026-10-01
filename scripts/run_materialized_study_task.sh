#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STUDY_ROOT="${1:?usage: run_materialized_study_task.sh STUDY_ROOT RUN_ORDINAL REFERENCE_D [SCRATCH_ROOT]}"
RUN_ORDINAL="${2:?missing run ordinal}"
REFERENCE_D="${3:?missing reference .d}"
SCRATCH_BASE="${4:-${TMPDIR:-/tmp}}"

[[ "$RUN_ORDINAL" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: RUN_ORDINAL must be positive" >&2; exit 2; }
[[ -d "$STUDY_ROOT" ]] || { echo "ERROR: study root missing: $STUDY_ROOT" >&2; exit 2; }
[[ -d "$REFERENCE_D" && -s "$REFERENCE_D/analysis.tdf" && -s "$REFERENCE_D/analysis.tdf_bin" ]] || {
  echo "ERROR: reference .d incomplete: $REFERENCE_D" >&2; exit 2;
}

PLAN="$STUDY_ROOT/OpenSwathTimSim.materialized_study_plan.json"
MANIFEST="$STUDY_ROOT/OpenSwathTimSim.study_manifest.tsv"
DESIGN="$STUDY_ROOT/OpenSwathTimSim.study_design_truth.tsv"
[[ -s "$PLAN" && -s "$MANIFEST" && -s "$DESIGN" ]] || { echo "ERROR: study plan artifacts missing" >&2; exit 1; }

readarray -t META < <(python - "$PLAN" "$MANIFEST" "$RUN_ORDINAL" <<'PY'
import csv, json, sys
plan=json.load(open(sys.argv[1], encoding='utf-8'))
ordinal=int(sys.argv[3])
with open(sys.argv[2], newline='', encoding='utf-8') as h:
    rows=[r for r in csv.DictReader(h, delimiter='\t') if int(r['RunOrdinal']) == ordinal]
if len(rows) != 1:
    raise SystemExit(f'expected one manifest row for ordinal {ordinal}, found {len(rows)}')
r=rows[0]
print(plan['blueprint_db'])
print(plan['selected_precursors'])
print(plan['fasta'])
print(r['RunName'])
print(r['RunId'])
print(r['Condition'])
PY
)
BLUEPRINT_DB="${META[0]}"
SELECTED="${META[1]}"
FASTA="${META[2]}"
RUN_NAME="${META[3]}"
RUN_ID="${META[4]}"
CONDITION="${META[5]}"

for path in "$BLUEPRINT_DB" "$SELECTED" "$FASTA"; do
  [[ -e "$path" ]] || { echo "ERROR: plan input missing: $path" >&2; exit 1; }
done

TIMSIM_THREADS="${TIMSIM_THREADS:-8}"
BATCH_SIZE="${BATCH_SIZE:-128}"
FRAME_BATCH_SIZE="${FRAME_BATCH_SIZE:-100}"
USE_GPU="${USE_GPU:-1}"
COMPACT_OUTPUT="${COMPACT_OUTPUT:-1}"

TASK_SCRATCH="$SCRATCH_BASE/openms-timsim-${RUN_NAME}-${RUN_ORDINAL}"
SOURCE_DIR="$TASK_SCRATCH/source"
CONFIG="$TASK_SCRATCH/${RUN_ID}.toml"
INPUT_TRUTH="$STUDY_ROOT/run_input_truth/$(printf '%03d' "$RUN_ORDINAL")_${RUN_ID}.tsv.gz"
REALIZED_TRUTH="$STUDY_ROOT/run_realized_truth/$(printf '%03d' "$RUN_ORDINAL")_${RUN_ID}.tsv.gz"
STATS="$STUDY_ROOT/run_stats/$(printf '%03d' "$RUN_ORDINAL")_${RUN_ID}.json"
PROVENANCE="$STUDY_ROOT/run_provenance/$(printf '%03d' "$RUN_ORDINAL")_${RUN_ID}.txt"
RUN_DIR="$STUDY_ROOT/$RUN_NAME"
OUT_DB="$RUN_DIR/synthetic_data.db"
OUT_D="$RUN_DIR/$RUN_NAME.d"

rm -rf "$TASK_SCRATCH" "$RUN_DIR"
mkdir -p "$SOURCE_DIR" "$STUDY_ROOT"/{run_input_truth,run_realized_truth,run_stats,run_provenance}

cleanup() { rm -rf "$TASK_SCRATCH"; }
trap cleanup EXIT

python "$ROOT/tools/materialize_study_run.py" \
  --blueprint-db "$BLUEPRINT_DB" \
  --study-manifest "$MANIFEST" \
  --design-truth "$DESIGN" \
  --selected-precursors "$SELECTED" \
  --run-ordinal "$RUN_ORDINAL" \
  --output-db "$SOURCE_DIR/synthetic_data.db" \
  --selected-input-truth-out "$INPUT_TRUTH" \
  --stats-out "$STATS"

GPU_ARGS=()
if [[ "$USE_GPU" == "1" || "$USE_GPU" == "true" ]]; then
  GPU_ARGS+=(--use-gpu)
fi
python "$ROOT/tools/render_materialized_run_config.py" \
  --repo-root "$ROOT" \
  --study-root "$STUDY_ROOT" \
  --study-manifest "$MANIFEST" \
  --plan "$PLAN" \
  --run-ordinal "$RUN_ORDINAL" \
  --existing-path "$SOURCE_DIR" \
  --reference "$REFERENCE_D" \
  --fasta "$FASTA" \
  --out "$CONFIG" \
  --timsim-threads "$TIMSIM_THREADS" \
  --batch-size "$BATCH_SIZE" \
  --frame-batch-size "$FRAME_BATCH_SIZE" \
  "${GPU_ARGS[@]}"

if [[ -z "${FIXTURE_VENV:-}" && -x "$ROOT/.venv/bin/python" ]]; then
  export FIXTURE_VENV="$ROOT/.venv"
fi
"$ROOT/scripts/run_timsim_config.sh" "$CONFIG"

[[ -s "$OUT_DB" && -d "$OUT_D" && -s "$OUT_D/analysis.tdf" && -s "$OUT_D/analysis.tdf_bin" ]] || {
  echo "ERROR: generated run incomplete: $RUN_NAME" >&2
  exit 1
}

python "$ROOT/tools/extract_materialized_run_truth.py" \
  --run-db "$OUT_DB" \
  --selected-precursors "$SELECTED" \
  --selected-input-truth "$INPUT_TRUTH" \
  --condition "$CONDITION" \
  --run-name "$RUN_NAME" \
  --out "$REALIZED_TRUTH"

# Input truth is fully represented in the realized partition after extraction.
rm -f "$INPUT_TRUTH"

{
  printf 'run_ordinal=%s\nrun_id=%s\nrun_name=%s\n' "$RUN_ORDINAL" "$RUN_ID" "$RUN_NAME"
  printf 'host=%s\n' "$(hostname -f 2>/dev/null || hostname)"
  printf 'started_from_blueprint=%s\n' "$BLUEPRINT_DB"
  printf 'reference_d=%s\n' "$REFERENCE_D"
  printf 'timsim_threads=%s\nuse_gpu=%s\ncompact_output=%s\n' "$TIMSIM_THREADS" "$USE_GPU" "$COMPACT_OUTPUT"
  printf 'slurm_job_id=%s\nslurm_array_task_id=%s\n' "${SLURM_JOB_ID:-}" "${SLURM_ARRAY_TASK_ID:-}"
  printf 'finished_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$PROVENANCE"

if [[ "$COMPACT_OUTPUT" == "1" || "$COMPACT_OUTPUT" == "true" ]]; then
  find "$RUN_DIR" -mindepth 1 -maxdepth 1 ! -name "$RUN_NAME.d" -exec rm -rf {} +
fi

printf 'PASS run=%s raw=%s truth=%s\n' "$RUN_NAME" "$OUT_D" "$REALIZED_TRUTH"
