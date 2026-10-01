#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck source=../common_paths.sh
source "$ROOT/scripts/common_paths.sh"

REFERENCE_D="${1:?usage: prepare_scale_reference.sh REFERENCE.d OUTPUT_ROOT}"
OUTPUT_ROOT="${2:?missing output root}"
PRECURSORS="${PRECURSORS:-150000}"
ENTRAPMENTS="${ENTRAPMENTS:-0}"
SIMULATED_PEPTIDES="${SIMULATED_PEPTIDES:-300000}"
FASTA_PEPTIDES="${FASTA_PEPTIDES:-600000}"
PEPTIDES_PER_PROTEIN="${PEPTIDES_PER_PROTEIN:-20}"
TIMSIM_THREADS="${TIMSIM_THREADS:-8}"
GRADIENT_LENGTH="${GRADIENT_LENGTH:-30}"
TRANSITIONS_PER_PRECURSOR="${TRANSITIONS_PER_PRECURSOR:-8}"
BLUEPRINT_SEED="${BLUEPRINT_SEED:-2026100201}"
FASTA_SEED="${FASTA_SEED:-1731}"
ENTRAPMENT_SEED="${ENTRAPMENT_SEED:-2026100202}"
COMPACT_REFERENCE="${COMPACT_REFERENCE:-1}"

REFERENCE_D="$(realpath "$REFERENCE_D")"
OUTPUT_ROOT="$(realpath -m "$OUTPUT_ROOT")"
require_external_output_path "$OUTPUT_ROOT" "scale reference output directory"
[[ -d "$REFERENCE_D" && -s "$REFERENCE_D/analysis.tdf" && -s "$REFERENCE_D/analysis.tdf_bin" ]] || {
  echo "ERROR: reference .d is incomplete: $REFERENCE_D" >&2; exit 2;
}
for value in "$PRECURSORS" "$SIMULATED_PEPTIDES" "$FASTA_PEPTIDES" "$TIMSIM_THREADS"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: count/thread options must be positive integers" >&2; exit 2; }
done
[[ "$ENTRAPMENTS" =~ ^[0-9]+$ ]] || { echo "ERROR: ENTRAPMENTS must be a non-negative integer" >&2; exit 2; }
(( SIMULATED_PEPTIDES >= PRECURSORS )) || { echo "ERROR: SIMULATED_PEPTIDES must be >= PRECURSORS" >&2; exit 2; }
(( FASTA_PEPTIDES >= SIMULATED_PEPTIDES )) || { echo "ERROR: FASTA_PEPTIDES must be >= SIMULATED_PEPTIDES" >&2; exit 2; }

VENV="$FIXTURE_VENV"
[[ -x "$VENV/bin/python" && -x "$VENV/bin/timsim" ]] || { echo "ERROR: generation environment unavailable under $VENV" >&2; exit 1; }
PYTHON="$VENV/bin/python"
TIMSIM="$VENV/bin/timsim"

"$PYTHON" - <<'PY'
import torch
if not torch.cuda.is_available():
    raise SystemExit("ERROR: scale reference requires GPU but torch.cuda.is_available() is false")
print(f"Scale-reference GPU: {torch.cuda.get_device_name(0)}; torch={torch.__version__}; CUDA={torch.version.cuda}")
PY

rm -rf "$OUTPUT_ROOT"
mkdir -p "$OUTPUT_ROOT"/{generated_inputs,rendered_configs,selection_qc}
FASTA="$OUTPUT_ROOT/generated_inputs/synthetic_proteome.fasta"
BLUEPRINT_NAME=OpenSwathTimSim_blueprint
BLUEPRINT_DB="$OUTPUT_ROOT/$BLUEPRINT_NAME/synthetic_data.db"
SELECTED="$OUTPUT_ROOT/OpenSwathTimSim.high_signal_selection.tsv"
N_PROTEINS=$(( (FASTA_PEPTIDES + PEPTIDES_PER_PROTEIN - 1) / PEPTIDES_PER_PROTEIN ))

"$PYTHON" "$ROOT/tools/generate_synthetic_fasta.py" \
  --out "$FASTA" --peptides "$FASTA_PEPTIDES" --seed "$FASTA_SEED" --peptides-per-protein "$PEPTIDES_PER_PROTEIN"

"$PYTHON" "$ROOT/tools/render_study_configs.py" \
  --repo-root "$ROOT" \
  --output-root "$OUTPUT_ROOT" \
  --rendered-dir "$OUTPUT_ROOT/rendered_configs" \
  --reference "$REFERENCE_D" \
  --fasta "$FASTA" \
  --blueprint-name "$BLUEPRINT_NAME" \
  --blueprint-sample-seed "$BLUEPRINT_SEED" \
  --n-proteins "$N_PROTEINS" \
  --num-peptides-total "$FASTA_PEPTIDES" \
  --num-sample-peptides "$SIMULATED_PEPTIDES" \
  --gradient-length "$GRADIENT_LENGTH" \
  --timsim-threads "$TIMSIM_THREADS" \
  --batch-size 128 --frame-batch-size 100 --use-gpu --blueprint-only

CONFIG="$OUTPUT_ROOT/rendered_configs/000_blueprint.toml"
if "$TIMSIM" --help 2>&1 | grep -q -- '--config'; then "$TIMSIM" --config "$CONFIG"; else "$TIMSIM" "$CONFIG"; fi
[[ -s "$BLUEPRINT_DB" ]] || { echo "ERROR: blueprint DB missing: $BLUEPRINT_DB" >&2; exit 1; }

"$PYTHON" "$ROOT/tools/select_reference_precursors.py" \
  --blueprint-db "$BLUEPRINT_DB" \
  --selected-out "$SELECTED" \
  --reference-truth-out "$OUTPUT_ROOT/OpenSwathTimSim.reference_truth.tsv" \
  --candidate-qc-out "$OUTPUT_ROOT/OpenSwathTimSim.selection_candidates.tsv" \
  --qc-json "$OUTPUT_ROOT/selection_qc/selection_summary.json" \
  --qc-report "$OUTPUT_ROOT/selection_qc/selection_report.md" \
  --precursors "$PRECURSORS" \
  --selection-mode global_stratified \
  --min-realized-event-proxy 50000 \
  --min-frame-abundance-sum 0.90 \
  --min-scan-abundance-sum 0.95 \
  --min-ion-relative-abundance 0.25 \
  --rt-edge-margin-fraction 0.05 \
  --rt-edge-min-seconds 6 \
  --minimum-fragments 8 \
  --mz-window-margin 1 \
  --im-window-margin 0.005

"$PYTHON" "$ROOT/tools/export_openswath_tsv.py" \
  --db "$BLUEPRINT_DB" \
  --selected-precursors "$SELECTED" \
  --out "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" \
  --truth-out "$OUTPUT_ROOT/OpenSwathTimSim.ground_truth.tsv" \
  --max-precursors "$PRECURSORS" --require-exact \
  --transitions-per-precursor "$TRANSITIONS_PER_PRECURSOR" --minimum-fragments 8 \
  --min-product-mz 350 --max-product-mz 2000 \
  --mz-window-margin 1 --im-window-margin 0.005 \
  --gradient-length "$GRADIENT_LENGTH" --rt-mode percent

# Entrapments do not affect raw generation. The existing small-fixture independent
# matcher intentionally solves a global donor assignment and is not used against
# the 150k donor population by default. Scale external-null generation is a
# separate downstream library step after raw-generation throughput is validated.
if (( ENTRAPMENTS > 0 )); then
  "$PYTHON" "$ROOT/tools/generate_entrapment_library.py" \
    --target-transitions "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" \
    --baseline-db "$BLUEPRINT_DB" \
    --out "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent.transitions.tsv" \
    --truth-out "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent_truth.tsv" \
    --mode independent --count "$ENTRAPMENTS" --seed "$ENTRAPMENT_SEED" \
    --fragment-collision-ppm 25 --max-paired-fragment-collisions 1 \
    --independent-min-rt-separation 2 --independent-min-im-separation 0.03 \
    --mz-window-margin 1 --im-window-margin 0.005
  cp "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent.transitions.tsv" "$OUTPUT_ROOT/OpenSwathTimSim.transitions.tsv"
  cp "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent_truth.tsv" "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_truth.tsv"
else
  cp "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" "$OUTPUT_ROOT/OpenSwathTimSim.transitions.tsv"
fi

"$PYTHON" - "$OUTPUT_ROOT" "$PRECURSORS" "$SIMULATED_PEPTIDES" "$FASTA_PEPTIDES" "$ENTRAPMENTS" <<'PY'
import csv, json, sys
from pathlib import Path
root=Path(sys.argv[1])
with (root/'OpenSwathTimSim.high_signal_selection.tsv').open(newline='', encoding='utf-8') as h:
    rows=list(csv.DictReader(h, delimiter='\t'))
counts={}
for row in rows: counts[row['protein_id']]=counts.get(row['protein_id'],0)+1
values=sorted(counts.values())
def pct(p):
    if not values: return None
    return values[min(len(values)-1, int(round((len(values)-1)*p)))]
payload={
  'reference_kind':'large_variable_proteome_v1',
  'targets':int(sys.argv[2]),
  'simulated_peptides_requested':int(sys.argv[3]),
  'fasta_peptides':int(sys.argv[4]),
  'independent_entrapments':int(sys.argv[5]),
  'selection_mode':'global_stratified_variable_proteome_contract',
  'selected_proteins':len(counts),
  'precursors_per_protein':{'min':min(values), 'median':pct(0.5), 'p90':pct(0.9), 'p95':pct(0.95), 'max':max(values)},
}
(root/'scale_reference_manifest.json').write_text(json.dumps(payload, indent=2)+'\n', encoding='utf-8')
print(json.dumps(payload, indent=2))
PY

if [[ "$COMPACT_REFERENCE" == "1" || "$COMPACT_REFERENCE" == "true" ]]; then
  BLUEPRINT_DIR="$OUTPUT_ROOT/$BLUEPRINT_NAME"
  find "$BLUEPRINT_DIR" -mindepth 1 -maxdepth 1 ! -name synthetic_data.db -exec rm -rf {} +
fi
printf '\nScale reference ready: %s\n' "$OUTPUT_ROOT"
printf 'Inspect %s/selection_qc/selection_summary.json and scale_reference_manifest.json before planning biological runs.\n' "$OUTPUT_ROOT"
