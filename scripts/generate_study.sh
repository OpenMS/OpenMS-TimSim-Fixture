#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

REFERENCE_D=""
OUTPUT_ROOT="$OPENMS_TIMSIM_DATA_ROOT/studies/control25_treatment25_targets1000"
CONTROL_RUNS=25
TREATMENT_RUNS=25
PRECURSORS=1000
TARGET_PROTEINS=250
PRECURSORS_PER_PROTEIN=4
SELECTION_MODE=protein_balanced
ENTRAPMENTS=1000
SIMULATED_PEPTIDES=10000
FASTA_PEPTIDES=20000
PEPTIDES_PER_PROTEIN=20
GRADIENT_LENGTH=30
TRANSITIONS_PER_PRECURSOR=8
MINIMUM_FRAGMENTS=8
BLUEPRINT_SEED=2026093001
FASTA_SEED=1729
STUDY_SEED=2026093101
SAMPLE_SEED_BASE=2026094000
ENTRAPMENT_SEED=2026093201
TIMSIM_THREADS=1
BATCH_SIZE=128
FRAME_BATCH_SIZE=100
RUN_LOG2_SD=0.05
PROTEIN_LOG2_SD=0.15
PEPTIDE_LOG2_SD=0.08
TREATMENT_UP_FRACTION=0.15
TREATMENT_DOWN_FRACTION=0.15
TREATMENT_MIN_ABS_LOG2FC=0.5
TREATMENT_MAX_ABS_LOG2FC=2.0
MIN_REALIZED_EVENT_PROXY=50000
MIN_FRAME_ABUNDANCE_SUM=0.90
MIN_SCAN_ABUNDANCE_SUM=0.95
MIN_ION_RELATIVE_ABUNDANCE=0.25
RT_EDGE_MARGIN_FRACTION=0.05
RT_EDGE_MIN_SECONDS=6.0
MZ_WINDOW_MARGIN=1.0
IM_WINDOW_MARGIN=0.005
ENTRAPMENT_MODE=independent
VENV="$FIXTURE_VENV"

usage() {
  cat >&2 <<EOF
Usage:
  $0 --reference /path/to/reference.d [options]

Generate a reproducible multi-run TimSim study. Targets are frozen from a
blueprint simulation before condition effects are applied. Control/treatment
runs then reuse the same molecular blueprint with independent abundance effects.

Study options:
  --output PATH                       Output directory
  --control-runs N                    Control biological replicates (default: 25)
  --treatment-runs N                  Treatment biological replicates (default: 25)
  --precursors N                      Frozen true precursor groups (default: 1000)
  --target-proteins N                 Protein-balanced target proteins (default: 250)
  --precursors-per-protein N          Frozen precursors per selected protein (default: 4)
  --selection-mode MODE               protein_balanced (default) or global_stratified
  --entrapments N                     Independent external-null precursors (default: 1000)
  --simulated-peptides N              Blueprint TimSim peptide population (default: 10000)
  --fasta-peptides N                  Synthetic FASTA peptide pool (default: 20000)
  --blueprint-seed N                  Blueprint TimSim sample seed
  --fasta-seed N                      Synthetic proteome seed
  --study-seed N                      Protein effects / replicate abundance seed
  --sample-seed-base N                First child TimSim seed base
  --timsim-threads N                  TimSim threads per run (default: 1)

Biological variation:
  --run-log2-sd X                     Run-global abundance SD (default: 0.05)
  --protein-log2-sd X                 Per-run protein abundance SD (default: 0.15)
  --peptide-log2-sd X                 Per-run peptide residual SD (default: 0.08)
  --treatment-up-fraction X           Up-regulated protein fraction (default: 0.15)
  --treatment-down-fraction X         Down-regulated protein fraction (default: 0.15)
  --treatment-min-abs-log2fc X        Minimum |design log2FC| (default: 0.5)
  --treatment-max-abs-log2fc X        Maximum |design log2FC| (default: 2.0)

Target-selection thresholds are applied to the blueprint only. If fewer than N
precursors qualify, increase --simulated-peptides/--fasta-peptides rather than
weakening the correctness criteria.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference) REFERENCE_D="${2:?missing value}"; shift 2 ;;
    --output) OUTPUT_ROOT="${2:?missing value}"; shift 2 ;;
    --control-runs) CONTROL_RUNS="${2:?missing value}"; shift 2 ;;
    --treatment-runs) TREATMENT_RUNS="${2:?missing value}"; shift 2 ;;
    --precursors) PRECURSORS="${2:?missing value}"; shift 2 ;;
    --target-proteins) TARGET_PROTEINS="${2:?missing value}"; shift 2 ;;
    --precursors-per-protein) PRECURSORS_PER_PROTEIN="${2:?missing value}"; shift 2 ;;
    --selection-mode) SELECTION_MODE="${2:?missing value}"; shift 2 ;;
    --entrapments) ENTRAPMENTS="${2:?missing value}"; shift 2 ;;
    --simulated-peptides) SIMULATED_PEPTIDES="${2:?missing value}"; shift 2 ;;
    --fasta-peptides) FASTA_PEPTIDES="${2:?missing value}"; shift 2 ;;
    --blueprint-seed) BLUEPRINT_SEED="${2:?missing value}"; shift 2 ;;
    --fasta-seed) FASTA_SEED="${2:?missing value}"; shift 2 ;;
    --study-seed) STUDY_SEED="${2:?missing value}"; shift 2 ;;
    --sample-seed-base) SAMPLE_SEED_BASE="${2:?missing value}"; shift 2 ;;
    --timsim-threads) TIMSIM_THREADS="${2:?missing value}"; shift 2 ;;
    --run-log2-sd) RUN_LOG2_SD="${2:?missing value}"; shift 2 ;;
    --protein-log2-sd) PROTEIN_LOG2_SD="${2:?missing value}"; shift 2 ;;
    --peptide-log2-sd) PEPTIDE_LOG2_SD="${2:?missing value}"; shift 2 ;;
    --treatment-up-fraction) TREATMENT_UP_FRACTION="${2:?missing value}"; shift 2 ;;
    --treatment-down-fraction) TREATMENT_DOWN_FRACTION="${2:?missing value}"; shift 2 ;;
    --treatment-min-abs-log2fc) TREATMENT_MIN_ABS_LOG2FC="${2:?missing value}"; shift 2 ;;
    --treatment-max-abs-log2fc) TREATMENT_MAX_ABS_LOG2FC="${2:?missing value}"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 2 ;;
  esac
done

[[ -n "$REFERENCE_D" ]] || { usage; exit 2; }
for value in "$CONTROL_RUNS" "$TREATMENT_RUNS" "$PRECURSORS" "$TARGET_PROTEINS" "$PRECURSORS_PER_PROTEIN" "$ENTRAPMENTS" "$SIMULATED_PEPTIDES" "$FASTA_PEPTIDES" "$TIMSIM_THREADS"; do
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: run/count/thread options must be positive integers" >&2; exit 2; }
done
for seed in "$BLUEPRINT_SEED" "$FASTA_SEED" "$STUDY_SEED" "$SAMPLE_SEED_BASE" "$ENTRAPMENT_SEED"; do
  [[ "$seed" =~ ^[0-9]+$ ]] || { echo "ERROR: seeds must be non-negative integers" >&2; exit 2; }
  (( seed <= 4294967295 )) || { echo "ERROR: seed exceeds uint32 range: $seed" >&2; exit 2; }
done
[[ "$SELECTION_MODE" == "protein_balanced" || "$SELECTION_MODE" == "global_stratified" ]] || {
  echo "ERROR: --selection-mode must be protein_balanced or global_stratified" >&2
  exit 2
}
if [[ "$SELECTION_MODE" == "protein_balanced" ]]; then
  EXPECTED_BALANCED_PRECURSORS=$(( TARGET_PROTEINS * PRECURSORS_PER_PROTEIN ))
  (( PRECURSORS == EXPECTED_BALANCED_PRECURSORS )) || {
    echo "ERROR: protein_balanced selection requires --precursors == --target-proteins * --precursors-per-protein ($PRECURSORS != $EXPECTED_BALANCED_PRECURSORS)" >&2
    exit 2
  }
fi
(( SIMULATED_PEPTIDES >= PRECURSORS )) || { echo "ERROR: --simulated-peptides must be >= --precursors" >&2; exit 2; }
(( FASTA_PEPTIDES >= SIMULATED_PEPTIDES )) || { echo "ERROR: --fasta-peptides must be >= --simulated-peptides" >&2; exit 2; }
[[ -f "$VENV/bin/activate" ]] || { echo "ERROR: run ./scripts/setup.sh first" >&2; exit 1; }

source "$VENV/bin/activate"
REFERENCE_D="$(realpath "$REFERENCE_D")"
OUTPUT_ROOT="$(realpath -m "$OUTPUT_ROOT")"
require_external_output_path "$OUTPUT_ROOT" "study output directory"
[[ -d "$REFERENCE_D" ]] || { echo "ERROR: reference .d directory does not exist: $REFERENCE_D" >&2; exit 1; }
[[ -s "$REFERENCE_D/analysis.tdf" && -s "$REFERENCE_D/analysis.tdf_bin" ]] || { echo "ERROR: reference .d is missing analysis.tdf/analysis.tdf_bin" >&2; exit 1; }

run_timsim() {
  local config="$1"
  if "$VENV/bin/timsim" --help 2>&1 | grep -q -- '--config'; then
    "$VENV/bin/timsim" --config "$config"
  else
    "$VENV/bin/timsim" "$config"
  fi
}

rm -rf "$OUTPUT_ROOT"
mkdir -p "$OUTPUT_ROOT"/{generated_inputs,rendered_configs,selection_qc,condition_inputs,tdfs}
FASTA="$OUTPUT_ROOT/generated_inputs/synthetic_proteome.fasta"
STUDY_MANIFEST="$OUTPUT_ROOT/OpenSwathTimSim.study_manifest.tsv"
DESIGN_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.study_design_truth.tsv"
ABUNDANCE_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.study_abundance_truth.tsv"
SELECTED_TSV="$OUTPUT_ROOT/OpenSwathTimSim.high_signal_selection.tsv"
REFERENCE_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.reference_truth.tsv"
REALIZED_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.realized_truth.tsv"
CANDIDATE_QC="$OUTPUT_ROOT/OpenSwathTimSim.selection_candidates.tsv"
QC_JSON="$OUTPUT_ROOT/selection_qc/selection_summary.json"
QC_REPORT="$OUTPUT_ROOT/selection_qc/selection_report.md"
BLUEPRINT_NAME="OpenSwathTimSim_blueprint"
BLUEPRINT_DB="$OUTPUT_ROOT/$BLUEPRINT_NAME/synthetic_data.db"
N_PROTEINS=$(( (FASTA_PEPTIDES + PEPTIDES_PER_PROTEIN - 1) / PEPTIDES_PER_PROTEIN ))

python "$ROOT/tools/generate_synthetic_fasta.py" \
  --out "$FASTA" \
  --peptides "$FASTA_PEPTIDES" \
  --seed "$FASTA_SEED" \
  --peptides-per-protein "$PEPTIDES_PER_PROTEIN"

python "$ROOT/tools/render_study_configs.py" \
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
  --batch-size "$BATCH_SIZE" \
  --frame-batch-size "$FRAME_BATCH_SIZE" \
  --blueprint-only

printf '\n=== TimSim blueprint ===\n'
run_timsim "$OUTPUT_ROOT/rendered_configs/000_blueprint.toml"
[[ -f "$BLUEPRINT_DB" ]] || { echo "ERROR: blueprint DB not found: $BLUEPRINT_DB" >&2; exit 1; }

# Freeze target identities before any control/treatment abundance is created.
python "$ROOT/tools/select_reference_precursors.py" \
  --blueprint-db "$BLUEPRINT_DB" \
  --selected-out "$SELECTED_TSV" \
  --reference-truth-out "$REFERENCE_TRUTH" \
  --candidate-qc-out "$CANDIDATE_QC" \
  --qc-json "$QC_JSON" \
  --qc-report "$QC_REPORT" \
  --precursors "$PRECURSORS" \
  --selection-mode "$SELECTION_MODE" \
  --target-proteins "$TARGET_PROTEINS" \
  --precursors-per-protein "$PRECURSORS_PER_PROTEIN" \
  --min-realized-event-proxy "$MIN_REALIZED_EVENT_PROXY" \
  --min-frame-abundance-sum "$MIN_FRAME_ABUNDANCE_SUM" \
  --min-scan-abundance-sum "$MIN_SCAN_ABUNDANCE_SUM" \
  --min-ion-relative-abundance "$MIN_ION_RELATIVE_ABUNDANCE" \
  --rt-edge-margin-fraction "$RT_EDGE_MARGIN_FRACTION" \
  --rt-edge-min-seconds "$RT_EDGE_MIN_SECONDS" \
  --minimum-fragments "$MINIMUM_FRAGMENTS" \
  --mz-window-margin "$MZ_WINDOW_MARGIN" \
  --im-window-margin "$IM_WINDOW_MARGIN"

python "$ROOT/tools/export_openswath_tsv.py" \
  --db "$BLUEPRINT_DB" \
  --selected-precursors "$SELECTED_TSV" \
  --out "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" \
  --truth-out "$OUTPUT_ROOT/OpenSwathTimSim.ground_truth.tsv" \
  --max-precursors "$PRECURSORS" \
  --require-exact \
  --transitions-per-precursor "$TRANSITIONS_PER_PRECURSOR" \
  --minimum-fragments "$MINIMUM_FRAGMENTS" \
  --min-product-mz 350 \
  --max-product-mz 2000 \
  --mz-window-margin "$MZ_WINDOW_MARGIN" \
  --im-window-margin "$IM_WINDOW_MARGIN" \
  --gradient-length "$GRADIENT_LENGTH" \
  --rt-mode percent

python "$ROOT/tools/generate_entrapment_library.py" \
  --target-transitions "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" \
  --baseline-db "$BLUEPRINT_DB" \
  --out "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent.transitions.tsv" \
  --truth-out "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent_truth.tsv" \
  --mode "$ENTRAPMENT_MODE" \
  --count "$ENTRAPMENTS" \
  --seed "$ENTRAPMENT_SEED" \
  --fragment-collision-ppm 25 \
  --max-paired-fragment-collisions 1 \
  --independent-min-rt-separation 2.0 \
  --independent-min-im-separation 0.03 \
  --mz-window-margin "$MZ_WINDOW_MARGIN" \
  --im-window-margin "$IM_WINDOW_MARGIN"
cp "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent.transitions.tsv" "$OUTPUT_ROOT/OpenSwathTimSim.transitions.tsv"
cp "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_independent_truth.tsv" "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_truth.tsv"

python "$ROOT/tools/prepare_study_databases.py" \
  --blueprint-db "$BLUEPRINT_DB" \
  --output-root "$OUTPUT_ROOT/condition_inputs" \
  --study-manifest-out "$STUDY_MANIFEST" \
  --design-truth-out "$DESIGN_TRUTH" \
  --abundance-truth-out "$ABUNDANCE_TRUTH" \
  --control-runs "$CONTROL_RUNS" \
  --treatment-runs "$TREATMENT_RUNS" \
  --study-seed "$STUDY_SEED" \
  --sample-seed-base "$SAMPLE_SEED_BASE" \
  --run-log2-sd "$RUN_LOG2_SD" \
  --protein-log2-sd "$PROTEIN_LOG2_SD" \
  --peptide-log2-sd "$PEPTIDE_LOG2_SD" \
  --treatment-up-fraction "$TREATMENT_UP_FRACTION" \
  --treatment-down-fraction "$TREATMENT_DOWN_FRACTION" \
  --treatment-min-abs-log2fc "$TREATMENT_MIN_ABS_LOG2FC" \
  --treatment-max-abs-log2fc "$TREATMENT_MAX_ABS_LOG2FC"

python "$ROOT/tools/render_study_configs.py" \
  --repo-root "$ROOT" \
  --output-root "$OUTPUT_ROOT" \
  --rendered-dir "$OUTPUT_ROOT/rendered_configs" \
  --reference "$REFERENCE_D" \
  --fasta "$FASTA" \
  --study-manifest "$STUDY_MANIFEST" \
  --blueprint-name "$BLUEPRINT_NAME" \
  --blueprint-sample-seed "$BLUEPRINT_SEED" \
  --n-proteins "$N_PROTEINS" \
  --num-peptides-total "$FASTA_PEPTIDES" \
  --num-sample-peptides "$SIMULATED_PEPTIDES" \
  --gradient-length "$GRADIENT_LENGTH" \
  --timsim-threads "$TIMSIM_THREADS" \
  --batch-size "$BATCH_SIZE" \
  --frame-batch-size "$FRAME_BATCH_SIZE"

printf '\n=== TimSim biological study runs ===\n'
mapfile -t RUN_CONFIGS < <(find "$OUTPUT_ROOT/rendered_configs" -maxdepth 1 -type f -name '[0-9][0-9][0-9]_[CT][0-9][0-9].toml' | sort)
EXPECTED_RUNS=$(( CONTROL_RUNS + TREATMENT_RUNS ))
(( ${#RUN_CONFIGS[@]} == EXPECTED_RUNS )) || { echo "ERROR: rendered ${#RUN_CONFIGS[@]} study configs; expected $EXPECTED_RUNS" >&2; exit 1; }
for config in "${RUN_CONFIGS[@]}"; do
  printf '\n--- %s ---\n' "$(basename "$config")"
  run_timsim "$config"
done

python "$ROOT/tools/build_study_realized_truth.py" \
  --fixture-root "$OUTPUT_ROOT" \
  --study-manifest "$STUDY_MANIFEST" \
  --selected-precursors "$SELECTED_TSV" \
  --abundance-truth "$ABUNDANCE_TRUTH" \
  --out "$REALIZED_TRUTH"

rm -rf "$OUTPUT_ROOT/tdfs"
mkdir -p "$OUTPUT_ROOT/tdfs"
mapfile -t RUN_NAMES < <(python - "$STUDY_MANIFEST" <<'PY'
import csv, sys
with open(sys.argv[1], newline='', encoding='utf-8') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        print(row['RunName'])
PY
)
for name in "${RUN_NAMES[@]}"; do
  source_d="$OUTPUT_ROOT/$name/$name.d"
  [[ -d "$source_d" ]] || { echo "ERROR: missing generated run TDF: $source_d" >&2; exit 1; }
  ln -s "../$name/$name.d" "$OUTPUT_ROOT/tdfs/$name.d"
done

python - "$OUTPUT_ROOT/fixture_manifest.json" <<PY
import hashlib, importlib.metadata, json, platform, sys
from pathlib import Path

path = Path(sys.argv[1])
root = Path(${OUTPUT_ROOT@Q})
# Read the authoritative TSV instead of reconstructing run order in Bash.
import csv
with open(root / "OpenSwathTimSim.study_manifest.tsv", newline="", encoding="utf-8") as handle:
    run_rows = list(csv.DictReader(handle, delimiter="\t"))
run_names = [row["RunName"] for row in run_rows]

def sha256_file(file_path):
    digest = hashlib.sha256()
    with open(file_path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

manifest = {
    "fixture_kind": "openms_timsim_multirun_study_v1",
    "selection_frozen_before_condition_effects": True,
    "selection_frozen_before_opendia": True,
    "software": {
        "python": platform.python_version(),
        "imspy-simulation": importlib.metadata.version("imspy-simulation"),
        "imspy-connector": importlib.metadata.version("imspy-connector"),
    },
    "reference_d": ${REFERENCE_D@Q},
    "blueprint_run": ${BLUEPRINT_NAME@Q},
    "requested_precursors": int(${PRECURSORS@Q}),
    "simulated_peptides": int(${SIMULATED_PEPTIDES@Q}),
    "fasta_peptides": int(${FASTA_PEPTIDES@Q}),
    "synthetic_proteins": int(${N_PROTEINS@Q}),
    "gradient_length_seconds": float(${GRADIENT_LENGTH@Q}),
    "transitions_per_precursor": int(${TRANSITIONS_PER_PRECURSOR@Q}),
    "runs": run_names,
    "study": {
        "control_runs": int(${CONTROL_RUNS@Q}),
        "treatment_runs": int(${TREATMENT_RUNS@Q}),
        "run_log2_sd": float(${RUN_LOG2_SD@Q}),
        "protein_log2_sd": float(${PROTEIN_LOG2_SD@Q}),
        "peptide_log2_sd": float(${PEPTIDE_LOG2_SD@Q}),
        "treatment_up_fraction": float(${TREATMENT_UP_FRACTION@Q}),
        "treatment_down_fraction": float(${TREATMENT_DOWN_FRACTION@Q}),
        "treatment_min_abs_log2fc": float(${TREATMENT_MIN_ABS_LOG2FC@Q}),
        "treatment_max_abs_log2fc": float(${TREATMENT_MAX_ABS_LOG2FC@Q}),
        "manifest_tsv": "OpenSwathTimSim.study_manifest.tsv",
        "design_truth_tsv": "OpenSwathTimSim.study_design_truth.tsv",
        "abundance_truth_tsv": "OpenSwathTimSim.study_abundance_truth.tsv",
        "realized_truth_tsv": "OpenSwathTimSim.realized_truth.tsv",
    },
    "randomization": {
        "seed_plan_version": 3,
        "blueprint_seed": int(${BLUEPRINT_SEED@Q}),
        "fasta_seed": int(${FASTA_SEED@Q}),
        "study_seed": int(${STUDY_SEED@Q}),
        "sample_seed_base": int(${SAMPLE_SEED_BASE@Q}),
        "entrapment_seed": int(${ENTRAPMENT_SEED@Q}),
    },
    "selection": {
        "scope": "blueprint_only_before_condition_effects",
        "mode": ${SELECTION_MODE@Q},
        "target_proteins": int(${TARGET_PROTEINS@Q}),
        "precursors_per_protein": int(${PRECURSORS_PER_PROTEIN@Q}),
        "collision_specificity": {
            "source": "complete_TimSim_blueprint_fragment_geometry",
            "uses_opendia": False,
            "uses_observed_raw_intensity": False,
            "max_rt_seconds": 6.0,
            "max_im_1_over_k0": 0.03,
            "max_fragment_ppm": 25.0,
            "top_fragments": 4,
            "signal_shortlist_multiplier": 2,
        },
        "min_realized_event_proxy_blueprint": float(${MIN_REALIZED_EVENT_PROXY@Q}),
        "min_frame_abundance_sum_blueprint": float(${MIN_FRAME_ABUNDANCE_SUM@Q}),
        "min_scan_abundance_sum_blueprint": float(${MIN_SCAN_ABUNDANCE_SUM@Q}),
        "min_ion_relative_abundance_blueprint": float(${MIN_ION_RELATIVE_ABUNDANCE@Q}),
        "selected_precursors_tsv": "OpenSwathTimSim.high_signal_selection.tsv",
        "reference_truth_tsv": "OpenSwathTimSim.reference_truth.tsv",
        "mz_window_margin_th": float(${MZ_WINDOW_MARGIN@Q}),
        "im_window_margin_1_over_k0": float(${IM_WINDOW_MARGIN@Q}),
    },
    "entrapment": {
        "enabled": True,
        "precursors": int(${ENTRAPMENTS@Q}),
        "seed": int(${ENTRAPMENT_SEED@Q}),
        "active_mode": "independent",
        "fragment_collision_ppm": 25.0,
        "max_paired_fragment_collisions": 1,
        "independent_min_rt_separation": 2.0,
        "independent_min_im_separation": 0.03,
        "role": "canonical_external_null",
        "calibration_eligible": True,
        "target_only_transition_tsv": "OpenSwathTimSim.target_only.transitions.tsv",
        "combined_transition_tsv": "OpenSwathTimSim.transitions.tsv",
        "truth_tsv": "OpenSwathTimSim.entrapment_truth.tsv",
    },
}
artifacts = {
    "synthetic_fasta": root / "generated_inputs/synthetic_proteome.fasta",
    "study_manifest": root / "OpenSwathTimSim.study_manifest.tsv",
    "study_design_truth": root / "OpenSwathTimSim.study_design_truth.tsv",
    "study_abundance_truth": root / "OpenSwathTimSim.study_abundance_truth.tsv",
    "selected_precursors": root / "OpenSwathTimSim.high_signal_selection.tsv",
    "reference_truth": root / "OpenSwathTimSim.reference_truth.tsv",
    "realized_truth": root / "OpenSwathTimSim.realized_truth.tsv",
    "target_only_transitions": root / "OpenSwathTimSim.target_only.transitions.tsv",
    "combined_transitions": root / "OpenSwathTimSim.transitions.tsv",
    "entrapment_truth": root / "OpenSwathTimSim.entrapment_truth.tsv",
}
manifest["artifact_sha256"] = {name: sha256_file(file_path) for name, file_path in artifacts.items()}
path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(f"Wrote {path}")
PY

python "$ROOT/tools/validate_study.py" "$OUTPUT_ROOT"
mapfile -t RUN_DBS < <(printf '%s\n' "${RUN_NAMES[@]}" | sed "s#^#$OUTPUT_ROOT/#; s#\$#/synthetic_data.db#")
python "$ROOT/tools/validate_experiment.py" "${RUN_DBS[@]}"

printf '\nMulti-run study generated under: %s\n' "$OUTPUT_ROOT"
printf 'Runs: %d control + %d treatment = %d\n' "$CONTROL_RUNS" "$TREATMENT_RUNS" "$EXPECTED_RUNS"
printf 'Targets: %d; independent entrapments: %d\n' "$PRECURSORS" "$ENTRAPMENTS"
if [[ "$SELECTION_MODE" == "protein_balanced" ]]; then
  printf 'Production target composition: %d proteins x %d precursors/protein\n' "$TARGET_PROTEINS" "$PRECURSORS_PER_PROTEIN"
fi
printf 'Study manifest: %s\n' "$STUDY_MANIFEST"
printf 'Design truth:   %s\n' "$DESIGN_TRUTH"
printf 'Realized truth: %s\n' "$REALIZED_TRUTH"
