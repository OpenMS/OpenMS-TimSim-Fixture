#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"
REFERENCE_D=""
OUTPUT_ROOT="$DEFAULT_FIXTURE_BUILD"
PRECURSORS="${MAX_PRECURSORS:-500}"
ENTRAPMENTS="${ENTRAPMENT_PRECURSORS:-}"
ENTRAPMENT_SEED="${ENTRAPMENT_SEED:-20260812}"
ENTRAPMENT_MODE="${ENTRAPMENT_MODE:-independent}"
ENTRAPMENT_INDEPENDENT_MIN_RT_SEPARATION="${ENTRAPMENT_INDEPENDENT_MIN_RT_SEPARATION:-2.0}"
ENTRAPMENT_INDEPENDENT_MIN_IM_SEPARATION="${ENTRAPMENT_INDEPENDENT_MIN_IM_SEPARATION:-0.03}"
ENTRAPMENT_FRAGMENT_COLLISION_PPM="${ENTRAPMENT_FRAGMENT_COLLISION_PPM:-25}"
ENTRAPMENT_MAX_PAIRED_FRAGMENT_COLLISIONS="${ENTRAPMENT_MAX_PAIRED_FRAGMENT_COLLISIONS:-1}"
SIMULATED_PEPTIDES=""
FASTA_PEPTIDES=""
PEPTIDES_PER_PROTEIN=20
GRADIENT_LENGTH=30
TRANSITIONS_PER_PRECURSOR=8
MINIMUM_FRAGMENTS=8
# Reproducibility seed plan. The FASTA seed is intentionally a separate fixed
# study invariant; SIMULATION_SEED controls stochastic realization of one DIA
# experiment. Child offsets are explicit and recorded in fixture_manifest.json.
SIMULATION_SEED="${SIMULATION_SEED:-2026081301}"
FASTA_SEED="${FASTA_SEED:-1729}"
SAMPLE_SEED=""
CONDITION_SEED=""
TIMSIM_THREADS=1
BATCH_SIZE=128
FRAME_BATCH_SIZE=100
VENV="$FIXTURE_VENV"

MIN_REALIZED_EVENT_PROXY=50000
MIN_FRAME_ABUNDANCE_SUM=0.90
MIN_SCAN_ABUNDANCE_SUM=0.95
MIN_ION_RELATIVE_ABUNDANCE=0.25
RT_EDGE_MARGIN_FRACTION=0.05
RT_EDGE_MIN_SECONDS=6.0
MZ_WINDOW_MARGIN=1.0
IM_WINDOW_MARGIN=0.005
RT_BINS=5
MZ_BINS=5
IM_BINS=5

REPLICATE_LOG2_SD=0.10
TREATMENT_LOG2_SD=0.08
TREATMENT_UP_FRACTION=0.15
TREATMENT_DOWN_FRACTION=0.15
TREATMENT_UP_FOLD=2.0
TREATMENT_DOWN_FOLD=0.5

usage() {
  cat >&2 <<EOF
Usage:
  $0 --reference /path/to/reference.d [options]

Build a reproducible TimSim/OpenDIA synthetic DIA fixture. Target selection
is based only on TimSim truth and acquisition/library geometry; OpenDIA results
are not consulted during selection.

Required:
  --reference PATH                 Bruker DIA-PASEF reference .d directory

Main options:
  --output PATH                    Output directory (default: $DEFAULT_FIXTURE_BUILD; must be outside source tree)
  --precursors N                   Exact selected true target groups (default: 500)
  --entrapments N                  Known-absent matched entrapment groups (default: same as --precursors)
  --entrapment-seed N              Deterministic entrapment shuffle seed (default: 20260812)
  --entrapment-mode MODE           independent (canonical FDR null) or paired_hard (adversarial stress)
                                   (default: independent)
  --simulation-seed N                Outer TimSim/DIA realization seed
                                   (default: 2026081301; may also use SIMULATION_SEED)
  --fasta-seed N                     Synthetic-proteome universe seed (default: 1729;
                                   intentionally fixed across outer simulations)
  --simulated-peptides N           Candidate peptides requested from TimSim
                                   (default: max(8 x precursors, precursors + 2500))
  --fasta-peptides N               Designed FASTA peptide pool
                                   (default: max(2 x simulated peptides, 4000))
  --gradient-length SECONDS        Simulated LC gradient (default: 30)
  --timsim-threads N               TimSim worker threads (default: 1)

Correctness-selection thresholds:
  --min-realized-event-proxy X     Minimum in every run (default: 50000)
  --min-frame-abundance-sum X      Minimum retained RT mass in every run (default: 0.90)
  --min-scan-abundance-sum X       Minimum retained IM mass in every run (default: 0.95)
  --min-ion-relative-abundance X   Minimum selected charge fraction in every run (default: 0.25)
  --rt-edge-margin-fraction X      Fraction of acquisition duration excluded at each RT edge (default: 0.05)
  --rt-edge-min-seconds X          Absolute RT-edge margin floor in seconds (default: 6.0)
                                   Effective margin = max(fractional margin, this floor)
  --mz-window-margin X             Safe DIA m/z margin in Th (default: 1.0)
  --im-window-margin X             Safe DIA IM margin in 1/K0 (default: 0.005)

Condition-effect options:
  --replicate-log2-sd X            Run-2 peptide variation (default: 0.10)
  --treatment-log2-sd X            Run-3 residual sample variation (default: 0.08)
  --treatment-up-fraction X        Fraction of proteins increased (default: 0.15)
  --treatment-down-fraction X      Fraction of proteins decreased (default: 0.15)
  --treatment-up-fold X            Increased-protein fold change (default: 2.0)
  --treatment-down-fold X          Decreased-protein fold change (default: 0.5)

If fewer than N candidates satisfy the fixed correctness criteria, generation
fails after writing selection QC. Increase --simulated-peptides and --fasta-peptides;
do not weaken the thresholds simply to reach N.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference) REFERENCE_D="${2:?missing value for --reference}"; shift 2 ;;
    --output) OUTPUT_ROOT="${2:?missing value for --output}"; shift 2 ;;
    --precursors) PRECURSORS="${2:?missing value for --precursors}"; shift 2 ;;
    --entrapments) ENTRAPMENTS="${2:?missing value for --entrapments}"; shift 2 ;;
    --entrapment-seed) ENTRAPMENT_SEED="${2:?missing value for --entrapment-seed}"; shift 2 ;;
    --entrapment-mode) ENTRAPMENT_MODE="${2:?missing value for --entrapment-mode}"; shift 2 ;;
    --simulation-seed) SIMULATION_SEED="${2:?missing value for --simulation-seed}"; shift 2 ;;
    --fasta-seed) FASTA_SEED="${2:?missing value for --fasta-seed}"; shift 2 ;;
    --simulated-peptides) SIMULATED_PEPTIDES="${2:?missing value for --simulated-peptides}"; shift 2 ;;
    --fasta-peptides) FASTA_PEPTIDES="${2:?missing value for --fasta-peptides}"; shift 2 ;;
    --gradient-length) GRADIENT_LENGTH="${2:?missing value for --gradient-length}"; shift 2 ;;
    --timsim-threads) TIMSIM_THREADS="${2:?missing value for --timsim-threads}"; shift 2 ;;
    --min-realized-event-proxy) MIN_REALIZED_EVENT_PROXY="${2:?missing value}"; shift 2 ;;
    --min-frame-abundance-sum) MIN_FRAME_ABUNDANCE_SUM="${2:?missing value}"; shift 2 ;;
    --min-scan-abundance-sum) MIN_SCAN_ABUNDANCE_SUM="${2:?missing value}"; shift 2 ;;
    --min-ion-relative-abundance) MIN_ION_RELATIVE_ABUNDANCE="${2:?missing value}"; shift 2 ;;
    --rt-edge-margin-fraction) RT_EDGE_MARGIN_FRACTION="${2:?missing value}"; shift 2 ;;
    --rt-edge-min-seconds) RT_EDGE_MIN_SECONDS="${2:?missing value}"; shift 2 ;;
    --mz-window-margin) MZ_WINDOW_MARGIN="${2:?missing value}"; shift 2 ;;
    --im-window-margin) IM_WINDOW_MARGIN="${2:?missing value}"; shift 2 ;;
    --replicate-log2-sd) REPLICATE_LOG2_SD="${2:?missing value}"; shift 2 ;;
    --treatment-log2-sd) TREATMENT_LOG2_SD="${2:?missing value}"; shift 2 ;;
    --treatment-up-fraction) TREATMENT_UP_FRACTION="${2:?missing value}"; shift 2 ;;
    --treatment-down-fraction) TREATMENT_DOWN_FRACTION="${2:?missing value}"; shift 2 ;;
    --treatment-up-fold) TREATMENT_UP_FOLD="${2:?missing value}"; shift 2 ;;
    --treatment-down-fold) TREATMENT_DOWN_FOLD="${2:?missing value}"; shift 2 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "ERROR: Unknown option: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$REFERENCE_D" ]]; then
  usage
  exit 2
fi
if [[ ! "$PRECURSORS" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --precursors must be a positive integer" >&2
  exit 2
fi
if [[ ! "$SIMULATION_SEED" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --simulation-seed must be a non-negative integer" >&2
  exit 2
fi
if [[ ! "$FASTA_SEED" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --fasta-seed must be a non-negative integer" >&2
  exit 2
fi
# Keep TimSim's own sample seed and the biological-condition perturbation seed
# explicit and reproducible. TimSim/rustims determinism fixes are provided by
# upstream imspy-simulation 0.4.2; no local patching or process-level RNG shim is used.
if (( SIMULATION_SEED > 4294967275 )); then
  echo "ERROR: --simulation-seed must be <= 4294967275" >&2
  exit 2
fi
if (( FASTA_SEED > 4294967295 )); then
  echo "ERROR: --fasta-seed must be <= 4294967295" >&2
  exit 2
fi

SAMPLE_SEED=$(( SIMULATION_SEED + 10 ))
CONDITION_SEED=$(( SIMULATION_SEED + 20 ))

printf 'Seed plan: simulation=%s, TimSim sample=%s, condition=%s, FASTA=%s\n' \
  "$SIMULATION_SEED" "$SAMPLE_SEED" "$CONDITION_SEED" "$FASTA_SEED"

run_timsim() {
  local config="$1"
  if "$VENV/bin/timsim" --help 2>&1 | grep -q -- '--config'; then
    "$VENV/bin/timsim" --config "$config"
  else
    "$VENV/bin/timsim" "$config"
  fi
}

if [[ -z "$ENTRAPMENTS" ]]; then
  ENTRAPMENTS="$PRECURSORS"
fi
if [[ ! "$ENTRAPMENTS" =~ ^[0-9]+$ ]]; then
  echo "ERROR: --entrapments must be a non-negative integer" >&2
  exit 2
fi
if [[ -z "$SIMULATED_PEPTIDES" ]]; then
  SIMULATED_PEPTIDES=$(( PRECURSORS * 8 ))
  if (( SIMULATED_PEPTIDES < PRECURSORS + 2500 )); then
    SIMULATED_PEPTIDES=$(( PRECURSORS + 2500 ))
  fi
fi
if [[ ! "$SIMULATED_PEPTIDES" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --simulated-peptides must be a positive integer" >&2
  exit 2
fi
if (( SIMULATED_PEPTIDES < PRECURSORS )); then
  echo "ERROR: --simulated-peptides must be >= --precursors" >&2
  exit 2
fi
if [[ -z "$FASTA_PEPTIDES" ]]; then
  FASTA_PEPTIDES=$(( SIMULATED_PEPTIDES * 2 ))
  if (( FASTA_PEPTIDES < 4000 )); then
    FASTA_PEPTIDES=4000
  fi
fi
if [[ ! "$FASTA_PEPTIDES" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --fasta-peptides must be a positive integer" >&2
  exit 2
fi
if (( FASTA_PEPTIDES < SIMULATED_PEPTIDES )); then
  echo "ERROR: --fasta-peptides must be >= --simulated-peptides" >&2
  exit 2
fi

if [[ ! -f "$VENV/bin/activate" ]]; then
  echo "ERROR: TimSim environment not found at $VENV" >&2
  echo "Run: ./scripts/setup.sh" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"

REFERENCE_D="$(realpath "$REFERENCE_D")"
require_external_output_path "$OUTPUT_ROOT" "--output"
if [[ ! -f "$REFERENCE_D/analysis.tdf" || ! -f "$REFERENCE_D/analysis.tdf_bin" ]]; then
  echo "ERROR: Reference must contain analysis.tdf and analysis.tdf_bin: $REFERENCE_D" >&2
  exit 1
fi
mkdir -p "$OUTPUT_ROOT"
OUTPUT_ROOT="$(realpath "$OUTPUT_ROOT")"

RUN1="OpenSwathTimSim_01_baseline"
RUN2="OpenSwathTimSim_02_biological_replicate"
RUN3="OpenSwathTimSim_03_treatment"
RUN_NAMES=("$RUN1" "$RUN2" "$RUN3")

GENERATED_INPUTS="$OUTPUT_ROOT/generated_inputs"
CONDITION_INPUTS="$OUTPUT_ROOT/condition_inputs"
RENDERED="$OUTPUT_ROOT/rendered_configs"
FASTA="$GENERATED_INPUTS/synthetic_proteome.fasta"
REPLICATE_INPUT_DIR="$CONDITION_INPUTS/control_replicate"
TREATMENT_INPUT_DIR="$CONDITION_INPUTS/treatment"
REPLICATE_DB="$REPLICATE_INPUT_DIR/synthetic_data.db"
TREATMENT_DB="$TREATMENT_INPUT_DIR/synthetic_data.db"
CONDITION_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.condition_truth.tsv"
SELECTED_TSV="$OUTPUT_ROOT/OpenSwathTimSim.high_signal_selection.tsv"
REALIZED_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.realized_truth.tsv"
CANDIDATE_QC="$OUTPUT_ROOT/OpenSwathTimSim.selection_candidates.tsv"
QC_DIR="$OUTPUT_ROOT/selection_qc"
QC_JSON="$QC_DIR/selection_summary.json"
QC_REPORT="$QC_DIR/selection_report.md"

for name in "${RUN_NAMES[@]}"; do
  rm -rf "$OUTPUT_ROOT/$name"
done
rm -rf "$OUTPUT_ROOT/tdfs" "$GENERATED_INPUTS" "$CONDITION_INPUTS" "$RENDERED" "$QC_DIR"
rm -f \
  "$OUTPUT_ROOT/OpenSwathTimSim.transitions.tsv" \
  "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" \
  "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_truth.tsv" \
  "$OUTPUT_ROOT/OpenSwathTimSim.ground_truth.tsv" \
  "$OUTPUT_ROOT/OpenSwathTimSim.condition_truth.tsv" \
  "$SELECTED_TSV" "$REALIZED_TRUTH" "$CANDIDATE_QC" \
  "$OUTPUT_ROOT/fixture_manifest.json"
mkdir -p "$GENERATED_INPUTS" "$REPLICATE_INPUT_DIR" "$TREATMENT_INPUT_DIR" "$RENDERED" "$QC_DIR"

python "$ROOT/tools/generate_synthetic_fasta.py" \
  --out "$FASTA" \
  --peptides "$FASTA_PEPTIDES" \
  --peptides-per-protein "$PEPTIDES_PER_PROTEIN" \
  --seed "$FASTA_SEED"

N_PROTEINS=$(( (FASTA_PEPTIDES + PEPTIDES_PER_PROTEIN - 1) / PEPTIDES_PER_PROTEIN ))

python "$ROOT/tools/render_configs.py" \
  --root "$ROOT" \
  --reference "$REFERENCE_D" \
  --output-root "$OUTPUT_ROOT" \
  --rendered-dir "$RENDERED" \
  --fasta "$FASTA" \
  --replicate-existing-path "$REPLICATE_INPUT_DIR" \
  --treatment-existing-path "$TREATMENT_INPUT_DIR" \
  --n-proteins "$N_PROTEINS" \
  --num-peptides-total "$FASTA_PEPTIDES" \
  --num-sample-peptides "$SIMULATED_PEPTIDES" \
  --gradient-length "$GRADIENT_LENGTH" \
  --sample-seed "$SAMPLE_SEED" \
  --timsim-threads "$TIMSIM_THREADS" \
  --batch-size "$BATCH_SIZE" \
  --frame-batch-size "$FRAME_BATCH_SIZE"

run_timsim "$RENDERED/run_01_baseline.toml"
BASE_DB="$OUTPUT_ROOT/$RUN1/synthetic_data.db"
[[ -f "$BASE_DB" ]] || { echo "ERROR: Expected baseline database was not found: $BASE_DB" >&2; exit 1; }

python "$ROOT/tools/prepare_condition_databases.py" \
  --baseline-db "$BASE_DB" \
  --replicate-db "$REPLICATE_DB" \
  --treatment-db "$TREATMENT_DB" \
  --truth-out "$CONDITION_TRUTH" \
  --seed "$CONDITION_SEED" \
  --replicate-log2-sd "$REPLICATE_LOG2_SD" \
  --treatment-log2-sd "$TREATMENT_LOG2_SD" \
  --treatment-up-fraction "$TREATMENT_UP_FRACTION" \
  --treatment-down-fraction "$TREATMENT_DOWN_FRACTION" \
  --treatment-up-fold "$TREATMENT_UP_FOLD" \
  --treatment-down-fold "$TREATMENT_DOWN_FOLD"

run_timsim "$RENDERED/run_02_biological_replicate.toml"
run_timsim "$RENDERED/run_03_treatment.toml"

RUN2_DB="$OUTPUT_ROOT/$RUN2/synthetic_data.db"
RUN3_DB="$OUTPUT_ROOT/$RUN3/synthetic_data.db"
for db in "$BASE_DB" "$RUN2_DB" "$RUN3_DB"; do
  [[ -f "$db" ]] || { echo "ERROR: Expected TimSim output DB not found: $db" >&2; exit 1; }
done

# Freeze the correctness set before any OpenDIA invocation. This script has no
# OpenDIA result inputs and fails if the simulator-only criteria yield < N groups.
python "$ROOT/tools/select_high_signal_precursors.py" \
  --baseline-db "$BASE_DB" \
  --control-db "$RUN2_DB" \
  --treatment-db "$RUN3_DB" \
  --selected-out "$SELECTED_TSV" \
  --truth-out "$REALIZED_TRUTH" \
  --candidate-qc-out "$CANDIDATE_QC" \
  --qc-json "$QC_JSON" \
  --qc-report "$QC_REPORT" \
  --precursors "$PRECURSORS" \
  --min-realized-event-proxy "$MIN_REALIZED_EVENT_PROXY" \
  --min-frame-abundance-sum "$MIN_FRAME_ABUNDANCE_SUM" \
  --min-scan-abundance-sum "$MIN_SCAN_ABUNDANCE_SUM" \
  --min-ion-relative-abundance "$MIN_ION_RELATIVE_ABUNDANCE" \
  --rt-edge-margin-fraction "$RT_EDGE_MARGIN_FRACTION" \
  --rt-edge-min-seconds "$RT_EDGE_MIN_SECONDS" \
  --minimum-fragments "$MINIMUM_FRAGMENTS" \
  --min-product-mz 350 \
  --max-product-mz 2000 \
  --mz-window-margin "$MZ_WINDOW_MARGIN" \
  --im-window-margin "$IM_WINDOW_MARGIN" \
  --rt-bins "$RT_BINS" \
  --mz-bins "$MZ_BINS" \
  --im-bins "$IM_BINS"

python "$ROOT/tools/export_openswath_tsv.py" \
  --db "$BASE_DB" \
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

case "$ENTRAPMENT_MODE" in
  paired_hard|independent) ;;
  *) echo "ERROR: --entrapment-mode must be paired_hard or independent" >&2; exit 2 ;;
esac

ENTRAPMENT_MODE_COMBINED="$OUTPUT_ROOT/OpenSwathTimSim.entrapment_${ENTRAPMENT_MODE}.transitions.tsv"
ENTRAPMENT_MODE_TRUTH="$OUTPUT_ROOT/OpenSwathTimSim.entrapment_${ENTRAPMENT_MODE}_truth.tsv"
python "$ROOT/tools/generate_entrapment_library.py" \
  --target-transitions "$OUTPUT_ROOT/OpenSwathTimSim.target_only.transitions.tsv" \
  --baseline-db "$BASE_DB" \
  --out "$ENTRAPMENT_MODE_COMBINED" \
  --truth-out "$ENTRAPMENT_MODE_TRUTH" \
  --mode "$ENTRAPMENT_MODE" \
  --count "$ENTRAPMENTS" \
  --seed "$ENTRAPMENT_SEED" \
  --fragment-collision-ppm "$ENTRAPMENT_FRAGMENT_COLLISION_PPM" \
  --max-paired-fragment-collisions "$ENTRAPMENT_MAX_PAIRED_FRAGMENT_COLLISIONS" \
  --independent-min-rt-separation "$ENTRAPMENT_INDEPENDENT_MIN_RT_SEPARATION" \
  --independent-min-im-separation "$ENTRAPMENT_INDEPENDENT_MIN_IM_SEPARATION" \
  --mz-window-margin "$MZ_WINDOW_MARGIN" \
  --im-window-margin "$IM_WINDOW_MARGIN"
cp "$ENTRAPMENT_MODE_COMBINED" "$OUTPUT_ROOT/OpenSwathTimSim.transitions.tsv"
cp "$ENTRAPMENT_MODE_TRUTH" "$OUTPUT_ROOT/OpenSwathTimSim.entrapment_truth.tsv"

mkdir -p "$OUTPUT_ROOT/tdfs"
for name in "${RUN_NAMES[@]}"; do
  source_d="$OUTPUT_ROOT/$name/$name.d"
  [[ -d "$source_d" ]] || { echo "ERROR: Expected generated TDF directory was not found: $source_d" >&2; exit 1; }
  ln -s "../$name/$name.d" "$OUTPUT_ROOT/tdfs/$name.d"
done

python - "$OUTPUT_ROOT/fixture_manifest.json" <<PY
import hashlib
import importlib.metadata
import json
import platform
import sys
from pathlib import Path

path = Path(sys.argv[1])
manifest = {
    "fixture_kind": "openms_timsim_identification_quantification_entrapment",
    "selection_frozen_before_opendia": True,
    "randomization": {
        "seed_plan_version": 2,
        "simulation_seed": int(${SIMULATION_SEED@Q}),
        "fasta_seed": int(${FASTA_SEED@Q}),
        "fasta_seed_role": "fixed_synthetic_proteome_universe",
        "timsim_sample_seed": int(${SAMPLE_SEED@Q}),
        "condition_seed": int(${CONDITION_SEED@Q}),
        "timsim_execution": "upstream_timsim_cli",
        "local_determinism_patches": False,
        "child_seed_offsets": {
            "timsim_sample_seed": 10,
            "condition_seed": 20,
        },
    },
    "software": {
        "python": platform.python_version(),
        "imspy-simulation": importlib.metadata.version("imspy-simulation"),
        "imspy-connector": importlib.metadata.version("imspy-connector"),
    },
    "reference_d": ${REFERENCE_D@Q},
    "requested_precursors": int(${PRECURSORS@Q}),
    "simulated_peptides": int(${SIMULATED_PEPTIDES@Q}),
    "fasta_peptides": int(${FASTA_PEPTIDES@Q}),
    "synthetic_proteins": int(${N_PROTEINS@Q}),
    "gradient_length_seconds": float(${GRADIENT_LENGTH@Q}),
    "transitions_per_precursor": int(${TRANSITIONS_PER_PRECURSOR@Q}),
    "minimum_fragments": int(${MINIMUM_FRAGMENTS@Q}),
    "entrapment": {
        "enabled": int(${ENTRAPMENTS@Q}) > 0,
        "precursors": int(${ENTRAPMENTS@Q}),
        "seed": int(${ENTRAPMENT_SEED@Q}),
        "active_mode": ${ENTRAPMENT_MODE@Q},
        "role": "canonical_external_null" if ${ENTRAPMENT_MODE@Q} == "independent" else "adversarial_interference_stress",
        "calibration_eligible": ${ENTRAPMENT_MODE@Q} == "independent",
        "fragment_collision_ppm": float(${ENTRAPMENT_FRAGMENT_COLLISION_PPM@Q}),
        "max_paired_fragment_collisions": int(${ENTRAPMENT_MAX_PAIRED_FRAGMENT_COLLISIONS@Q}),
        "independent_min_rt_separation": float(${ENTRAPMENT_INDEPENDENT_MIN_RT_SEPARATION@Q}),
        "independent_min_im_separation": float(${ENTRAPMENT_INDEPENDENT_MIN_IM_SEPARATION@Q}),
        "target_only_transition_tsv": "OpenSwathTimSim.target_only.transitions.tsv",
        "combined_transition_tsv": "OpenSwathTimSim.transitions.tsv",
        "truth_tsv": "OpenSwathTimSim.entrapment_truth.tsv",
        "mode_specific_transition_tsv": "OpenSwathTimSim.entrapment_" + ${ENTRAPMENT_MODE@Q} + ".transitions.tsv",
        "mode_specific_truth_tsv": "OpenSwathTimSim.entrapment_" + ${ENTRAPMENT_MODE@Q} + "_truth.tsv",
    },
    "runs": [${RUN1@Q}, ${RUN2@Q}, ${RUN3@Q}],
    "selection": {
        "min_realized_event_proxy_all_runs": float(${MIN_REALIZED_EVENT_PROXY@Q}),
        "min_frame_abundance_sum_all_runs": float(${MIN_FRAME_ABUNDANCE_SUM@Q}),
        "min_scan_abundance_sum_all_runs": float(${MIN_SCAN_ABUNDANCE_SUM@Q}),
        "min_ion_relative_abundance_all_runs": float(${MIN_ION_RELATIVE_ABUNDANCE@Q}),
        "rt_edge_margin_fraction": float(${RT_EDGE_MARGIN_FRACTION@Q}),
        "rt_edge_min_seconds": float(${RT_EDGE_MIN_SECONDS@Q}),
        "rt_edge_policy": "max(fraction_of_acquisition, absolute_seconds)",
        "mz_window_margin_th": float(${MZ_WINDOW_MARGIN@Q}),
        "im_window_margin_1_over_k0": float(${IM_WINDOW_MARGIN@Q}),
        "selected_precursors_tsv": "OpenSwathTimSim.high_signal_selection.tsv",
        "realized_truth_tsv": "OpenSwathTimSim.realized_truth.tsv",
        "candidate_qc_tsv": "OpenSwathTimSim.selection_candidates.tsv",
        "qc_summary_json": "selection_qc/selection_summary.json",
        "qc_report_md": "selection_qc/selection_report.md",
    },
    "conditions": {
        ${RUN1@Q}: "baseline/control",
        ${RUN2@Q}: "control replicate with peptide-level abundance variation",
        ${RUN3@Q}: "treatment with coherent protein-level fold changes",
    },
    "replicate_log2_sd": float(${REPLICATE_LOG2_SD@Q}),
    "treatment_log2_sd": float(${TREATMENT_LOG2_SD@Q}),
    "treatment_up_fraction": float(${TREATMENT_UP_FRACTION@Q}),
    "treatment_down_fraction": float(${TREATMENT_DOWN_FRACTION@Q}),
    "treatment_up_fold": float(${TREATMENT_UP_FOLD@Q}),
    "treatment_down_fold": float(${TREATMENT_DOWN_FOLD@Q}),
}
def sha256_file(file_path: Path) -> str:
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

artifact_paths = {
    "synthetic_fasta": Path(${FASTA@Q}),
    "condition_truth": Path(${CONDITION_TRUTH@Q}),
    "selected_precursors": Path(${SELECTED_TSV@Q}),
    "realized_truth": Path(${REALIZED_TRUTH@Q}),
    "target_only_transitions": Path(${OUTPUT_ROOT@Q}) / "OpenSwathTimSim.target_only.transitions.tsv",
    "combined_transitions": Path(${OUTPUT_ROOT@Q}) / "OpenSwathTimSim.transitions.tsv",
    "entrapment_truth": Path(${OUTPUT_ROOT@Q}) / "OpenSwathTimSim.entrapment_truth.tsv",
    "ground_truth": Path(${OUTPUT_ROOT@Q}) / "OpenSwathTimSim.ground_truth.tsv",
}
manifest["artifact_sha256"] = {
    name: sha256_file(file_path) for name, file_path in artifact_paths.items()
}

path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(f"Wrote {path}")
PY

"$ROOT/scripts/check_fixture.sh" "$OUTPUT_ROOT"

printf '\nOpenMS-TimSim fixture generated under: %s\n' "$OUTPUT_ROOT"
printf 'Frozen selection: %s\n' "$SELECTED_TSV"
printf 'Realized truth:   %s\n' "$REALIZED_TRUTH"
printf 'Selection QC:     %s\n' "$QC_REPORT"
printf '\nOpenDIA input symlinks:\n'
for link in "$OUTPUT_ROOT"/tdfs/*.d; do
  printf '  %s -> %s\n' "$link" "$(readlink "$link")"
done
