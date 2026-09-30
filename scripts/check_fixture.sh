#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-}"
if [[ -z "$ROOT" ]]; then
  echo "Usage: $0 /path/to/generated/fixture" >&2
  exit 2
fi
ROOT="$(realpath "$ROOT")"
MANIFEST="$ROOT/fixture_manifest.json"
[[ -s "$MANIFEST" ]] || { echo "Missing $MANIFEST" >&2; exit 1; }

mapfile -t SETTINGS < <(python3 - "$MANIFEST" <<'PY'
import json
import sys
manifest = json.load(open(sys.argv[1], encoding="utf-8"))
print(manifest["requested_precursors"])
print(manifest["transitions_per_precursor"])
print(int(manifest.get("entrapment", {}).get("precursors", 0)))
for run in manifest["runs"]:
    print(run)
PY
)
EXPECTED_PRECURSORS="${SETTINGS[0]}"
EXPECTED_TRANSITIONS="${SETTINGS[1]}"
EXPECTED_ENTRAPMENTS="${SETTINGS[2]}"
RUN_NAMES=("${SETTINGS[3]}" "${SETTINGS[4]}" "${SETTINGS[5]}")

DBS=()
for name in "${RUN_NAMES[@]}"; do
  d="$ROOT/$name/$name.d"
  db="$ROOT/$name/synthetic_data.db"
  link="$ROOT/tdfs/$name.d"
  [[ -s "$d/analysis.tdf" ]] || { echo "Missing $d/analysis.tdf" >&2; exit 1; }
  [[ -s "$d/analysis.tdf_bin" ]] || { echo "Missing $d/analysis.tdf_bin" >&2; exit 1; }
  [[ -s "$db" ]] || { echo "Missing $db" >&2; exit 1; }
  [[ -L "$link" ]] || { echo "Missing TDF symlink $link" >&2; exit 1; }
  [[ "$(realpath "$link")" == "$(realpath "$d")" ]] || {
    echo "TDF symlink does not point to current generated run: $link" >&2
    exit 1
  }
  DBS+=("$db")
done

for required in \
  "$ROOT/OpenSwathTimSim.condition_truth.tsv" \
  "$ROOT/OpenSwathTimSim.high_signal_selection.tsv" \
  "$ROOT/OpenSwathTimSim.realized_truth.tsv" \
  "$ROOT/OpenSwathTimSim.selection_candidates.tsv" \
  "$ROOT/selection_qc/selection_summary.json" \
  "$ROOT/selection_qc/selection_report.md" \
  "$ROOT/OpenSwathTimSim.target_only.transitions.tsv"; do
  [[ -s "$required" ]] || { echo "Missing $required" >&2; exit 1; }
done
if (( EXPECTED_ENTRAPMENTS > 0 )); then
  [[ -s "$ROOT/OpenSwathTimSim.entrapment_truth.tsv" ]] || { echo "Missing entrapment truth TSV" >&2; exit 1; }
fi

TARGET_ONLY="$ROOT/OpenSwathTimSim.target_only.transitions.tsv"
COMBINED="$ROOT/OpenSwathTimSim.transitions.tsv"
ENTRAPMENT_TRUTH="$ROOT/OpenSwathTimSim.entrapment_truth.tsv"
BASE_DB="$ROOT/${RUN_NAMES[0]}/synthetic_data.db"
[[ -s "$COMBINED" ]] || { echo "Missing $COMBINED" >&2; exit 1; }

python3 - "$TARGET_ONLY" "$COMBINED" "$ENTRAPMENT_TRUTH" "$BASE_DB" \
  "$EXPECTED_PRECURSORS" "$EXPECTED_ENTRAPMENTS" "$EXPECTED_TRANSITIONS" <<'PY'
import csv
import sqlite3
import sys
from collections import Counter
from pathlib import Path

target_path, combined_path, ent_truth_path, db_path = map(Path, sys.argv[1:5])
expected_targets = int(sys.argv[5])
expected_entrapments = int(sys.argv[6])
expected_transitions = int(sys.argv[7])

def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

target = read(target_path)
combined = read(combined_path)
if not target or not combined:
    raise SystemExit("Transition TSV has no rows")

target_groups = Counter(row["TransitionGroupId"] for row in target)
if len(target_groups) != expected_targets:
    raise SystemExit(f"Target-only library has {len(target_groups)} groups; expected {expected_targets}")
wrong = {group: count for group, count in target_groups.items() if count != expected_transitions}
if wrong:
    raise SystemExit(f"Target-only transition counts are wrong; examples: {list(wrong.items())[:10]}")

combined_groups = Counter(row["TransitionGroupId"] for row in combined)
expected_groups = expected_targets + expected_entrapments
if len(combined_groups) != expected_groups:
    raise SystemExit(f"Combined library has {len(combined_groups)} groups; expected {expected_groups}")
wrong = {group: count for group, count in combined_groups.items() if count != expected_transitions}
if wrong:
    raise SystemExit(f"Combined transition counts are wrong; examples: {list(wrong.items())[:10]}")
if any(float(row["LibraryIntensity"]) <= 0 for row in combined):
    raise SystemExit("All library intensities must be positive")
if len({row["TransitionId"] for row in combined}) != len(combined):
    raise SystemExit("TransitionId values are not unique")

label_sequences = {}
for row in combined:
    label = str(row.get("PeptideGroupLabel", "")).strip()
    sequence = str(row.get("PeptideSequence", "")).strip()
    if not label:
        raise SystemExit("PeptideGroupLabel must be populated for every transition")
    label_sequences.setdefault(label, set()).add(sequence)
bad_labels = {label: sorted(sequences) for label, sequences in label_sequences.items() if len(sequences) != 1}
if bad_labels:
    raise SystemExit(
        "Each PeptideGroupLabel must map to exactly one PeptideSequence; examples: "
        f"{list(bad_labels.items())[:10]}"
    )

if any(str(row.get("Decoy", "0")).strip().lower() not in {"0", "false"} for row in combined):
    raise SystemExit("Input target/entrapment library must contain Decoy=0 only; OpenDIA appends decoys downstream")

if expected_entrapments:
    ent_truth = read(ent_truth_path)
    if len(ent_truth) != expected_entrapments:
        raise SystemExit(f"Entrapment truth has {len(ent_truth)} rows; expected {expected_entrapments}")
    ent_groups = {row["EntrapmentTransitionGroupId"] for row in ent_truth}
    ent_sequences = {row["EntrapmentSequence"].strip().upper() for row in ent_truth}
    if len(ent_groups) != expected_entrapments or len(ent_sequences) != expected_entrapments:
        raise SystemExit("Entrapment groups/sequences are not unique")
    observed_ent_groups = {group for group in combined_groups if group.startswith("ENTRAPMENT_")}
    if observed_ent_groups != ent_groups:
        raise SystemExit("Combined entrapment groups do not match OpenSwathTimSim.entrapment_truth.tsv")

    if any(not str(row.get("EntrapmentMode", "")).strip() for row in ent_truth):
        raise SystemExit("Entrapment truth is missing EntrapmentMode")
    modes = {str(row["EntrapmentMode"]).strip() for row in ent_truth}
    if len(modes) != 1:
        raise SystemExit(f"Entrapment truth mixes construction modes: {sorted(modes)}")
    mode = next(iter(modes))
    if mode == "independent":
        bad_same = [row for row in ent_truth if row.get("CoordinateDonorTransitionGroupId") == row.get("SourceTargetTransitionGroupId")]
        if bad_same:
            raise SystemExit("Independent entrapments contain source==coordinate-donor rows")
        bad_geometry = [row for row in ent_truth if str(row.get("IndependentCoordinatePasefCompatible", "")).lower() not in {"true", "1", "yes"}]
        if bad_geometry:
            raise SystemExit("Independent entrapments contain m/z/IM coordinates outside safe diaPASEF geometry")

    with sqlite3.connect(db_path) as con:
        cols = [row[1] for row in con.execute('PRAGMA table_info("peptides")')]
        seq_col = next((c for c in ("sequence", "peptide") if c in cols), None)
        if seq_col is None:
            raise SystemExit("Could not locate peptide sequence column in baseline DB")
        simulated = {str(row[0]).strip().upper() for row in con.execute(f'SELECT "{seq_col}" FROM peptides')}
    overlap = ent_sequences & simulated
    if overlap:
        raise SystemExit(f"Entrapment sequences occur in TimSim data; examples: {sorted(overlap)[:10]}")

print(
    f"Transition library OK: {expected_targets} true targets + {expected_entrapments} entrapments, "
    f"each with {expected_transitions} input transitions = {len(combined)} rows"
)
PY

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
python3 "$REPO_ROOT/tools/validate_experiment.py" "${DBS[@]}"
python3 "$REPO_ROOT/tools/validate_high_signal_selection.py" "$ROOT"
