#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

FIXTURE="${1:-}"
MODE="${2:-independent}"
if [[ -z "$FIXTURE" || "$FIXTURE" == "--help" || "$FIXTURE" == "-h" ]]; then
  cat >&2 <<EOF
Usage: $0 FIXTURE_DIR [independent|paired_hard]

Regenerate only the target+entrapment library for an existing frozen fixture.
TimSim runs, selected true targets, and quantitative truth are not modified.

  independent  canonical external null for q-value/FDR calibration (default)
  paired_hard  adversarial interference stress test; not the canonical FDR null
EOF
  [[ -n "$FIXTURE" ]] && exit 0 || exit 2
fi
case "$MODE" in
  independent|paired_hard) ;;
  *) echo "ERROR: mode must be independent or paired_hard" >&2; exit 2 ;;
esac

FIXTURE="$(realpath "$FIXTURE")"
MANIFEST="$FIXTURE/fixture_manifest.json"
TARGETS="$FIXTURE/OpenSwathTimSim.target_only.transitions.tsv"
[[ -f "$MANIFEST" ]] || { echo "ERROR: missing fixture manifest: $MANIFEST" >&2; exit 1; }
[[ -f "$TARGETS" ]] || { echo "ERROR: missing target-only transitions: $TARGETS" >&2; exit 1; }
[[ -x "$FIXTURE_VENV/bin/python" ]] || { echo "ERROR: fixture environment not found: $FIXTURE_VENV" >&2; exit 1; }

manifest_value() {
  "$FIXTURE_VENV/bin/python" - "$MANIFEST" "$1" "$2" <<'PY'
import json, sys
m=json.load(open(sys.argv[1]))
path=sys.argv[2].split('.')
value=m
for key in path:
    value=value.get(key) if isinstance(value, dict) else None
    if value is None:
        break
print(sys.argv[3] if value is None else value)
PY
}

FIXTURE_KIND="$(manifest_value fixture_kind unknown)"
BLUEPRINT_RUN="$(manifest_value blueprint_run '')"
if [[ -n "$BLUEPRINT_RUN" ]]; then
  BASE_DB="$FIXTURE/$BLUEPRINT_RUN/synthetic_data.db"
else
  BASE_DB="$FIXTURE/OpenSwathTimSim_01_baseline/synthetic_data.db"
fi
[[ -f "$BASE_DB" ]] || { echo "ERROR: missing TimSim blueprint/baseline DB: $BASE_DB" >&2; exit 1; }

COUNT="$(manifest_value entrapment.precursors 0)"
SEED="$(manifest_value entrapment.seed 20260812)"
COLLISION_PPM="$(manifest_value entrapment.fragment_collision_ppm 25)"
MAX_COLLISIONS="$(manifest_value entrapment.max_paired_fragment_collisions 1)"
MIN_RT_SEP="$(manifest_value entrapment.independent_min_rt_separation 2.0)"
MIN_IM_SEP="$(manifest_value entrapment.independent_min_im_separation 0.03)"
MZ_MARGIN="$(manifest_value selection.mz_window_margin_th 1.0)"
IM_MARGIN="$(manifest_value selection.im_window_margin_1_over_k0 0.005)"

MODE_LIBRARY="$FIXTURE/OpenSwathTimSim.entrapment_${MODE}.transitions.tsv"
MODE_TRUTH="$FIXTURE/OpenSwathTimSim.entrapment_${MODE}_truth.tsv"

"$FIXTURE_VENV/bin/python" "$ROOT/tools/generate_entrapment_library.py" \
  --target-transitions "$TARGETS" \
  --baseline-db "$BASE_DB" \
  --out "$MODE_LIBRARY" \
  --truth-out "$MODE_TRUTH" \
  --mode "$MODE" \
  --count "$COUNT" \
  --seed "$SEED" \
  --fragment-collision-ppm "$COLLISION_PPM" \
  --max-paired-fragment-collisions "$MAX_COLLISIONS" \
  --independent-min-rt-separation "$MIN_RT_SEP" \
  --independent-min-im-separation "$MIN_IM_SEP" \
  --mz-window-margin "$MZ_MARGIN" \
  --im-window-margin "$IM_MARGIN"

cp "$MODE_LIBRARY" "$FIXTURE/OpenSwathTimSim.transitions.tsv"
cp "$MODE_TRUTH" "$FIXTURE/OpenSwathTimSim.entrapment_truth.tsv"

"$FIXTURE_VENV/bin/python" - "$MANIFEST" "$MODE" <<'PY'
import hashlib, json, sys
from pathlib import Path
path=Path(sys.argv[1]); mode=sys.argv[2]; root=path.parent
m=json.loads(path.read_text())
ent=m.setdefault('entrapment', {})
ent['active_mode']=mode
ent['role']='canonical_external_null' if mode == 'independent' else 'adversarial_interference_stress'
ent['calibration_eligible']=mode == 'independent'
ent['combined_transition_tsv']='OpenSwathTimSim.transitions.tsv'
ent['truth_tsv']='OpenSwathTimSim.entrapment_truth.tsv'
ent['mode_specific_transition_tsv']=f'OpenSwathTimSim.entrapment_{mode}.transitions.tsv'
ent['mode_specific_truth_tsv']=f'OpenSwathTimSim.entrapment_{mode}_truth.tsv'

def sha(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as fh:
        for block in iter(lambda: fh.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()
art=m.setdefault('artifact_sha256', {})
art['combined_transitions']=sha(root/'OpenSwathTimSim.transitions.tsv')
art['entrapment_truth']=sha(root/'OpenSwathTimSim.entrapment_truth.tsv')
path.write_text(json.dumps(m, indent=2)+'\n')
print(f'Updated {path}: active entrapment mode={mode}')
PY

if [[ "$FIXTURE_KIND" == "openms_timsim_multirun_study_v1" ]]; then
  "$FIXTURE_VENV/bin/python" "$ROOT/tools/validate_study.py" "$FIXTURE"
else
  "$ROOT/scripts/check_fixture.sh" "$FIXTURE"
fi
printf '\nActive entrapment mode: %s\n' "$MODE"
printf 'Library: %s\n' "$FIXTURE/OpenSwathTimSim.transitions.tsv"
printf 'Truth:   %s\n' "$FIXTURE/OpenSwathTimSim.entrapment_truth.tsv"
