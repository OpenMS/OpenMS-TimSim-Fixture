#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

CONFIG="${1:?usage: run_timsim_config.sh CONFIG.toml}"
[[ -s "$CONFIG" ]] || { echo "ERROR: config does not exist: $CONFIG" >&2; exit 2; }

VENV="${FIXTURE_VENV:-/opt/venv}"
TIMSIM="$VENV/bin/timsim"
PYTHON="$VENV/bin/python"
[[ -x "$TIMSIM" && -x "$PYTHON" ]] || { echo "ERROR: TimSim runtime not found under $VENV" >&2; exit 1; }

if grep -Eq '^use_gpu[[:space:]]*=[[:space:]]*true[[:space:]]*$' "$CONFIG"; then
  "$PYTHON" - <<'PY'
import torch
if not torch.cuda.is_available():
    raise SystemExit("ERROR: config requests use_gpu=true but CUDA is unavailable; use singularity exec --nv")
print(f"CUDA available: {torch.cuda.get_device_name(0)}; torch={torch.__version__}; torch CUDA={torch.version.cuda}")
PY
fi

if "$TIMSIM" --help 2>&1 | grep -q -- '--config'; then
  exec "$TIMSIM" --config "$CONFIG"
else
  exec "$TIMSIM" "$CONFIG"
fi
