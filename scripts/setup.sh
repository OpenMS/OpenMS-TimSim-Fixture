#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV="${FIXTURE_VENV:-$ROOT/.venv}"
PYTHON_VERSION="${PYTHON_VERSION:-3.11}"

command -v uv >/dev/null 2>&1 || {
  echo "ERROR: uv is required (https://docs.astral.sh/uv/)." >&2
  exit 1
}

uv venv --python "$PYTHON_VERSION" "$VENV"
uv pip install --python "$VENV/bin/python" -e "$ROOT[benchmark,test]"
"$VENV/bin/python" -m pip check 2>/dev/null || true
"$VENV/bin/python" - <<'PY'
import importlib.metadata
import sys
print("Python", sys.version.split()[0])
print("imspy-simulation", importlib.metadata.version("imspy-simulation"))
print("imspy-connector", importlib.metadata.version("imspy-connector"))
PY
"$VENV/bin/timsim" --help >/dev/null
printf '\nEnvironment ready: %s\n' "$VENV"
