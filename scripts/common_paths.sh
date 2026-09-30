#!/usr/bin/env bash
# Shared paths for OpenMS-TimSim-Fixture.
if [[ -z "${ROOT:-}" ]]; then
  ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

OPENMS_TIMSIM_DATA_ROOT="${OPENMS_TIMSIM_DATA_ROOT:-$HOME/Documents/datasets/OpenMS-TimSim-Fixture}"
FIXTURE_VENV="${FIXTURE_VENV:-$ROOT/.venv}"
DEFAULT_FIXTURE_BUILD="$OPENMS_TIMSIM_DATA_ROOT/fixtures/default"
DEFAULT_BENCHMARK_ROOT="$OPENMS_TIMSIM_DATA_ROOT/benchmarks"

normalize_path() { realpath -m -- "$1"; }
path_is_within() {
  local candidate root
  candidate="$(normalize_path "$1")"
  root="$(normalize_path "$2")"
  [[ "$candidate" == "$root" || "$candidate" == "$root"/* ]]
}
require_external_output_path() {
  local candidate="$(normalize_path "$1")"
  local label="${2:-output}"
  if path_is_within "$candidate" "$ROOT"; then
    echo "ERROR: $label must be outside the repository." >&2
    echo "  requested: $candidate" >&2
    echo "  repository: $(normalize_path "$ROOT")" >&2
    echo "  default data root: $(normalize_path "$OPENMS_TIMSIM_DATA_ROOT")" >&2
    return 1
  fi
}
