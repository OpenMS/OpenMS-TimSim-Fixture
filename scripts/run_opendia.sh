#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=common_paths.sh
source "$ROOT/scripts/common_paths.sh"

BUILD_DIR="${1:-$DEFAULT_FIXTURE_BUILD}"
OUT_DIR="${2:-$DEFAULT_BENCHMARK_ROOT/opendia}"
OPENMS_BUILD="${OPENMS_BUILD:-}"
THREADS="${THREADS:-12}"
LIBRARY_TSV="${3:-${OPENDIA_LIBRARY:-$BUILD_DIR/OpenSwathTimSim.transitions.tsv}}"

[[ -n "$OPENMS_BUILD" ]] || {
  echo "ERROR: set OPENMS_BUILD to an OpenMS build containing bin/OpenDIA" >&2
  exit 2
}
[[ -x "$OPENMS_BUILD/bin/OpenDIA" ]] || { echo "ERROR: OpenDIA not found under $OPENMS_BUILD/bin" >&2; exit 1; }
[[ -s "$LIBRARY_TSV" ]] || { echo "ERROR: OpenDIA transition list not found: $LIBRARY_TSV" >&2; exit 1; }
require_external_output_path "$OUT_DIR" "OpenDIA output directory"

mapfile -t TDF_INPUTS < <(find "$BUILD_DIR/tdfs" -maxdepth 1 -type l -name '*.d' | sort)
(( ${#TDF_INPUTS[@]} == 3 )) || {
  echo "ERROR: expected exactly three fixture .d inputs under $BUILD_DIR/tdfs" >&2
  exit 1
}

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR/OpenDIA_intermediates"

INTERMEDIATE_DIR="$OUT_DIR/OpenDIA_intermediates"

env LD_LIBRARY_PATH="$OPENMS_BUILD/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  "$OPENMS_BUILD/bin/OpenDIA" \
  -in "${TDF_INPUTS[@]}" \
  -tr "$LIBRARY_TSV" \
  -workflow:library_mode auto \
  -workflow:working_format parquet \
  -workflow:keep_intermediate_files true \
  -workflow:intermediate_dir "$INTERMEDIATE_DIR" \
  -out_dir "$OUT_DIR" \
  -TargetedDataExtraction:readOptions cacheWorkingInMemory \
  -threads "$THREADS"

printf '\nOpenDIA library: %s\n' "$LIBRARY_TSV"
printf 'OpenDIA output: %s\n' "$OUT_DIR"
printf 'Persistent OpenDIA intermediates: %s\n' "$INTERMEDIATE_DIR"
