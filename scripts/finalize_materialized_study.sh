#!/usr/bin/env bash
set -Eeuo pipefail
trap 'echo "ERROR on line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STUDY_ROOT="${1:?usage: finalize_materialized_study.sh STUDY_ROOT}"
python "$ROOT/tools/finalize_materialized_study.py" --study-root "$STUDY_ROOT"
