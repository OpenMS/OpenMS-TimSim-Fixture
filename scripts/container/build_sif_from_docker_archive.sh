#!/usr/bin/env bash
set -Eeuo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 /path/to/image.docker.tar [/path/to/output.sif]" >&2
  exit 2
fi

DOCKER_ARCHIVE="$1"
if [[ ! -f "$DOCKER_ARCHIVE" ]]; then
  echo "Docker archive does not exist: $DOCKER_ARCHIVE" >&2
  exit 2
fi

if [[ $# -eq 2 ]]; then
  SIF="$2"
else
  SIF="${DOCKER_ARCHIVE%.docker.tar}.sif"
fi

if command -v apptainer >/dev/null 2>&1; then
  ENGINE=apptainer
elif command -v singularity >/dev/null 2>&1; then
  ENGINE=singularity
else
  echo "Neither apptainer nor singularity is available in PATH." >&2
  exit 2
fi

mkdir -p "$(dirname "$SIF")"
rm -f "$SIF" "${SIF}.sha256"

"$ENGINE" build --force "$SIF" "docker-archive:$DOCKER_ARCHIVE"
sha256sum "$SIF" > "${SIF}.sha256"

"$ENGINE" exec "$SIF" \
  python -c 'import importlib.metadata, torch; assert torch.version.cuda is not None; print("imspy-simulation", importlib.metadata.version("imspy-simulation")); print("torch", torch.__version__, "CUDA build", torch.version.cuda)'

printf 'SIF=%s\n' "$SIF"
printf 'SIF_SHA256=%s\n' "$(awk '{print $1}' "${SIF}.sha256")"
