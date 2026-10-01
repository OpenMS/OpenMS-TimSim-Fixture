#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

GIT_SHA="${GIT_SHA:-$(git rev-parse HEAD)}"
SHORT_SHA="${GIT_SHA:0:12}"
IMAGE="${IMAGE:-openms-timsim-fixture:sha-${SHORT_SHA}}"
BUILD_UTC="${BUILD_UTC:-$(date -u +%Y-%m-%dT%H:%M:%SZ)}"

docker build \
  --file docker/Dockerfile \
  --build-arg "GIT_SHA=$GIT_SHA" \
  --build-arg "BUILD_UTC=$BUILD_UTC" \
  --tag "$IMAGE" \
  .

docker run --rm "$IMAGE" \
  python -c 'import importlib.metadata, torch; print("imspy-simulation", importlib.metadata.version("imspy-simulation")); print("torch", torch.__version__, "CUDA build", torch.version.cuda)'

echo "IMAGE=$IMAGE"
