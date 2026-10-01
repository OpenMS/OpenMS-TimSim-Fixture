#!/usr/bin/env bash
set -Eeuo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

GIT_SHA="${GIT_SHA:-$(git rev-parse HEAD)}"
SHORT_SHA="${GIT_SHA:0:12}"
IMAGE="${IMAGE:-openms-timsim-fixture:sha-${SHORT_SHA}}"
OUT="${1:-$ROOT/openms-timsim-fixture_sha-${SHORT_SHA}.docker.tar}"

docker image inspect "$IMAGE" >/dev/null
docker save --output "$OUT" "$IMAGE"
sha256sum "$OUT" > "${OUT}.sha256"
echo "$OUT"
