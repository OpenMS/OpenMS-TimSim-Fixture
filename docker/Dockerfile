# CUDA-runtime-only image for TimSim generation. Development compilers/toolchains are not retained.
FROM ghcr.io/astral-sh/uv:0.8.22 AS uv
FROM nvidia/cuda:12.4.1-cudnn-runtime-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_SHA=unknown
ARG BUILD_UTC=unknown

LABEL org.opencontainers.image.title="OpenMS-TimSim-Fixture" \
      org.opencontainers.image.source="https://github.com/OpenMS/OpenMS-TimSim-Fixture" \
      org.opencontainers.image.revision="$GIT_SHA"

COPY --from=uv /uv /usr/local/bin/uv

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ca-certificates \
      libglib2.0-0 \
      libgl1 \
      libgomp1 \
      libsqlite3-0 \
 && rm -rf /var/lib/apt/lists/*

ENV UV_PYTHON_INSTALL_DIR=/opt/uv-python \
    UV_PYTHON_PREFERENCE=only-managed \
    UV_COMPILE_BYTECODE=1 \
    VIRTUAL_ENV=/opt/venv \
    FIXTURE_VENV=/opt/venv \
    OPENMS_TIMSIM_DATA_ROOT=/work \
    PATH=/opt/venv/bin:/usr/local/bin:/usr/bin:/bin \
    MPLCONFIGDIR=/tmp/matplotlib \
    XDG_CACHE_HOME=/tmp/.cache

RUN uv python install 3.11

WORKDIR /opt/openms-timsim-fixture
COPY . /opt/openms-timsim-fixture

# generation extra is intentionally smaller than the full benchmark environment.
RUN uv venv --python 3.11 /opt/venv \
 && uv pip install --python /opt/venv/bin/python ".[generation]" \
 && /opt/venv/bin/python - <<'PY'
import importlib.metadata
import torch
assert torch.version.cuda is not None, "CUDA-enabled PyTorch was not installed"
print("Python runtime OK")
print("imspy-simulation", importlib.metadata.version("imspy-simulation"))
print("torch", torch.__version__, "CUDA", torch.version.cuda)
PY

RUN { \
      printf 'git_sha=%s\n' "$GIT_SHA"; \
      printf 'build_utc=%s\n' "$BUILD_UTC"; \
      /opt/venv/bin/python - <<'PY'; \
import importlib.metadata
import platform
import torch
print(f"python={platform.python_version()}")
for name in ("openms-timsim-fixture", "imspy-simulation", "imspy-connector", "imspy-core", "imspy-predictors", "torch", "pandas"):
    try:
        print(f"{name}={importlib.metadata.version(name)}")
    except importlib.metadata.PackageNotFoundError:
        pass
print(f"torch_cuda={torch.version.cuda}")
PY
    } > /opt/openms-timsim-fixture/container-build-info.txt

CMD ["/bin/bash"]
