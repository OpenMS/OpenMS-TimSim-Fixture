# syntax=docker/dockerfile:1.7
# CUDA-runtime-only image for TimSim generation. Development compilers/toolchains are not retained.
FROM ghcr.io/astral-sh/uv:0.8.22 AS uv
FROM nvidia/cuda:12.8.1-cudnn-runtime-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_SHA=unknown
ARG BUILD_UTC=unknown
ARG PYTORCH_VERSION=2.11.0+cu128
ARG PYTORCH_CUDA_INDEX=https://download.pytorch.org/whl/cu128

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

# The generation extra is intentionally smaller than the full benchmark environment.
# Pin PyTorch to a CUDA 12.8 wheel so the image remains compatible with hosts whose
# NVIDIA driver exposes CUDA 12.8 but cannot initialize CUDA 13.x binaries.
RUN uv venv --python 3.11 /opt/venv \
 && UV_EXTRA_INDEX_URL="$PYTORCH_CUDA_INDEX" UV_INDEX_STRATEGY=unsafe-best-match \
      uv pip install --python /opt/venv/bin/python \
        "torch==$PYTORCH_VERSION" ".[generation]" \
 && /opt/venv/bin/python -c 'import importlib.metadata as m, torch; assert torch.version.cuda == "12.8", f"Expected CUDA 12.8 PyTorch, got {torch.version.cuda}"; print("Python runtime OK"); print("imspy-simulation", m.version("imspy-simulation")); print("torch", torch.__version__, "CUDA", torch.version.cuda)'

RUN printf 'git_sha=%s\nbuild_utc=%s\n' "$GIT_SHA" "$BUILD_UTC" > /opt/openms-timsim-fixture/container-build-info.txt \
 && /opt/venv/bin/python -c 'import importlib.metadata as m, platform, torch; names=("openms-timsim-fixture","imspy-simulation","imspy-connector","imspy-core","imspy-predictors","torch","pandas"); print("python=" + platform.python_version()); [print(name + "=" + m.version(name)) for name in names if name in {d.metadata["Name"] for d in m.distributions()}]; print("torch_cuda=" + str(torch.version.cuda))' >> /opt/openms-timsim-fixture/container-build-info.txt

CMD ["/bin/bash"]
