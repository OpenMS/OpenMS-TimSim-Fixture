# Container images

The repository provides a portable CUDA-enabled generation image for TimSim-based synthetic DIA-PASEF studies.
Cluster-specific schedulers, filesystem paths, usernames, aliases, and deployment wrappers are intentionally not part of this repository.

## Docker image

`docker/Dockerfile` contains only the generation environment. It uses a CUDA/cuDNN runtime base, Python 3.11, the pinned `imspy-simulation==0.4.2`, and the fixture generation code. OpenMS/OpenDIA and the plotting/scientific-analysis stack are not installed in this image.

Build locally from a clean Git revision:

```bash
./scripts/container/build_docker.sh
```

The default local tag is:

```text
openms-timsim-fixture:sha-<12-character Git SHA>
```

Export exactly that image as a Docker archive when an offline/HPC conversion is needed:

```bash
./scripts/container/export_docker_archive.sh
```

## SIF image

A Docker archive can be converted to a SIF with either Apptainer or SingularityCE:

```bash
./scripts/container/build_sif_from_docker_archive.sh \
  openms-timsim-fixture_sha-<sha>.docker.tar \
  openms-timsim-fixture_sha-<sha>.sif
```

The helper writes a sibling `.sha256` file and verifies that the resulting image contains a CUDA-enabled PyTorch build. GPU availability itself is a property of the execution host and should be checked at runtime, for example:

```bash
apptainer exec --nv openms-timsim-fixture_sha-<sha>.sif \
  python -c 'import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))'
```

Replace `apptainer` with `singularity` on systems that provide SingularityCE.

## GitHub Actions distribution

`.github/workflows/container-images.yml` builds the Docker image on `main`, version tags, or manual dispatch.

For each Git revision the workflow:

1. builds and validates the CUDA-enabled Docker image;
2. pushes an immutable `sha-<12-character SHA>` image to GitHub Container Registry;
3. updates the `main` tag on `main`, and the version/`latest` tags for version-tag builds;
4. exports the exact locally built Docker image to a Docker archive;
5. installs Apptainer on the GitHub-hosted Linux runner;
6. builds a SIF from that Docker archive;
7. validates the SIF payload and records Docker/SIF checksums and build provenance;
8. uploads the SIF, its SHA-256 file, and CI provenance as a GitHub Actions artifact.

The registry image is therefore the normal Docker distribution surface, while the workflow artifact is the portable SIF distribution surface.

## Deployment boundary

The public repository deliberately stops at portable image production. Site-specific Slurm scripts, storage roots, GPU resource requests, module configuration, and retry launchers should live in a separate local/private deployment directory. This prevents one site's assumptions or personal filesystem paths from becoming part of the project API.

## GitHub Actions troubleshooting

The container workflow builds with BuildKit `--progress=plain` so a failed image layer is visible directly in the Actions log. To inspect only failed steps from the command line:

```bash
gh run view <run-id> --log-failed
```

The Dockerfile deliberately avoids shell/Dockerfile heredoc blocks in build-critical `RUN` instructions. This keeps the image build compatible with the BuildKit frontend used by GitHub-hosted runners. The project metadata also declares the `generation` extra explicitly; the container installs `.[generation]` and therefore fails early if the generation dependency contract is accidentally removed.
