# OpenMS-TimSim-Fixture current state

## Repository

Public repository: `OpenMS/OpenMS-TimSim-Fixture`.

The published smoke baseline is tag `v0.1.0-smoke` at commit `a62b878`. The current
development state implements the replicate-aware multi-run study architecture and
post-processing recovery path.

## Canonical multi-run design

The production design is 25 control + 25 treatment runs, 1,000 frozen true
precursors arranged as exactly 250 proteins × 4 precursors/protein, 1,000 independent
entrapments, a 14,000-peptide TimSim blueprint sampled from a 28,000-peptide synthetic
FASTA candidate pool, and eight frozen transitions per precursor.

Targets are frozen from a blueprint simulation before treatment effects or
replicate-specific abundance are generated. Treatment design defaults to 15%
up-regulated proteins, 15% down-regulated proteins, heterogeneous |log2FC| between
0.5 and 2.0, and independent run-, protein-, and peptide-level abundance variation.

Three quantitative truth layers are explicit:

1. `OpenSwathTimSim.study_design_truth.tsv`: declared population-level protein effects.
2. `OpenSwathTimSim.study_abundance_truth.tsv`: actual per-run TimSim input abundance.
3. `OpenSwathTimSim.realized_truth.tsv`: post-simulation target/run precursor truth.

The target-only OpenDIA lane is authoritative for identification and quantification.
The target+independent-entrapment lane is authoritative for external-null FDR/FDP and
PEP calibration.

## 3+3 full-library validation — completed

Generation completed successfully with 3 control + 3 treatment runs, 1,000 targets,
and 1,000 independent entrapments.

- six TimSim DIA-PASEF runs completed;
- each run contains 8,363 peptides and 12,614 precursor ions;
- 1,000 frozen targets and 6,000 target/run realized-truth rows;
- RT/mz/IM identities are unchanged across all six runs;
- abundance differs for >8,300 of 8,363 peptides in every non-reference run.

OpenDIA target-only lane:

- 1,000/1,000 target precursors recovered in every run at q <= 0.01;
- 6,000/6,000 target/run observations recovered;
- per-run RT MAE 0.76–0.81 s;
- per-run IM MAE approximately 0.00185–0.00189;
- target-only OpenDIA runtime approximately 68 s wall and 1.1 GB peak memory.

Quantitative truth decomposition:

- peptide design vs realized: Pearson r=0.9796, MAE=0.1138;
- protein design vs realized: Pearson r=0.9807, MAE=0.1116;
- peptide realized vs OpenDIA: Pearson r=0.7957, MAE=0.1776, RMSE=0.4658;
- protein realized vs OpenDIA: Pearson r=0.7837, MAE=0.2096, RMSE=0.4745.

The high design-vs-realized correlations validate the study perturbation model. The
larger realized-vs-OpenDIA RMSE relative to MAE and individual examples with >1 log2
unit error indicate a heavy-tail quantification component that must be diagnosed
before the 25+25 production generation.

Differential abundance at 3+3 using BH q <= 0.05 and |log2FC| >= 0.5:

- peptide sensitivity 0.656, precision 0.888, false-positive rate 0.0339;
- protein sensitivity 0.618, precision 0.840, false-positive rate 0.0458.

This 3+3 differential sensitivity is not the production power target; the main purpose
of this stage is architecture, truth, measurement, and calibration validation.

## Independent-entrapment calibration — validated

At nominal q <= 0.01:

- run-level: TP=5,949, entrapment FP=48, empirical FDP=0.0080, recall=0.9915;
- global peptide: TP=995, entrapment FP=12, empirical FDP=0.0119, recall=0.9950;
- final OpenDIA export: TP=5,953, entrapment FP=8, empirical FDP=0.00134;
- unique exported precursor identities: 995 true + 2 entrapments, FDP=0.00201.

Discrimination/calibration:

- run-level ROC AUC 0.9982;
- global-peptide ROC AUC 0.9986;
- entrapment-vs-generated-decoy KS statistic 0.0667;
- run-level PEP Brier 0.0076 / ECE 0.0087;
- global-peptide PEP Brier 0.0099 / ECE 0.0107.

The independent external-null construction is therefore behaving close to the nominal
1% scale in the 1,000-target/1,000-entrapment validation and remains the canonical FDR
lane. The paired-hard construction remains an interference stress test only.

## Quantification diagnostic gate — completed

The 3+3 target-only outputs were diagnosed without rerunning TimSim or OpenDIA.

Per-run absolute abundance tracking:

- mean log2 TimSim-input vs OpenDIA-intensity Pearson r = 0.8921;
- minimum run Pearson r = 0.8873;
- mean regression slope = 0.8197, range 0.8086–0.8303.

Effect recovery:

- peptide realized vs OpenDIA: r=0.7957, slope=0.8108, MAE=0.1776, RMSE=0.4658;
- peptide after excluding the worst 1% residuals: r=0.8903, slope=0.8409, MAE=0.1479;
- protein realized vs OpenDIA: r=0.7837, slope=0.7788, MAE=0.2096, RMSE=0.4745;
- protein after excluding the worst 1% residuals: r=0.8567, slope=0.8121, MAE=0.1813.

Residual tails are substantial rather than confined to a single extreme observation:

- peptide median absolute residual 0.0208 log2, p95 0.9399, p99 1.8076;
- 103/1000 peptides exceed 0.5 log2 residual and 47/1000 exceed 1.0;
- protein median absolute residual 0.0339 log2, p95 1.0635, p99 1.8520;
- 76/607 proteins exceed 0.5 log2 residual and 36/607 exceed 1.0.

Large peptide effect residuals are associated with poor localization and lower signal:

- Spearman |effect residual| vs mean absolute RT error = 0.5954;
- Spearman |effect residual| vs mean absolute IM error = 0.4984;
- Spearman |effect residual| vs median realized signal = -0.5312.

Protein aggregation is strongly dependent on selected precursor count:

- one selected precursor: n=302, r=0.7001, MAE=0.2942, slope=0.7400;
- two selected precursors: n=236, r=0.8586, MAE=0.1471, slope=0.7936;
- three or more selected precursors: n=69, r=0.9899, MAE=0.0528, slope=0.9054.

Interpretation: the multi-run fixture, design truth, realized truth, identification lane,
and independent-null FDR lane are validated. The remaining issue is quantitative
measurement fidelity. The nearly identical ~0.82 per-run slopes indicate a systematic
response compression that additional biological replicates alone will not remove. The
large residual tail is additionally enriched for low-signal and poorly localized peak
groups.

## Raw-signal oracle — completed

The truth-centered raw DIA-PASEF oracle was run on the existing 3+3 study without
regenerating TimSim data or rerunning OpenDIA. The reader backend was
`imspy_core.timstof.TimsDatasetDIA`. All 6,000 target/run observations had positive raw
signal and all eight target transitions were matched for every observation.

Per-run absolute-signal scaling was highly reproducible across all six runs:

- mean TimSim input -> raw Pearson r = 0.8123; mean slope = 0.6482;
- mean realized event-proxy -> raw Pearson r = 0.8207; mean slope = 0.6574;
- mean raw -> OpenDIA Pearson r = 0.9011; mean slope = 1.0376;
- raw -> OpenDIA slopes ranged only from 1.0296 to 1.0481.

Condition-effect decomposition:

- realized input effect -> raw oracle effect: r=0.7302, slope=0.5646, MAE=0.2710;
- after trimming the worst 1% residuals: r=0.7817, slope=0.6159, MAE=0.2504;
- raw oracle effect -> OpenDIA effect: r=0.7361, slope=0.9700, MAE=0.2153;
- after trimming the worst 1% residuals: r=0.8342, slope=0.9838, MAE=0.1841.

Interpretation: OpenDIA is approximately unit-response relative to the total raw fragment
signal measured in the truth-centered oracle. The systematic amplitude compression seen
against TimSim input is already present before OpenDIA quantification. However, the raw
oracle integrates all signal falling in the target's fragment/RT/IM windows, including
cofragmented/interfering ions, so it is an observable-channel diagnostic rather than a
strict target-attributable ground truth. It must not replace design truth or simulator
input truth.

The benchmark should therefore retain three quantitative layers:

1. biological/design truth and realized TimSim input abundance;
2. truth-centered observable raw fragment signal;
3. OpenDIA reported intensity.

The current results support that OpenDIA is not the primary cause of the ~0.82 response
compression relative to simulator input. The remaining scientific question is whether
the raw attenuation is expected DIA cofragmentation/interference captured by the oracle,
or a nonlinear/raw-response property of the simulator itself.

## Interference sensitivity sweep — completed

The raw-oracle sensitivity/interference sweep completed successfully on the existing 3+3
study after the `PrecursorMz` merge hotfix. No TimSim or OpenDIA rerun was required. The
sweep evaluated 144 fixed diagnostic configurations in a single pass per raw run:

- RT half-windows: 1, 2, 4, and 6 s;
- IM half-windows: 0.01, 0.02, and 0.03 1/K0;
- fragment tolerances: 10, 15, and 25 ppm;
- transition subsets: all 8, top 6, top 4, and top 3 by simulator-only collision specificity.

Primary broad-oracle baseline (`±6 s`, `±0.03 1/K0`, `±25 ppm`, all 8 transitions):

- event-proxy -> raw Pearson r = 0.8207; slope = 0.6574;
- realized-condition effect -> raw effect r = 0.7302; slope = 0.5646; MAE = 0.2710;
- raw -> OpenDIA per-run slope = 1.0376;
- raw-effect -> OpenDIA effect slope = 0.9700.

Predeclared tight/high-specificity configuration (`±1 s`, `±0.01 1/K0`, `±10 ppm`, top 4
transitions):

- event-proxy -> raw Pearson r = 0.9754; slope = 0.9567;
- realized-condition effect -> raw effect r = 0.9717; slope = 0.9499; MAE = 0.0388;
- raw -> OpenDIA per-run slope = 0.8505;
- raw-effect -> OpenDIA effect slope = 0.8516;
- raw positive fraction = 1.0000 for all 6,000 target/run observations.

The movement toward unit response is systematic across sweep factors:

- RT tightening from 6 s to 1 s raises mean proxy->raw slope from 0.8336 to 0.9115;
- IM tightening from 0.03 to 0.01 raises mean proxy->raw slope from 0.8581 to 0.8971;
- fragment tolerance tightening from 25 to 10 ppm raises mean proxy->raw slope from 0.8351 to 0.9093;
- reducing from all 8 to top 3 low-collision transitions raises mean proxy->raw slope from 0.8175 to 0.9203.

The diagnostic gate therefore resolves to **`interference_supported`**. The broad raw-oracle
attenuation is predominantly explained by DIA cofragmentation/interference captured when
wide RT/IM/mass windows and all transitions are summed. TimSim target-attributable response
is near unit-scale under tight, low-collision extraction, and OpenDIA is not the primary source
of the original ~0.82 input-vs-observed compression.

This sweep remains diagnostic only. Canonical benchmark extraction or OpenDIA parameters must
not be tuned to maximize agreement with the oracle. The raw oracle continues to serve only as
a measurement-model diagnostic alongside biological/design truth and OpenDIA output.

## Production target-selection contract — implemented

The canonical multi-run selector now supports a hard protein-balanced production contract. `generate_study.sh` defaults to:

- selection mode: `protein_balanced`;
- 250 selected proteins;
- exactly 4 frozen precursor groups per selected protein;
- 1,000 true target precursors total;
- 1,000 independent entrapments.

Blueprint high-signal eligibility thresholds are unchanged. For each protein, the selector first forms a 2× shortlist using simulator signal quality, then chooses the final four by simulator-only fragment-collision specificity. The collision model uses the complete TimSim blueprint fragment universe and the fixed geometry established by the diagnostic sweep (6 s RT, 0.03 1/K0, 25 ppm; top four fragments summarized). OpenDIA outputs and observed raw intensities are not selection inputs.

Protein selection is balanced across RT, precursor m/z, ion mobility, charge, and SWATH window geometry. `validate_study.py` now fails closed unless the selected table exactly matches the manifest-declared protein count and precursor multiplicity, and requires the simulator-only specificity/protein-rank provenance columns.

The implementation includes focused tests for:

- signal gating before collision-specificity selection;
- geometric collision counting;
- exact protein-balanced precursor multiplicity;
- fail-closed study validation;
- production defaults/provenance in `generate_study.sh`.

The implementation adds four focused production-selector tests. The complete current checkout should report **30 passing tests** after the overlay is applied; local source-overlay validation also includes Python compilation, shell syntax, and whitespace checks.

## Final 3+3 production-design generation — completed

The first production-contract attempt used a 10,000-peptide blueprint / 20,000-peptide FASTA candidate pool. TimSim completed the blueprint successfully, but the selector found only 198 proteins with at least four eligible precursor groups. The production selector therefore failed closed before any biological runs were generated, as intended. No high-signal threshold was relaxed.

The production validation was rerun with a 14,000-peptide blueprint / 28,000-peptide FASTA candidate pool. This completed successfully and produced the fresh study:

- 3 control + 3 treatment runs;
- 1,000 frozen target precursors;
- exactly 250 selected proteins × 4 precursor groups/protein;
- 1,000 independent entrapments;
- 1,400 proteins in the synthetic design universe;
- 11,869 realized peptides/run;
- 17,867 precursor ions/run;
- 276 frames/run.

The generated fixture passed all structural/truth checks:

- `Study contract OK: 6 runs (3 control + 3 treatment), 1000 frozen targets, 1400 design proteins`;
- `Production target contract OK: 250 proteins x 4 precursors/protein`;
- design truth, realized abundance truth, and per-run realized precursor truth are internally consistent;
- RT/mz/IM identities are unchanged across all six runs;
- abundance changes versus the first run affect 11,822–11,847 of 11,869 realized peptides (median 11,841).

Canonical fresh production-validation paths:

- study root: `<DATA_ROOT>/studies/production_validation_3x3`;
- target truth: `OpenSwathTimSim.realized_truth.tsv` (6,000 target/run rows);
- study manifest: `OpenSwathTimSim.study_manifest.tsv`;
- design truth: `OpenSwathTimSim.study_design_truth.tsv`.

The candidate-pool requirement is empirically established for this seed/configuration: 10k/20k is insufficient (198/250 proteins), while 14k/28k satisfies the fixed 250×4 contract. Do not lower scientific eligibility thresholds to reduce this requirement.

## Final 3+3 production-design OpenDIA validation — completed

The fresh protein-balanced production fixture was benchmarked with the isolated target-only and independent-entrapment OpenDIA lanes.

Target-only identification/recovery:

- 1,000/1,000 target peptides at global 1% FDR;
- 250/250 selected proteins at global 1% FDR;
- 6,000/6,000 target/run observations exported;
- 100% target/run recovery in every one of the six runs;
- per-run RT MAE approximately 0.99–1.13 s;
- per-run IM MAE approximately 0.00244–0.00256 1/K0.

Production quantitative truth decomposition:

- peptide design vs realized: r=0.9819, MAE=0.1072;
- protein design vs realized: r=0.9839, MAE=0.1001;
- peptide realized vs OpenDIA: r=0.7551, MAE=0.2362;
- protein realized vs OpenDIA: r=0.8638, MAE=0.1854.

Effect-recovery regression:

- peptide slope=0.7616; trimmed worst 1% slope=0.7851;
- protein slope=0.7390; trimmed worst 1% slope=0.7863;
- per-run TimSim-input vs OpenDIA mean r=0.8443 and mean slope=0.7557.

The 4-precursor/protein contract therefore improved protein-level correlation relative to the earlier mixed-multiplicity study, but did not remove systematic effect-amplitude attenuation.

Independent-entrapment calibration remains strong and conservative at nominal q<=0.01:

- run-level TP=5,817, entrapment FP=19, empirical FDP=0.00326;
- global peptide TP=986, entrapment FP=4, empirical FDP=0.00404;
- final export TP=5,875, entrapment FP=14, empirical FDP=0.00238;
- run-level ROC AUC=0.9976; global-peptide ROC AUC=0.9991.

This validates the production fixture structurally and as an external-null FDR benchmark.

## Production raw-signal oracle / interference confirmation — completed

The truth-centered raw oracle and the full 144-configuration interference sweep were rerun on the final 14,000-blueprint / 250×4 production target set. All 6,000 target/run observations had positive raw signal and all eight target transitions were observed.

Broad oracle (`±6 s`, `±0.03 1/K0`, `±25 ppm`, all 8 transitions):

- mean event-proxy -> raw slope=0.5682;
- realized-condition effect -> raw effect r=0.6744, slope=0.5094, MAE=0.3204;
- raw -> OpenDIA mean per-run slope=1.0735;
- raw-effect -> OpenDIA effect slope=0.9722.

Predeclared tight/high-specificity oracle (`±1 s`, `±0.01 1/K0`, `±10 ppm`, top 4 simulator-ranked low-collision transitions):

- mean event-proxy -> raw r=0.9730, slope=0.9582;
- realized-condition effect -> raw effect r=0.9785, slope=0.9558, MAE=0.0387;
- raw -> OpenDIA mean per-run slope=0.7857;
- raw-effect -> OpenDIA effect r=0.7777, slope=0.8030, MAE=0.2167;
- raw positive fraction=1.0000.

The sweep again resolves to **`interference_supported`**. Tightening RT, IM, fragment-mass windows and reducing to low-collision transitions systematically moves TimSim event-proxy/raw and realized-effect/raw response toward unit slope. This closes the concern that the production fixture or TimSim fragment response itself is intrinsically compressing the declared biological effects.

The production result also refines the earlier interpretation: relative to the tight, low-collision observable signal, current OpenDIA quantification retains roughly 20% effect-amplitude attenuation (`raw-effect -> OpenDIA` slope approximately 0.80). The broad-oracle raw->OpenDIA slope near one must therefore not be interpreted as proving unit-response OpenDIA quantification; broad extraction itself contains strong DIA interference.

This remaining attenuation is now a property to be measured by the benchmark rather than a fixture-generation defect. The canonical fixture must not be tuned against OpenDIA or against the diagnostic oracle to remove it.

## Current state / next step

The benchmark architecture is now frozen for production. No further target-selector, high-signal threshold, OpenDIA parameter, or canonical assay-library tuning is required before scaling.

Canonical production design:

- 25 control + 25 treatment runs;
- 14,000-peptide TimSim blueprint sampled from a 28,000-peptide synthetic FASTA candidate pool;
- 250 selected proteins × exactly 4 frozen precursors/protein = 1,000 true targets;
- 1,000 independent entrapments;
- eight frozen library transitions per precursor;
- design/realized truth frozen before OpenDIA;
- independent-entrapment lane retained as the canonical external-null FDR/FDP benchmark;
- raw oracle and tight/high-specificity sweep retained as diagnostic-only measurement-model layers.

Next action: generate the canonical 25+25 study with the same policy family and unchanged production-selection thresholds. GPU execution may be used through the portable Docker/SIF path, but site-specific scheduling and storage remain external deployment concerns. First prove CPU-to-GPU frozen-assay parity, then generate/finalize/validate the fixture, and only then run the isolated target-only and independent-entrapment OpenDIA lanes. Use the 25+25 cohort to measure differential-abundance power, effect-size dependence, empirical false-positive behavior, replicate-count subsampling, and final quantitative/FDR performance.

## Portable GPU container distribution — implemented

The production architecture remains frozen. GPU execution is treated only as an execution optimization; target selection, high-signal thresholds, assay construction, truth, and OpenDIA parameters remain unchanged.

The public repository now contains only portable container infrastructure:

- `docker/Dockerfile`: CUDA/cuDNN runtime generation image with Python 3.11 and pinned `imspy-simulation==0.4.2`;
- `scripts/container/build_docker.sh`: local immutable Docker build keyed by Git SHA;
- `scripts/container/export_docker_archive.sh`: exact Docker archive plus SHA-256;
- `scripts/container/build_sif_from_docker_archive.sh`: generic Apptainer/SingularityCE conversion from that Docker archive;
- `.github/workflows/container-images.yml`: builds and validates Docker, publishes the immutable Docker image to GHCR, converts the exact built image to SIF, validates the SIF, and uploads SIF/checksum/provenance as a workflow artifact;
- `docs/CONTAINERS.md`: portable container build/distribution contract.

Site-specific HPC deployment has deliberately been removed from the public Git tree. Personal filesystem roots, cluster aliases, scheduler resource policy, Slurm launchers, and retry wrappers belong in a separate local/private deployment bundle. The scientific repository therefore exposes the portable execution contract (`generate_study.sh --use-gpu`, container image, SIF image) without making any one site's scheduler layout part of the project API.

Current next step: validate the GitHub Actions container build at the current revision, download or deploy the resulting immutable SIF to the chosen GPU cluster using site-local tooling, verify CPU-to-GPU frozen-assay parity, then generate the canonical 25+25 study. Do not change the frozen production selection/assay policy in response to cluster execution details.

## Container CI first-run hardening — 2026-10-01

The first public GitHub Actions container run (`36872487807`, revision `6873bab`) reached the `Build Docker image` step and exited with code 1 before any image or SIF artifact was published. The public tree exposed two CI-fragile issues that are now hardened:

- the Dockerfile used multiline heredoc-style Python `RUN` blocks without an explicit Dockerfile frontend contract; these build-time checks/provenance writers are now single-line Python commands and the Dockerfile pins `# syntax=docker/dockerfile:1.7`;
- the public `pyproject.toml` had lost the `generation` optional dependency group even though the container installs `.[generation]`; the extra is restored with `pandas>=2.1`.

The workflow now builds with BuildKit `--progress=plain` so subsequent failures expose the exact failing command in `gh run view <run-id> --log-failed`. A regression test enforces the Dockerfile/frontend, no-heredoc, generation-extra, and plain-progress contracts.

Current next step: push this CI-hardening commit, let `container-images.yml` rerun, then inspect/publish the immutable GHCR Docker image and CI-derived SIF. Do not start the 25+25 GPU generation until the image passes the SIF validation and the site-local CUDA probe / CPU-to-GPU frozen-assay parity gate.

## Container CI success and cluster deployment — 2026-10-01

GitHub Actions run `36874548664` completed successfully for the portable container workflow. The single `Docker + SIF` job completed every stage: checkout, Buildx setup, metadata resolution, GHCR login, Docker build, GHCR push, exact Docker export, Apptainer installation, Docker-archive-to-SIF conversion, SIF payload verification, provenance writing, artifact staging, and SIF artifact upload.

The portable image is therefore ready for cluster deployment. On a SingularityCE cluster, the GHCR Docker image can be pulled directly with a `docker://ghcr.io/...` URI; SingularityCE converts OCI/Docker layers into a local SIF. This route is convenient for site deployment but may not yield a byte-identical SIF to the CI artifact because the conversion occurs again on the cluster. If byte-identical SIF provenance is required, use the CI-uploaded SIF artifact (or publish the CI-built SIF itself through an OCI/ORAS registry in a future workflow iteration).

Site-specific cluster launchers remain outside the public repository. The preferred execution layout is still preparation -> GPU Slurm array -> finalization, with one biological run per array task and retry isolation by array index. Before any large production run, perform a short `singularity exec --nv` CUDA probe and preserve the container/source revision in the run specification.

## Large-scale realism benchmark proposal — discussion state

The validated 1,000-target / 250-protein × 4-precursor fixture remains the frozen correctness/FDR benchmark and should not be replaced. A second, larger scale benchmark is now under consideration to exercise realistic DIA library density and raw-data complexity: approximately 150,000 target precursor groups across 25 control + 25 treatment runs.

This scale-up should be treated as a new benchmark tier rather than a modification of the validated 1,000-target scientific contract. The existing 250 × 4 composition cannot simply be scaled to 150,000 targets because that would imply 37,500 proteins. The large benchmark needs its own protein/precursor composition, informed by a blueprint-only eligibility census, while retaining the same anti-circularity rules (selection before OpenDIA; no raw-oracle or OpenDIA feedback into target choice).

A 150,000-target library with 8 transitions/precursor already contains 1.2 million target transitions. Using 150,000 independent entrapments would double the target-labelled library to 2.4 million transitions before OpenDIA-generated decoys and would add substantial compute/storage cost with little calibration benefit. The preferred starting point is therefore approximately 150,000 true targets plus 25,000 independent entrapments; the entrapment count can be revisited after a scale pilot.

Do not launch the full 50-run 150k study immediately. First run a GPU scale pilot that measures: blueprint eligibility and feasible per-protein multiplicity; GPU memory and utilization; TimSim wall time; output `.d` size; frame-assembly scaling; and frozen-assay reproducibility. A sensible first candidate universe is on the order of 150,000 simulated peptides from roughly 300,000 FASTA candidate peptides, but the exact target protein count / precursors-per-protein contract must be chosen from the simulator-only eligibility census rather than assumed in advance.

Current next step: pull the immutable GHCR image to the target GPU cluster, run the CUDA probe, then run a blueprint-only / one-run large-scale pilot before deciding the final 150k composition and 25+25 scheduler concurrency. Keep the validated 1k benchmark as the correctness reference throughout.

## GPU-cluster CUDA compatibility probe and container correction — 2026-10-01

The first successfully published GPU container (`e35cf5307f622abf847e187677aa048377ed2129`) was pulled from GHCR on a target SingularityCE GPU cluster and converted successfully to a local SIF. Embedded provenance reported Python 3.11.13, `imspy-simulation=0.4.2`, `imspy-core=0.4.2`, `imspy-predictors=0.5.2`, `torch=2.14.1`, and `torch_cuda=13.0`.

A GPU allocation with Singularity `--nv` exposed the host NVIDIA driver correctly, but `torch.cuda.is_available()` returned false. PyTorch reported that the host driver supports CUDA 12.8 while the container's PyTorch binary was compiled for CUDA 13.0. This is a container dependency-resolution mismatch, not a TimSim, Singularity, Slurm, or GPU-allocation failure. Do not use the `e35cf5307f62` image for GPU production generation.

The portable container contract is now corrected to target CUDA 12.8 explicitly:

- CUDA base image: `nvidia/cuda:12.8.1-cudnn-runtime-ubuntu22.04`;
- PyTorch: `torch==2.11.0+cu128` from the official PyTorch cu128 index;
- `imspy-predictors 0.5.2` is compatible because it requires `torch>=2.0.0`;
- Docker build validation and SIF validation both fail unless `torch.version.cuda == "12.8"`;
- CI provenance explicitly records `pytorch_cuda=12.8`.

Current next step: build/publish a new immutable container at the corrected revision, pull that new SHA-tagged image to the GPU cluster, and repeat the short `singularity exec --nv` CUDA probe. Only after `torch.cuda.is_available()` is true and the GPU name is reported should the 1,000-target 25+25 parallel generation be submitted. Keep the prior CUDA-13 SIF only as failed-deployment provenance or remove it to reclaim space; never alias it as the production image.
