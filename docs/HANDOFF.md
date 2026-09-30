# OpenMS-TimSim-Fixture current state

## Repository

Public repository: `OpenMS/OpenMS-TimSim-Fixture`.

The published smoke baseline is tag `v0.1.0-smoke` at commit `a62b878`. The current
development state implements the replicate-aware multi-run study architecture and
post-processing recovery path.

## Canonical multi-run design

The production design is 25 control + 25 treatment runs, 1,000 frozen true
precursors, 1,000 independent entrapments, a 10,000-peptide TimSim blueprint sampled
from a 20,000-peptide synthetic FASTA, and eight transitions per precursor.

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

## Current next step

Before generating the full 25+25 production study, run an oracle sensitivity/interference
decomposition on the existing 3+3 raw files. No TimSim or OpenDIA rerun is required.

The next diagnostic should evaluate the same raw data across multiple truth-centered
extraction windows and transition subsets, for example:

- RT half-windows: 1, 2, 4, and 6 s;
- IM half-windows: 0.01, 0.02, and 0.03 1/K0;
- fragment ppm windows: 10, 15, and 25 ppm;
- all eight transitions versus high-specificity / low-collision transition subsets.

For each configuration, report event-proxy -> raw and raw -> OpenDIA slopes/correlations,
plus condition-effect slopes. If event-proxy -> raw approaches unit slope as extraction
becomes more specific, the attenuation is primarily cofragmentation/interference and the
current biological design is suitable for scaling. If the slope remains near ~0.65 across
tight, specific windows, inspect TimSim's target-attributable fragment generation / raw
response before producing the 25+25 canonical benchmark.

The production protein-level benchmark should also use a simulator-only target-selection
contract that yields multiple high-quality precursors per protein. The current 3+
precursor stratum already shows near-ideal protein effect recovery and should guide the
final 1,000-precursor composition.
