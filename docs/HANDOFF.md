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

## Raw-signal oracle implementation — ready to run

The next diagnostic implementation is complete in source and does not require TimSim or
OpenDIA to be rerun.

New entry points:

- `tools/extract_raw_signal_oracle.py`: reads generated Bruker DIA-PASEF `.d` data using
  the TimSim-compatible `imspy`/`imspy-core` raw reader; assigns each frozen precursor
  to its simulator DIA window group; and integrates all frozen quantifying fragments in
  fixed truth-centered RT/IM/mass windows.
- `tools/benchmark_raw_signal_oracle.py`: decomposes per-run absolute signal and
  condition-effect response into TimSim -> raw and raw -> OpenDIA components.
- `scripts/run_raw_signal_oracle.sh`: probe/extract/benchmark convenience wrapper using
  the repository `.venv/bin/python` interpreter.

Canonical initial oracle parameters are ±6 s RT, ±0.03 1/K0 IM, and ±25 ppm fragment
m/z. Extraction uses `RealizedRTApex`, `RealizedIMApexSQLite`, the frozen target-only
transition library, and the simulator DIA window-group geometry. OpenDIA-selected
coordinates, scores, q-values, and peak boundaries are excluded.

The oracle writes per-target/run and per-transition raw measurements, per-run response
regressions, condition-level raw log2FC values, a Markdown report, and a JSON summary.
A reader-probe mode validates the installed `imspy` frame API before scanning all six
runs.

Repository validation after this implementation: 21 fast tests pass, including raw
transition integration, raw-frame adapter normalization, decomposition of a synthetic
realized->raw->OpenDIA response, and a contract that the oracle runner does not invoke
OpenDIA or study generation.

## Current next step

Run the reader probe on the existing 3+3 study, then run the full raw oracle. Do not
regenerate TimSim and do not rerun either OpenDIA lane.

The decisive comparisons are:

1. per-run `RealizedEventProxy` vs raw-oracle intensity slope;
2. per-run raw-oracle intensity vs OpenDIA intensity slope;
3. realized condition log2FC vs raw-oracle log2FC;
4. raw-oracle log2FC vs OpenDIA log2FC.

If TimSim/event-proxy -> raw is close to linear while raw -> OpenDIA retains the ~0.82
compression, investigate OpenDIA extraction/integration before scaling. If the raw
response is itself compressed, promote raw realized signal to the measurement-level
quantitative truth and characterize the simulator response before changing OpenDIA.

The full 25+25 production study remains on hold until this gate is resolved. For the
eventual protein benchmark, separately evaluate a simulator-only target-selection
contract with at least three high-quality precursors per benchmark protein; this must
not depend on OpenDIA results.
