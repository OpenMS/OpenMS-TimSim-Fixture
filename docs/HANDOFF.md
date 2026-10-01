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

## Current next step — final 3+3 production-design validation

Generate a **fresh** 3 control + 3 treatment study using the final 250×4 selector. Do not reuse the earlier `validation_3x3` target set because it was frozen under the global precursor selector.

Required first checkpoint is the blueprint selection stage itself:

- exactly 1,000 targets;
- exactly 250 proteins;
- exactly 4 precursors/protein;
- zero singleton/doubleton/tripleton selected proteins;
- selection remains blueprint-only and OpenDIA-independent.

If fewer than 250 proteins have four eligible precursor groups, increase `--simulated-peptides` and `--fasta-peptides`; do not weaken high-signal thresholds. Once the fresh 3+3 study passes structural validation, run the existing isolated target-only and independent-entrapment OpenDIA lanes and compare identification, FDR, and protein quantitative recovery. If that passes, generate the canonical 25+25 study without further architecture changes.


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

- study root: `/home/sing/Documents/datasets/OpenMS-TimSim-Fixture/studies/production_validation_3x3`;
- target truth: `OpenSwathTimSim.realized_truth.tsv` (6,000 target/run rows);
- study manifest: `OpenSwathTimSim.study_manifest.tsv`;
- design truth: `OpenSwathTimSim.study_design_truth.tsv`.

The candidate-pool requirement is therefore now empirically established for this seed/configuration: 10k/20k is insufficient (198/250 proteins), while 14k/28k satisfies the fixed 250×4 contract. Do not lower scientific eligibility thresholds to reduce this requirement.

## Current next step — OpenDIA validation of the final production design

Do not change fixture generation further at this checkpoint. Run the existing isolated OpenDIA lanes on `production_validation_3x3`:

1. target-only lane for identification, localization, quantification, and differential effects;
2. independent-entrapment lane for external-null FDR/FDP calibration;
3. study benchmark/post-processing;
4. quantification diagnostics focused on protein realized-vs-OpenDIA recovery.

Primary acceptance questions:

- target/run recovery remains near 100% at q <= 0.01;
- independent-entrapment run/global FDP remains close to nominal 1%;
- design -> realized effects remain near the previous ~0.98 correlation;
- enforcing four precursor groups per selected protein materially improves protein realized -> OpenDIA effect recovery relative to the prior mixed-multiplicity study (previous r≈0.78 overall; prior 3+ precursor stratum r≈0.99, MAE≈0.053);
- no systematic RT/IM localization regression appears.

If this final 3+3 production-design validation passes, generate the canonical 25 control + 25 treatment benchmark using the same 14k/28k (or larger, but not smaller) candidate pool and the unchanged 250×4 selection contract. No additional architecture iteration is planned unless the final 3+3 OpenDIA validation exposes a new failure.
