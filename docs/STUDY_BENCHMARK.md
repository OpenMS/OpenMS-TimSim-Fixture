# Multi-run study benchmark

`generate_study.sh` builds the larger benchmark used for replicate-aware identification,
quantification, differential-abundance, and external-null FDR evaluation.

## Canonical study

The default study contains:

- 25 control runs;
- 25 treatment runs;
- 1,000 frozen true precursor groups;
- 1,000 independent known-absent entrapment groups;
- 10,000 peptides sampled into the TimSim blueprint from a 20,000-peptide synthetic FASTA;
- 8 transitions per precursor group.

The target set is frozen from a separate TimSim blueprint before treatment effects or
replicate-specific abundance variation are generated. This prevents target selection
from conditioning on treatment response or downstream OpenDIA results.

## Replicate model

Each study run reuses the same molecular blueprint: sequence, protein assignment,
charge states, precursor m/z, predicted RT, ion mobility, and acquisition geometry are
fixed. Peptide abundance changes according to independent deterministic random effects:

`log2 abundance = blueprint + condition + run + protein-within-run + peptide-within-run`

The default standard deviations are 0.05, 0.15, and 0.08 log2 units for the run,
protein, and peptide replicate components respectively.

Thirty percent of proteins are treatment-responsive by default: 15% up and 15% down.
Absolute treatment effects are sampled uniformly from 0.5 to 2.0 log2 fold-change.
The remaining proteins have a design effect of zero.

## Three quantitative truth layers

The study deliberately keeps three different quantities:

1. **Design truth** (`OpenSwathTimSim.study_design_truth.tsv`) — the predeclared
   treatment effect for each protein.
2. **Realized input abundance** (`OpenSwathTimSim.study_abundance_truth.tsv`) — the
   abundance actually written to every TimSim source database after biological
   replicate variation.
3. **Realized precursor truth** (`OpenSwathTimSim.realized_truth.tsv`) — post-simulation
   precursor signal/coordinate truth for every frozen target in every generated run.

OpenDIA quantification is compared primarily against realized condition abundance.
Design-vs-realized comparisons quantify how much biological replicate variation moves
the finite simulated cohort away from the declared population effect.

## Recommended validation sequence

Before generating 25+25 runs, exercise the architecture with 3+3 runs while retaining
the full 1,000-target/1,000-entrapment library:

```bash
./scripts/generate_study.sh \
  --reference /absolute/path/to/reference.d \
  --output "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/studies/validation_3x3" \
  --control-runs 3 \
  --treatment-runs 3 \
  --precursors 1000 \
  --entrapments 1000 \
  --simulated-peptides 10000 \
  --fasta-peptides 20000
```

Then run the isolated OpenDIA lanes:

```bash
export OPENMS_BUILD=/path/to/OpenMS-build
THREADS=12 ./scripts/run_study_benchmark_lanes.sh \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/studies/validation_3x3" \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/benchmarks/validation_3x3"
```

After that passes, generate the full study using the defaults:

```bash
./scripts/generate_study.sh \
  --reference /absolute/path/to/reference.d \
  --output "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/studies/control25_treatment25_targets1000"
```

## Benchmark outputs

The target-only OpenDIA lane is summarized by `tools/benchmark_study.py`:

- per-run target recovery and RT/IM accuracy;
- target/run measurement table;
- peptide design, realized, and observed effects;
- protein design, realized, and observed effects;
- Welch tests with Benjamini-Hochberg q-values;
- sensitivity, precision, and false-positive rate for differential calls.

The independent-entrapment lane continues to use `benchmark_entrapment.py`, which
reports run-level, global-peptide, and exported-result FDP plus q/PEP calibration.
For studies with many runs the per-run calibration tables are retained, while plots
avoid drawing dozens of overlapping individual-run curves.

## Re-benchmark completed OpenDIA lanes

OpenDIA execution and report generation are separate steps. If
`opendia_target_only/` and `opendia_entrapment/` already completed successfully,
regenerate only the study/FDR reports with:

```bash
./scripts/benchmark_study_lanes.sh STUDY_DIR BENCH_DIR
```

This helper uses the repository interpreter at `.venv/bin/python` explicitly and
does not depend on a system `python` executable being available on `PATH`.

## Quantification diagnostic gate before scaling

`benchmark_study_lanes.sh` also runs `tools/diagnose_study_quantification.py` on the
target-only study results. The diagnostic separates per-run abundance tracking from
condition-effect recovery and reports effect-regression slope, heavy-tail residuals,
per-run TimSim-input versus OpenDIA-intensity correlations, RT/IM associations, and
protein performance stratified by the number of selected precursors.

The report is written under:

```text
results/study/quantification_diagnostics/
```

Use the 3+3 diagnostic to decide whether errors are primarily finite-replicate/random
variation or systematic quantification failures before generating the 25+25 study.

## Raw-signal oracle gate

When the target-only study shows systematic quantitative compression or a heavy residual
tail, use `run_raw_signal_oracle.sh` before generating a larger cohort. The oracle reads
the already generated Bruker `.d` files and integrates the frozen target transition m/z
values inside fixed windows around simulator-realized RT and ion-mobility coordinates.
OpenDIA-selected coordinates, scores, and peak boundaries are not inputs to extraction.

First verify the installed `imspy` raw reader against one DIA frame:

```bash
RAW_ORACLE_PROBE_ONLY=1 ./scripts/run_raw_signal_oracle.sh STUDY_DIR BENCH_DIR
```

Then run the full oracle without rerunning TimSim or OpenDIA:

```bash
./scripts/run_raw_signal_oracle.sh STUDY_DIR BENCH_DIR
```

Default truth-centered integration windows are ±6 s in RT, ±0.03 1/K0 in ion mobility,
and ±25 ppm around each of the eight frozen library fragment m/z values. They can be
overridden with `RAW_ORACLE_RT_HALF_WINDOW`, `RAW_ORACLE_IM_HALF_WINDOW`, and
`RAW_ORACLE_FRAGMENT_PPM` for sensitivity analysis; the canonical result should retain
one declared parameter set rather than tuning windows against OpenDIA performance.

The raw-oracle report decomposes the measurement chain into:

1. TimSim realized input/event proxy -> raw integrated fragment signal;
2. raw integrated fragment signal -> OpenDIA intensity;
3. realized condition log2FC -> raw-oracle log2FC -> OpenDIA log2FC.

`RealizedEventProxy` is the preferred charge-state-specific absolute-signal comparator
because it includes peptide events, chromatographic frame mass, ion relative abundance,
and mobility scan mass. Condition effects remain the primary diagnostic because
sequence-specific fragment yield cancels within a precursor across runs.

The oracle is strictly diagnostic. It must never be used to select benchmark targets,
construct the assay library, choose entrapments, or tune OpenDIA on the canonical fixture.

