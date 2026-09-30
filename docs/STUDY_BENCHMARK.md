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
