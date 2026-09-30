# OpenMS-TimSim-Fixture

`OpenMS-TimSim-Fixture` generates reproducible synthetic DIA-PASEF data for testing OpenMS/OpenDIA.

The fixture provides explicit ground truth for three complementary benchmark questions:

1. **Identification** — recovery and localization of precursor groups known to be present.
2. **Quantification** — recovery of known abundance changes across baseline, biological-replicate, and treatment runs.
3. **Entrapment/FDR** — acceptance of known-absent target-labelled precursors while OpenDIA generates its own decoys.

Target selection is simulator-only. OpenDIA results are never used to decide which targets enter the fixture.

The simulator dependency is pinned to **`imspy-simulation==0.4.2`**. The canonical fixture disables real-reference signal and noise superposition:

```toml
add_real_data_noise = false
superimpose_on_reference = false
```

The reference Bruker `.d` therefore contributes acquisition geometry/layout rather than hidden analyte signal.

## Repository layout

```text
OpenMS-TimSim-Fixture/
├── configs/                       # TimSim condition templates
├── data/                          # TimSim modification configuration
├── scripts/
│   ├── setup.sh                   # create the uv environment
│   ├── generate_fixture.sh        # build a complete synthetic fixture
│   ├── rebuild_entrapments.sh     # regenerate only the external-null library
│   ├── run_opendia.sh             # run OpenDIA with a selected transition library
│   ├── run_benchmark_lanes.sh     # isolated target-only and entrapment OpenDIA lanes
│   ├── benchmark_fixture.sh       # ID/quantification/correctness/FDR reports
│   └── check_fixture.sh           # structural and truth validation
├── tools/                         # scientific generation and benchmark tools
├── tests/                         # fast unit/contract tests
├── docs/
│   ├── DESIGN.md
│   ├── REPRODUCIBILITY.md
│   └── HANDOFF.md
└── pyproject.toml
```

Generated data deliberately lives outside the repository. By default:

```text
~/Documents/datasets/OpenMS-TimSim-Fixture/
```

Override this with `OPENMS_TIMSIM_DATA_ROOT`.

## 1. Set up the environment

The project uses `uv` and Python 3.11 by default.

```bash
cd OpenMS-TimSim-Fixture
./scripts/setup.sh
```

This creates `.venv` and installs the project, `imspy-simulation==0.4.2`, benchmark dependencies, and pytest.

## 2. Generate the fixture

A Bruker DIA-PASEF `.d` reference is required so TimSim can reproduce a realistic acquisition layout.

```bash
./scripts/generate_fixture.sh \
  --reference /absolute/path/to/reference.d \
  --output "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/fixtures/default"
```

Default scientific design:

- 500 frozen true precursor groups;
- 500 known-absent independent entrapment groups;
- 8 input transitions per group;
- 3 runs: baseline, biological replicate, treatment;
- explicit treatment protein fold changes;
- no real-data signal superposition;
- simulator-only high-signal selection;
- OpenDIA decoy generation downstream.

If fewer than 500 candidates meet the fixed high-signal criteria, increase the candidate pool rather than weakening the correctness thresholds:

```bash
./scripts/generate_fixture.sh \
  --reference /absolute/path/to/reference.d \
  --output "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/fixtures/default" \
  --simulated-peptides 6000 \
  --fasta-peptides 12000
```

## 3. Fixture contract

The generated fixture contains:

```text
fixture/
├── OpenSwathTimSim.high_signal_selection.tsv
├── OpenSwathTimSim.realized_truth.tsv
├── OpenSwathTimSim.selection_candidates.tsv
├── OpenSwathTimSim.condition_truth.tsv
├── OpenSwathTimSim.ground_truth.tsv
├── OpenSwathTimSim.target_only.transitions.tsv
├── OpenSwathTimSim.transitions.tsv
├── OpenSwathTimSim.entrapment_truth.tsv
├── fixture_manifest.json
├── selection_qc/
├── tdfs/
├── OpenSwathTimSim_01_baseline/
├── OpenSwathTimSim_02_biological_replicate/
└── OpenSwathTimSim_03_treatment/
```

`fixture_manifest.json` records the TimSim version, seed plan, selection thresholds, condition design, input counts, and SHA-256 hashes of key truth/library artifacts.

Validate an existing fixture with:

```bash
./scripts/check_fixture.sh \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/fixtures/default"
```

## 4. Run OpenDIA

Build OpenMS/OpenDIA separately and export the build path:

```bash
export OPENMS_BUILD=/path/to/OpenMS-build
```

Use separate OpenDIA lanes so the external-null library cannot perturb the identification/quantification scoring population:

```bash
THREADS=12 ./scripts/run_benchmark_lanes.sh \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/fixtures/default" \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/benchmarks/default"
```

This creates:

- `opendia_target_only/` from `OpenSwathTimSim.target_only.transitions.tsv` for identification, localization, and quantification;
- `opendia_entrapment/` from `OpenSwathTimSim.transitions.tsv` for external-null FDR calibration;
- `results/` containing the benchmark reports for both lanes.

The separation is intentional. External target-labelled negatives change the OpenDIA/Percolator scoring population and can alter target peak-group selection. Identification and quantification therefore use the target-only scoring lane.

For a one-off invocation, `run_opendia.sh` accepts an explicit transition library:

```bash
./scripts/run_opendia.sh FIXTURE OUT_DIR LIBRARY.tsv
```

The runner retains `workflow.oswpq` under `OpenDIA_intermediates/` for run-level and global-peptide analysis.

## 5. Benchmark OpenDIA

`run_benchmark_lanes.sh` runs benchmarking automatically. To benchmark existing lane outputs manually:

```bash
./scripts/benchmark_fixture.sh \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/fixtures/default" \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/benchmarks/default/opendia_target_only" \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/benchmarks/default/results" \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/benchmarks/default/opendia_entrapment"
```

The benchmark reports:

- target recovery and RT/IM localization;
- control-replicate and treatment quantitative recovery;
- explicit correctness failures;
- run-level, global-peptide, and final-export entrapment FDP;
- target/entrapment discrimination and PEP calibration.

## Entrapment modes

The default `independent` entrapment construction is the canonical external null. Entrapment sequences are absent from the TimSim population, use their own precursor m/z and fragment ions, and receive deterministic target-like RT/IM coordinates that remain compatible with the DIA-PASEF acquisition geometry.

`paired_hard` is available as an optional adversarial interference stress test. It deliberately places an absent peptide annotation at the precursor m/z, RT, and IM of a real target. Its false-assignment fraction is therefore an interference-robustness metric, not a nominal FDR calibration measurement.

Entrapments remain `Decoy=0`; OpenDIA generates its own decoys independently.

## Validation status

A 50-target/50-entrapment DIA-PASEF smoke fixture has been exercised end-to-end with `imspy-simulation==0.4.2` and isolated OpenDIA benchmark lanes.

Current smoke-scale observations:

- target-only identification recovery: 50/50 in baseline, control replicate, and treatment;
- independent-entrapment lane: one accepted global entrapment at nominal q <= 0.01 (1/51 empirical FDP at global peptide level);
- run-level independent-entrapment FDP: 3/153 at nominal q <= 0.01;
- treatment fold-change recovery remains under investigation despite complete identification.

Fifty external negatives are sufficient for a smoke test but not for a precise 1% FDR estimate; larger fixtures and multiple independent realizations are required for calibration studies.

## Multi-run study benchmark

For a replicate-aware benchmark, use `generate_study.sh`. The canonical larger study is **25 control + 25 treatment runs**, with **1,000 frozen targets** and **1,000 independent entrapments**. Targets are selected from a separate TimSim blueprint before treatment effects or replicate abundance variation are created.

Validate the architecture first with 3+3 runs at the full library size:

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

Then benchmark the two isolated OpenDIA lanes:

```bash
export OPENMS_BUILD=/path/to/OpenMS-build
THREADS=12 ./scripts/run_study_benchmark_lanes.sh \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/studies/validation_3x3" \
  "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/benchmarks/validation_3x3"
```

The study records three quantitative truth layers: predeclared protein `DesignLog2FC`, per-run realized TimSim input abundance, and post-simulation realized precursor truth. This separates finite-cohort biological variation from OpenDIA measurement error. See `docs/STUDY_BENCHMARK.md` for the full design.

For the production 25+25 study, omit the run-count overrides:

```bash
./scripts/generate_study.sh \
  --reference /absolute/path/to/reference.d \
  --output "$HOME/Documents/datasets/OpenMS-TimSim-Fixture/studies/control25_treatment25_targets1000"
```

## Reproducibility policy

The synthetic-proteome seed is separate from the simulation seed so multiple synthetic DIA realizations can reuse the same peptide universe. TimSim's `sample_seed`, condition perturbation seed, software versions, and truth hashes are written to the fixture manifest.

The project uses the pinned upstream simulator directly. Reproducibility is tested as a property of that dependency. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Tests

Fast repository tests do not require generating a complete TimSim run:

```bash
.venv/bin/python -m pytest
```

They cover synthetic FASTA replay, tryptic peptide construction, condition-coordinate preservation, noise-free TimSim templates, fixture-manifest contracts, OpenDIA intermediate retention, entrapment construction, and isolated benchmark lanes.

## Design notes

See [`docs/DESIGN.md`](docs/DESIGN.md) for the scientific contract and [`docs/HANDOFF.md`](docs/HANDOFF.md) for the current development status and next validation steps.
