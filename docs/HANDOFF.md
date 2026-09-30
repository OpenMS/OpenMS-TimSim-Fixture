# OpenMS-TimSim-Fixture development handoff

Updated: 2026-09-30

## Project goal

`OpenMS-TimSim-Fixture` produces reproducible synthetic DIA-PASEF fixtures for OpenMS/OpenDIA with explicit truth for:

- precursor identification and localization;
- quantitative recovery across baseline, biological-replicate, and treatment runs;
- external-null error calibration with known-absent target-labelled entrapments.

Target selection is based only on TimSim truth and is frozen before OpenDIA runs. Entrapments remain `Decoy=0`; OpenDIA generates its own decoys.

## Dependency and fixture contract

- Python: 3.11+
- `imspy-simulation==0.4.2`
- canonical TimSim configs set `add_real_data_noise=false`
- canonical TimSim configs set `superimpose_on_reference=false`
- generated fixtures live outside the source repository
- fixture manifests record randomization, software versions, selection settings, and SHA-256 hashes for key truth/library artifacts

The default fixture design requests 500 true high-signal precursor groups and 500 independent entrapments, with eight transitions per group and three DIA-PASEF runs.

## Benchmark architecture

The benchmark uses two isolated OpenDIA scoring lanes over the same frozen TimSim data:

1. **Target-only lane** — `OpenSwathTimSim.target_only.transitions.tsv`; authoritative for identification, localization, correctness, and quantification.
2. **Independent-entrapment lane** — `OpenSwathTimSim.transitions.tsv`; authoritative for external-null FDR/FDP calibration.

`paired_hard` is available separately as an adversarial interference challenge. It is not used as the canonical FDR null.

`scripts/run_benchmark_lanes.sh` runs both lanes and dispatches each output to the appropriate benchmark stage.

## Validated smoke configuration

A 50-target/50-entrapment DIA-PASEF fixture has been exercised end-to-end with the isolated benchmark architecture.

### Repository validation

- pytest: 11 passed
- retained Python tools compile
- supported shell scripts pass `bash -n`

### Target-only identification/localization

All 50 frozen targets were identified in all three runs:

| Run | Recovery | RT MAE (s) | IM MAE |
|---|---:|---:|---:|
| baseline | 50/50 | 0.2509 | 0.000542 |
| control replicate | 50/50 | 0.2500 | 0.000540 |
| treatment | 50/50 | 0.3568 | 0.000497 |

Correctness postprocessing reports 50/50 complete across all three runs. Four rows exceed configured localization thresholds; none are missing identifications.

### Quantification

Control-replicate quantitative agreement is strong:

- complete precursor pairs: 50
- expected-vs-observed log2 ratio Pearson r: 0.9887
- MAE: 0.0027
- bias: -0.0022

Treatment-effect recovery is the main unresolved scientific issue:

- peptides quantified in all three runs: 50
- peptide expected-vs-observed log2FC Pearson r: 0.5118
- peptide MAE: 0.1078
- proteins quantified in all three runs: 49
- protein expected-vs-observed log2FC Pearson r: 0.5065
- protein MAE: 0.1573
- peptide three-class accuracy at |log2FC| >= 0.5: 0.98
- changed-protein direction accuracy: 1.0
- protein altered-vs-unchanged AUC: 0.973

Because the treatment correlation remains low in the isolated target-only lane, the next investigation should focus on quantitative truth/aggregation and OpenDIA intensity semantics rather than entrapment scoring effects.

### Independent external-null calibration

At nominal q <= 0.01:

| Context | TP | Entrapment FP | Empirical FDP |
|---|---:|---:|---:|
| run-level peak groups | 150 | 3 | 0.0196 |
| global peptide identities | 50 | 1 | 0.0196 |
| final export rows | 150 | 3 | 0.0196 |

Additional smoke metrics:

- run-level target-vs-entrapment ROC AUC: 0.9996
- global-peptide ROC AUC: 0.9996
- entrapment-vs-generated-decoy KS statistic: 0.1504
- Wasserstein distance: 0.2169
- run-level PEP Brier score: 0.0147
- run-level PEP ECE: 0.0145
- global-peptide PEP Brier score: 0.0322
- global-peptide PEP ECE: 0.0345

With 50 external negatives, a single global false positive already corresponds to 1/51 = 1.96% FDP. This fixture is therefore a functional smoke test, not a precise calibration study.

## Current next steps

1. Investigate treatment quantitative recovery in the target-only lane.
   - compare the TimSim abundance truth used by the benchmark with the exact quantity represented by OpenDIA exported intensity;
   - inspect the peptides/proteins with the largest log2FC residuals;
   - determine whether aggregation, peak-group selection, normalization, or treatment truth construction explains the reduced correlation;
   - repeat the identical target-only OpenDIA run if needed to measure scoring variability at this small library size.
2. Keep identification/localization and external-null calibration unchanged while investigating quantification.
3. Once treatment quantification is understood, scale to approximately 500 true targets + 500 independent entrapments.
4. Use multiple independent simulation/entrapment seeds for FDP variability when assessing nominal FDR calibration.
5. Keep `paired_hard` as a separate interference-robustness benchmark rather than mixing it with canonical FDR calibration.

## Core commands

Environment:

```bash
./scripts/setup.sh
```

Generate a fixture:

```bash
./scripts/generate_fixture.sh \
  --reference /absolute/path/to/reference.d \
  --output /absolute/path/to/fixture
```

Validate a fixture:

```bash
./scripts/check_fixture.sh /absolute/path/to/fixture
```

Run isolated OpenDIA lanes and benchmark them:

```bash
export OPENMS_BUILD=/absolute/path/to/OpenMS-build
THREADS=12 ./scripts/run_benchmark_lanes.sh \
  /absolute/path/to/fixture \
  /absolute/path/to/benchmark
```
