# Large-scale and single-cell TimSim benchmark tiers

The large-scale tiers are deliberately separate from the frozen 1,000-target correctness/calibration benchmark.

The correctness tier answers whether OpenDIA recovers known targets, known abundance changes, and an independent external-null entrapment set under a tightly controlled 250-protein x 4-precursor design. The scale tiers answer different questions: throughput under a large assay universe, variable protein multiplicity, and low-input/run-specific missingness.

## Shared scale reference

Both large-bulk and single-cell studies reuse one frozen simulator-only assay reference. The default census is:

- 600,000 synthetic FASTA candidate peptides;
- 300,000 TimSim blueprint peptides;
- 150,000 frozen target precursor groups;
- `variable_proteome` selection (the existing simulator-only stratified selector with no fixed precursor count per protein);
- scale raw generation starts target-only by default; external-null entrapments are added later as a separate library-only scale step because they do not alter synthetic raw generation.

Create the reference with:

```bash
TIMSIM_THREADS=8 \
./scripts/scale/prepare_scale_reference.sh \
  /path/to/reference.d \
  /path/to/scale_reference_150k
```

This is a census gate. Inspect `selection_qc/selection_summary.json` and `selection_qc/selection_report.md` before launching biological runs. If fewer than 150,000 targets qualify, increase the synthetic/blueprint peptide population; do not relax the validated high-signal filters just to hit the target count.

The reference preset runs TimSim on GPU and, after selection/library generation succeeds, compacts the blueprint directory to retain the authoritative `synthetic_data.db` rather than the large blueprint raw/intermediate payload.

## Scratch-local materialization

Large studies do not persist one perturbed source DB per biological run. Instead:

1. the immutable blueprint DB and deterministic run plan are persisted;
2. an array worker copies/materializes one perturbed source DB in node-local scratch;
3. TimSim generates one run;
4. selected-target input and realized truth are extracted immediately;
5. only the generated `.d`, compact truth/stats, and provenance are persisted;
6. the temporary source DB and nonessential TimSim run intermediates are deleted.

This architecture exists because the 25+25 / 1,000-target validation study occupied about 17 GB while each final `.d` was only about 80-84 MB. At scale, persistent source DB copies would dominate storage and NFS traffic.

The generic worker is:

```bash
USE_GPU=1 TIMSIM_THREADS=8 COMPACT_OUTPUT=1 \
./scripts/run_materialized_study_task.sh \
  /path/to/study \
  1 \
  /path/to/reference.d \
  /local/scratch
```

Schedulers should map one run ordinal to one array task.

## Large bulk v1

The recommended first production scale is 100 runs rather than immediately jumping to 500:

- 50 control + 50 treatment;
- 150,000 frozen assay targets;
- naturally variable precursor multiplicity across proteins;
- same biological effect model as the correctness tier;
- scratch-local run materialization.

Create the run plan after the scale reference passes:

```bash
./scripts/scale/plan_large_bulk_study.sh \
  /path/to/scale_reference_150k \
  /path/to/bulk_50x50_150k
```

After all run tasks finish:

```bash
./scripts/finalize_materialized_study.sh \
  /path/to/bulk_50x50_150k
```

The 100-run study is the first throughput/storage checkpoint. Expand toward 200-500 runs only after measuring raw size, node-local scratch usage, GPU/CPU utilization, queue efficiency, and downstream OpenDIA resource use.

## Single-cell calibration tier

Single-cell simulation uses the same frozen 150k assay universe but changes run-level abundance generation. The model intentionally supports true zero-event peptides:

1. blueprint abundance;
2. a global low-input log2 offset;
3. cell-size variation;
4. protein biological variation;
5. peptide residual variation;
6. Poisson event sampling;
7. zero events are retained rather than clamped to one.

Truth distinguishes:

- `BiologicalPresentInRun`: `RealizedInputEvents > 0`;
- `ObservableInSimulation`: TimSim generated positive realized precursor signal for the selected assay target.

The first calibration panel is 48 cells: 24 control + 24 treatment, balanced across four predeclared relative input offsets (-6, -8, -10, -12 by default). These are relative TimSim input levels, not picogram claims.

```bash
./scripts/scale/plan_single_cell_calibration.sh \
  /path/to/scale_reference_150k \
  /path/to/single_cell_calibration_48
```

After generation and OpenDIA analysis, choose the production input offset based on observed identification depth and missingness. The fixture must not choose the offset to maximize OpenDIA agreement or effect recovery.

## Single-cell production tier

Only after calibration, create the production cell plan. The default is 100 control + 100 treatment cells:

```bash
./scripts/scale/plan_single_cell_study.sh \
  /path/to/scale_reference_150k \
  /path/to/single_cell_100x100 \
  <CALIBRATED_LOG2_OFFSET>
```

The production tier should preserve the selected input level and all other abundance parameters. A later 500-cell performance tier can be generated from the same reference and model without changing the assay definition.

## Truth and anti-circularity

Target selection remains frozen from the high-input blueprint before condition effects and before OpenDIA. Neither raw-oracle measurements nor OpenDIA outputs influence the 150k target set.

Single-cell run-specific zeros are part of biological truth, not entrapments. Independent entrapments remain globally absent target-labelled external negatives and serve a different FDR-calibration role.


## Scale-reference selection contract

The scale/realism tier deliberately uses a different target-selection contract from the frozen correctness tier.

`protein_balanced` and `global_stratified` retain the established high-signal thresholds used for correctness/calibration fixtures. `scale_variable` is intended for dense realistic assay libraries: it requires blueprint presence, safe DIA geometry, the configured minimum usable fragment count, positive realized simulator signal, and RT-edge safety, then freezes one charge state per peptide while preserving the joint RT/mz/IM/realized-abundance distribution. It does **not** require every library precursor to exceed the correctness-tier minimum event proxy or high ion-fraction thresholds. This is intentional: a 150k realism library must contain difficult and lower-abundance assays rather than 150k preselected easy identifications.

The selector remains simulator-only and runs before any biological condition effects, raw-oracle inspection, or OpenDIA analysis. Selection is deterministic within joint RT/mz/IM/abundance strata.

A completed large blueprint may be reselected without rerunning TimSim by setting `REUSE_BLUEPRINT=1` when invoking `scripts/scale/prepare_scale_reference.sh`. The reuse path requires the existing blueprint SQLite DB and synthetic FASTA and changes only target/library/QC products.
