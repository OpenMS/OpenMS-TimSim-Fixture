# OpenMS-TimSim-Fixture current state

## Repository

Public repository: `OpenMS/OpenMS-TimSim-Fixture`.

The published smoke baseline is tag `v0.1.0-smoke` at commit `a62b878`. It provides a
noise-free TimSim/OpenDIA fixture with simulator-only target selection, isolated
OpenDIA target/FDR lanes, independent entrapments for external-null calibration, and a
paired-hard interference stress mode.

## Validated smoke result

The 50-target smoke fixture reached 50/50 target identification in every isolated
OpenDIA target-only run. The independent-entrapment lane had one accepted global
entrapment among 51 accepted global target-labelled identities at nominal q <= 0.01;
50 entrapments are therefore adequate for smoke validation but too small for precise
1% calibration.

The three-run quantitative benchmark exposed that a single treatment replicate is not
sufficient to distinguish finite-realization variation from measurement behavior. A
larger replicate-aware design is now implemented.

## Multi-run study implementation

The current working tree adds a generalized multi-run study path:

- `scripts/generate_study.sh`
- `scripts/run_study_benchmark_lanes.sh`
- `configs/study_blueprint.toml.in`
- `configs/study_from_existing.toml.in`
- `tools/prepare_study_databases.py`
- `tools/render_study_configs.py`
- `tools/select_reference_precursors.py`
- `tools/build_study_realized_truth.py`
- `tools/benchmark_study.py`
- `tools/validate_study.py`
- `docs/STUDY_BENCHMARK.md`

Canonical defaults are 25 control + 25 treatment runs, 1,000 frozen true precursors,
1,000 independent entrapments, a 10,000-peptide TimSim blueprint sampled from a
20,000-peptide synthetic FASTA, and eight transitions per precursor.

Targets are frozen from the blueprint before treatment effects or replicate-specific
abundance are created. This intentionally allows condition-specific missingness after
selection instead of selecting it away.

Treatment design defaults to 15% up-regulated proteins, 15% down-regulated proteins,
and heterogeneous |log2FC| between 0.5 and 2.0. Independent replicate variation is
modeled with run-, protein-, and peptide-level log2 SDs of 0.05, 0.15, and 0.08.

Three truth layers are explicit:

1. `OpenSwathTimSim.study_design_truth.tsv`: population/design protein effects.
2. `OpenSwathTimSim.study_abundance_truth.tsv`: actual per-run TimSim input abundance.
3. `OpenSwathTimSim.realized_truth.tsv`: post-simulation target/run precursor truth.

`run_opendia.sh` and `benchmark_entrapment.py` no longer assume exactly three runs.
The study benchmark performs replicate-aware peptide/protein effect analysis and Welch
+ Benjamini-Hochberg differential calls while preserving design-vs-realized and
realized-vs-observed comparisons separately.

## Validation status

Repository tests after the multi-run implementation: 15 passed. Python compilation and
shell syntax validation also pass in the assistant validation environment.

The next runtime validation should be a 3-control + 3-treatment study using the full
1,000-target/1,000-entrapment library. This validates target availability, TimSim
multi-run rendering, OpenDIA multi-input behavior, and study statistics before spending
compute on the full 25+25 study.

## Next execution

Use the existing PXD017703 DIA-PASEF reference used by the smoke benchmark. Generate
`validation_3x3` with `generate_study.sh --control-runs 3 --treatment-runs 3` while
keeping 1,000 targets, 1,000 entrapments, 10,000 simulated peptides, and 20,000 FASTA
peptides. If that passes fixture validation and both OpenDIA lanes, generate the full
25+25 study with the default run counts.
