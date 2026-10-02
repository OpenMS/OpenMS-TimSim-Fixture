from __future__ import annotations

import csv
import json
import sqlite3
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
SCRIPTS = ROOT / "scripts"


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def make_blueprint(path: Path) -> None:
    with sqlite3.connect(path) as con:
        con.execute("CREATE TABLE peptides (peptide_id INTEGER, sequence TEXT, protein TEXT, events INTEGER)")
        con.executemany(
            "INSERT INTO peptides VALUES (?, ?, ?, ?)",
            [
                (1, "PEPTIDEK", "P1", 1000),
                (2, "PEPTIDER", "P1", 500),
                (3, "ANOTHERK", "P2", 200),
                (4, "SEQUENCEK", "P3", 100),
            ],
        )
        con.commit()


def make_selected(path: Path) -> None:
    write_tsv(
        path,
        [
            {
                "selection_rank": 1,
                "precursor_key": "PEPTIDEK/2",
                "sequence": "PEPTIDEK",
                "protein_id": "P1",
                "charge": 2,
                "precursor_mz": 500.0,
                "assay_rt": 100.0,
                "assay_im": 1.0,
            },
            {
                "selection_rank": 2,
                "precursor_key": "ANOTHERK/2",
                "sequence": "ANOTHERK",
                "protein_id": "P2",
                "charge": 2,
                "precursor_mz": 600.0,
                "assay_rt": 200.0,
                "assay_im": 1.1,
            },
        ],
    )


def test_materialized_bulk_plan_preserves_events_without_noise(tmp_path: Path) -> None:
    blueprint = tmp_path / "blueprint.db"
    selected = tmp_path / "selected.tsv"
    fasta = tmp_path / "synthetic.fasta"
    study = tmp_path / "study"
    make_blueprint(blueprint)
    make_selected(selected)
    fasta.write_text(">P1\nMPEPTIDEKAA\n", encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "plan_materialized_study.py"),
            "--blueprint-db",
            str(blueprint),
            "--output-root",
            str(study),
            "--selected-precursors",
            str(selected),
            "--fasta",
            str(fasta),
            "--control-runs",
            "1",
            "--treatment-runs",
            "1",
            "--profile",
            "bulk",
            "--run-log2-sd",
            "0",
            "--protein-log2-sd",
            "0",
            "--peptide-log2-sd",
            "0",
            "--treatment-up-fraction",
            "0",
            "--treatment-down-fraction",
            "0",
        ],
        check=True,
    )
    output_db = tmp_path / "run.db"
    truth = tmp_path / "truth.tsv"
    stats = tmp_path / "stats.json"
    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "materialize_study_run.py"),
            "--blueprint-db",
            str(blueprint),
            "--study-manifest",
            str(study / "OpenSwathTimSim.study_manifest.tsv"),
            "--design-truth",
            str(study / "OpenSwathTimSim.study_design_truth.tsv"),
            "--selected-precursors",
            str(selected),
            "--run-ordinal",
            "1",
            "--output-db",
            str(output_db),
            "--selected-input-truth-out",
            str(truth),
            "--stats-out",
            str(stats),
        ],
        check=True,
    )
    with sqlite3.connect(output_db) as con:
        observed = dict(con.execute("SELECT sequence, events FROM peptides"))
    assert observed["PEPTIDEK"] == 1000
    assert observed["ANOTHERK"] == 200
    rows = list(csv.DictReader(truth.open(newline="", encoding="utf-8"), delimiter="\t"))
    assert {int(row["BiologicalPresentInRun"]) for row in rows} == {1}


def test_single_cell_plan_supports_true_zero_event_targets(tmp_path: Path) -> None:
    blueprint = tmp_path / "blueprint.db"
    selected = tmp_path / "selected.tsv"
    fasta = tmp_path / "synthetic.fasta"
    study = tmp_path / "study"
    make_blueprint(blueprint)
    make_selected(selected)
    fasta.write_text(">P1\nMPEPTIDEKAA\n", encoding="utf-8")

    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "plan_materialized_study.py"),
            "--blueprint-db",
            str(blueprint),
            "--output-root",
            str(study),
            "--selected-precursors",
            str(selected),
            "--fasta",
            str(fasta),
            "--control-runs",
            "1",
            "--treatment-runs",
            "1",
            "--profile",
            "single_cell",
            "--input-log2-offset",
            "-100",
            "--event-sampling",
            "poisson",
            "--allow-zero-events",
            "--run-log2-sd",
            "0",
            "--cell-size-log2-sd",
            "0",
            "--protein-log2-sd",
            "0",
            "--peptide-log2-sd",
            "0",
            "--treatment-up-fraction",
            "0",
            "--treatment-down-fraction",
            "0",
        ],
        check=True,
    )
    output_db = tmp_path / "run.db"
    truth = tmp_path / "truth.tsv"
    stats = tmp_path / "stats.json"
    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "materialize_study_run.py"),
            "--blueprint-db",
            str(blueprint),
            "--study-manifest",
            str(study / "OpenSwathTimSim.study_manifest.tsv"),
            "--design-truth",
            str(study / "OpenSwathTimSim.study_design_truth.tsv"),
            "--selected-precursors",
            str(selected),
            "--run-ordinal",
            "1",
            "--output-db",
            str(output_db),
            "--selected-input-truth-out",
            str(truth),
            "--stats-out",
            str(stats),
        ],
        check=True,
    )
    with sqlite3.connect(output_db) as con:
        assert all(events == 0 for (events,) in con.execute("SELECT events FROM peptides"))
    rows = list(csv.DictReader(truth.open(newline="", encoding="utf-8"), delimiter="\t"))
    assert {int(row["BiologicalPresentInRun"]) for row in rows} == {0}
    payload = json.loads(stats.read_text(encoding="utf-8"))
    assert payload["zero_event_peptides"] == 4
    assert payload["selected_present_targets"] == 0


def test_single_cell_offset_levels_are_balanced_within_condition(tmp_path: Path) -> None:
    blueprint = tmp_path / "blueprint.db"
    selected = tmp_path / "selected.tsv"
    fasta = tmp_path / "synthetic.fasta"
    study = tmp_path / "study"
    make_blueprint(blueprint)
    make_selected(selected)
    fasta.write_text(">P1\nMPEPTIDEKAA\n", encoding="utf-8")
    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "plan_materialized_study.py"),
            "--blueprint-db",
            str(blueprint),
            "--output-root",
            str(study),
            "--selected-precursors",
            str(selected),
            "--fasta",
            str(fasta),
            "--control-runs",
            "8",
            "--treatment-runs",
            "8",
            "--profile",
            "single_cell",
            "--input-log2-offset-levels=-6,-8,-10,-12",
            "--event-sampling",
            "poisson",
            "--allow-zero-events",
        ],
        check=True,
    )
    rows = list(csv.DictReader((study / "OpenSwathTimSim.study_manifest.tsv").open(newline="", encoding="utf-8"), delimiter="\t"))
    for condition in ("control", "treatment"):
        levels = [row["InputLevel"] for row in rows if row["Condition"] == condition]
        assert {level: levels.count(level) for level in set(levels)} == {"L1": 2, "L2": 2, "L3": 2, "L4": 2}


def test_scale_presets_keep_reference_and_run_materialization_separate() -> None:
    reference = (SCRIPTS / "scale" / "prepare_scale_reference.sh").read_text(encoding="utf-8")
    assert "PRECURSORS:-150000" in reference
    assert "SIMULATED_PEPTIDES:-300000" in reference
    assert "SELECTION_MODE:-scale_variable" in reference
    assert '--selection-mode "$SELECTION_MODE"' in reference
    assert "REUSE_BLUEPRINT" in reference
    bulk = (SCRIPTS / "scale" / "plan_large_bulk_study.sh").read_text(encoding="utf-8")
    sc = (SCRIPTS / "scale" / "plan_single_cell_calibration.sh").read_text(encoding="utf-8")
    runner = (SCRIPTS / "run_materialized_study_task.sh").read_text(encoding="utf-8")
    assert "CONTROL_RUNS:-50" in bulk and "TREATMENT_RUNS:-50" in bulk
    assert "INPUT_LEVELS:--6,-8,-10,-12" in sc
    assert '--input-log2-offset-levels="$INPUT_LEVELS"' in sc
    assert "COMPACT_OUTPUT" in runner
    assert "materialize_study_run.py" in runner


def test_materialized_finalizer_keeps_partitioned_gzip_truth(tmp_path: Path) -> None:
    import gzip

    study = tmp_path / "study"
    study.mkdir()
    selected = study / "OpenSwathTimSim.high_signal_selection.tsv"
    make_selected(selected)
    write_tsv(
        study / "OpenSwathTimSim.study_design_truth.tsv",
        [
            {"ProteinId": "P1", "TreatmentClass": "unchanged", "DesignLog2FC": "0", "DesignFoldChange": "1"},
            {"ProteinId": "P2", "TreatmentClass": "unchanged", "DesignLog2FC": "0", "DesignFoldChange": "1"},
        ],
    )
    write_tsv(
        study / "OpenSwathTimSim.study_manifest.tsv",
        [
            {"RunOrdinal": 1, "RunId": "C001", "RunName": "OpenSwathTimSim_control_001", "Condition": "control", "Replicate": 1, "InputLevel": "L1"},
            {"RunOrdinal": 2, "RunId": "T001", "RunName": "OpenSwathTimSim_treatment_001", "Condition": "treatment", "Replicate": 1, "InputLevel": "L1"},
        ],
    )
    (study / "OpenSwathTimSim.materialized_study_plan.json").write_text(
        json.dumps({"profile": "single_cell"}) + "\n", encoding="utf-8"
    )
    (study / "run_realized_truth").mkdir()
    (study / "run_stats").mkdir()
    for ordinal, run_id, run_name, condition in [
        (1, "C001", "OpenSwathTimSim_control_001", "control"),
        (2, "T001", "OpenSwathTimSim_treatment_001", "treatment"),
    ]:
        raw = study / run_name / f"{run_name}.d"
        raw.mkdir(parents=True)
        (raw / "analysis.tdf").write_bytes(b"tdf")
        (raw / "analysis.tdf_bin").write_bytes(b"bin")
        truth = study / "run_realized_truth" / f"{ordinal:03d}_{run_id}.tsv.gz"
        rows = [
            {"RunName": run_name, "ProteinId": "P1", "BiologicalPresentInRun": 1, "ObservableInSimulation": 1},
            {"RunName": run_name, "ProteinId": "P2", "BiologicalPresentInRun": 0, "ObservableInSimulation": 0},
        ]
        with gzip.open(truth, "wt", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        (study / "run_stats" / f"{ordinal:03d}_{run_id}.json").write_text(
            json.dumps(
                {
                    "profile": "single_cell",
                    "input_level": "L1",
                    "input_log2_offset": -8,
                    "positive_peptides": 1,
                    "selected_present_targets": 1,
                }
            )
            + "\n",
            encoding="utf-8",
        )

    subprocess.run(
        [sys.executable, str(TOOLS / "finalize_materialized_study.py"), "--study-root", str(study)],
        check=True,
    )
    assert (study / "OpenSwathTimSim.realized_truth_index.tsv").is_file()
    assert not (study / "OpenSwathTimSim.realized_truth.tsv").exists()
    summary = json.loads((study / "materialized_study_summary.json").read_text(encoding="utf-8"))
    assert summary["truth_layout"] == "gzip_partition_per_run"
    assert summary["partitioned_truth_rows"] == 4
    assert summary["level_summary"][0]["median_selected_targets_observable"] == 1


def test_scale_variable_selector_is_distinct_from_correctness_contract() -> None:
    selector = (TOOLS / "select_reference_precursors.py").read_text(encoding="utf-8")
    assert '"scale_variable"' in selector
    assert "scale_structural_eligible" in selector
    assert "strict_high_signal_eligible" in selector
    assert "abundance_bin" in selector
    assert "scale_variable_v1" in selector
    assert "do not relax fragment/geometry safety criteria" in selector
