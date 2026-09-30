from __future__ import annotations

import importlib.util
import sqlite3
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"


def load_tool(name: str):
    path = TOOLS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_synthetic_fasta_is_reproducible(tmp_path: Path) -> None:
    script = TOOLS / "generate_synthetic_fasta.py"
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    c = tmp_path / "c.fasta"

    for output, seed in ((a, 1729), (b, 1729), (c, 1730)):
        subprocess.run(
            [sys.executable, str(script), "--out", str(output), "--peptides", "100", "--seed", str(seed)],
            check=True,
            capture_output=True,
            text=True,
        )

    assert a.read_bytes() == b.read_bytes()
    assert a.read_bytes() != c.read_bytes()


def test_generated_peptides_are_unique_and_tryptic() -> None:
    module = load_tool("generate_synthetic_fasta")
    import random

    rng = random.Random(1729)
    used: set[str] = set()
    peptides = [module.make_peptide(rng, i, used) for i in range(250)]

    assert len(peptides) == len(set(peptides))
    assert all(9 <= len(p) <= 18 for p in peptides)
    assert all(p[-1] in "KR" for p in peptides)
    assert all("K" not in p[:-1] and "R" not in p[:-1] for p in peptides)


def test_config_templates_disable_reference_noise() -> None:
    for path in sorted((ROOT / "configs").glob("*.toml.in")):
        text = path.read_text(encoding="utf-8")
        assert "add_real_data_noise = false" in text
        assert "superimpose_on_reference = false" in text


def test_project_pins_supported_timsim_release() -> None:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert '"imspy-simulation==0.4.2"' in text


def test_prepare_condition_database_preserves_coordinates(tmp_path: Path) -> None:
    module = load_tool("prepare_condition_databases")
    baseline = tmp_path / "baseline.db"
    replicate = tmp_path / "replicate.db"
    treatment = tmp_path / "treatment.db"
    truth = tmp_path / "truth.tsv"

    with sqlite3.connect(baseline) as con:
        con.execute(
            "CREATE TABLE peptides (peptide_id INTEGER PRIMARY KEY, sequence TEXT, protein TEXT, events REAL, retention_time REAL)"
        )
        con.execute(
            "CREATE TABLE ions (peptide_id INTEGER, charge INTEGER, mz REAL, ion_mobility REAL)"
        )
        for i in range(1, 21):
            con.execute(
                "INSERT INTO peptides VALUES (?, ?, ?, ?, ?)",
                (i, f"PEPTIDE{i}K", f"P{(i - 1) // 2}", 100000.0 + i, 10.0 + i),
            )
            con.execute("INSERT INTO ions VALUES (?, ?, ?, ?)", (i, 2, 500.0 + i, 1.0 + i / 1000))

    rc = module.main if False else None  # keep module import covered; CLI is tested below
    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "prepare_condition_databases.py"),
            "--baseline-db", str(baseline),
            "--replicate-db", str(replicate),
            "--treatment-db", str(treatment),
            "--truth-out", str(truth),
            "--seed", "1234",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    def coords(path: Path):
        with sqlite3.connect(path) as con:
            peptides = con.execute("SELECT peptide_id, sequence, retention_time FROM peptides ORDER BY peptide_id").fetchall()
            ions = con.execute("SELECT peptide_id, charge, mz, ion_mobility FROM ions ORDER BY peptide_id, charge").fetchall()
        return peptides, ions

    base_coords = coords(baseline)
    assert coords(replicate) == base_coords
    assert coords(treatment) == base_coords
    assert truth.is_file() and truth.stat().st_size > 0


def test_high_signal_validator_accepts_current_fixture_kind() -> None:
    module = load_tool("validate_high_signal_selection")
    assert module.FIXTURE_KIND == "openms_timsim_identification_quantification_entrapment"
    assert "high_signal_correctness" not in (TOOLS / "validate_high_signal_selection.py").read_text(encoding="utf-8")


def test_opendia_runner_persists_workflow_archive() -> None:
    text = (ROOT / "scripts" / "run_opendia.sh").read_text(encoding="utf-8")
    assert '-workflow:keep_intermediate_files true' in text
    assert '-workflow:intermediate_dir "$INTERMEDIATE_DIR"' in text
    assert 'INTERMEDIATE_DIR="$OUT_DIR/OpenDIA_intermediates"' in text


def test_independent_entrapments_are_canonical_default() -> None:
    text = (ROOT / "scripts" / "generate_fixture.sh").read_text(encoding="utf-8")
    assert 'ENTRAPMENT_MODE="${ENTRAPMENT_MODE:-independent}"' in text
    design = (ROOT / "docs" / "DESIGN.md").read_text(encoding="utf-8")
    assert "canonical external null" in design
    assert "adversarial interference stress test" in design


def test_entrapment_rebuild_helper_preserves_timsim_fixture() -> None:
    text = (ROOT / "scripts" / "rebuild_entrapments.sh").read_text(encoding="utf-8")
    assert "generate_entrapment_library.py" in text
    assert "OpenSwathTimSim.target_only.transitions.tsv" in text
    assert "OpenSwathTimSim_01_baseline/synthetic_data.db" in text
    assert "run_timsim" not in text
    assert "generate_fixture.sh" not in text


def test_opendia_runner_accepts_explicit_library() -> None:
    text = (ROOT / "scripts" / "run_opendia.sh").read_text(encoding="utf-8")
    assert 'LIBRARY_TSV="${3:-${OPENDIA_LIBRARY:-$BUILD_DIR/OpenSwathTimSim.transitions.tsv}}"' in text
    assert '-tr "$LIBRARY_TSV"' in text


def test_benchmark_lanes_isolate_target_quant_from_entrapment_scoring() -> None:
    text = (ROOT / "scripts" / "run_benchmark_lanes.sh").read_text(encoding="utf-8")
    assert "OpenSwathTimSim.target_only.transitions.tsv" in text
    assert "OpenSwathTimSim.transitions.tsv" in text
    assert "opendia_target_only" in text
    assert "opendia_entrapment" in text

    benchmark = (ROOT / "scripts" / "benchmark_fixture.sh").read_text(encoding="utf-8")
    assert 'TARGET_OPENDIA_DIR="${2:-$DEFAULT_BENCHMARK_ROOT/opendia_target_only}"' in benchmark
    assert 'ENTRAPMENT_OPENDIA_DIR="${4:-$TARGET_OPENDIA_DIR}"' in benchmark
    assert '--opendia-dir "$TARGET_OPENDIA_DIR"' in benchmark
    assert '--opendia-dir "$ENTRAPMENT_OPENDIA_DIR"' in benchmark


def _make_blueprint_db(path: Path, peptides: int = 30) -> None:
    with sqlite3.connect(path) as con:
        con.execute(
            "CREATE TABLE peptides (peptide_id INTEGER PRIMARY KEY, sequence TEXT, protein TEXT, events REAL, retention_time REAL)"
        )
        con.execute(
            "CREATE TABLE ions (peptide_id INTEGER, charge INTEGER, mz REAL, ion_mobility REAL)"
        )
        for i in range(1, peptides + 1):
            con.execute(
                "INSERT INTO peptides VALUES (?, ?, ?, ?, ?)",
                (i, f"STUDYPEPTIDE{i}K", f"P{(i - 1) // 3:03d}", 100000.0 + 100 * i, 10.0 + i / 10),
            )
            con.execute(
                "INSERT INTO ions VALUES (?, ?, ?, ?)",
                (i, 2 + (i % 2), 450.0 + i, 0.9 + i / 1000),
            )


def test_prepare_study_databases_is_deterministic_and_separates_truth_layers(tmp_path: Path) -> None:
    blueprint = tmp_path / "blueprint.db"
    _make_blueprint_db(blueprint)

    def run_once(root: Path) -> tuple[str, str, str]:
        manifest = root / "study.tsv"
        design = root / "design.tsv"
        abundance = root / "abundance.tsv"
        subprocess.run(
            [
                sys.executable,
                str(TOOLS / "prepare_study_databases.py"),
                "--blueprint-db", str(blueprint),
                "--output-root", str(root / "inputs"),
                "--study-manifest-out", str(manifest),
                "--design-truth-out", str(design),
                "--abundance-truth-out", str(abundance),
                "--control-runs", "2",
                "--treatment-runs", "2",
                "--study-seed", "12345",
                "--sample-seed-base", "20000",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        return manifest.read_text(), design.read_text(), abundance.read_text()

    first = run_once(tmp_path / "a")
    second = run_once(tmp_path / "b")

    import csv
    manifest_rows = list(csv.DictReader(first[0].splitlines(), delimiter="\t"))
    manifest_rows_2 = list(csv.DictReader(second[0].splitlines(), delimiter="\t"))
    stable_columns = ["RunOrdinal", "RunId", "RunName", "Condition", "Replicate", "TimSimSampleSeed", "AbundanceSeed"]
    assert [[row[column] for column in stable_columns] for row in manifest_rows] == [
        [row[column] for column in stable_columns] for row in manifest_rows_2
    ]
    assert first[1:] == second[1:]
    design_rows = list(csv.DictReader(first[1].splitlines(), delimiter="\t"))
    abundance_rows = list(csv.DictReader(first[2].splitlines(), delimiter="\t"))
    assert [row["RunId"] for row in manifest_rows] == ["C01", "C02", "T01", "T02"]
    assert {row["Condition"] for row in manifest_rows} == {"control", "treatment"}
    assert len(design_rows) == 10
    assert len(abundance_rows) == 30 * 4
    assert all(float(row["ConditionLog2Effect"]) == 0.0 for row in abundance_rows if row["Condition"] == "control")
    assert any(float(row["ConditionLog2Effect"]) != 0.0 for row in abundance_rows if row["Condition"] == "treatment")

    # Molecular coordinates are copied unchanged while abundance changes per run.
    base_coords = None
    event_vectors = []
    for row in manifest_rows:
        db = Path(row["InputDatabase"])
        with sqlite3.connect(db) as con:
            coords = con.execute("SELECT peptide_id, sequence, retention_time FROM peptides ORDER BY peptide_id").fetchall()
            ions = con.execute("SELECT peptide_id, charge, mz, ion_mobility FROM ions ORDER BY peptide_id").fetchall()
            events = tuple(value[0] for value in con.execute("SELECT events FROM peptides ORDER BY peptide_id"))
        if base_coords is None:
            base_coords = (coords, ions)
        assert (coords, ions) == base_coords
        event_vectors.append(events)
    assert len(set(event_vectors)) == 4


def test_render_study_configs_supports_arbitrary_run_manifest(tmp_path: Path) -> None:
    reference = tmp_path / "reference.d"
    reference.mkdir()
    fasta = tmp_path / "synthetic.fasta"
    fasta.write_text(">P\nPEPTIDEK\n", encoding="utf-8")
    manifest = tmp_path / "study.tsv"
    rows = [
        (1, "C01", "OpenSwathTimSim_control_01", "control", 1, 1001),
        (2, "C02", "OpenSwathTimSim_control_02", "control", 2, 1002),
        (3, "T01", "OpenSwathTimSim_treatment_01", "treatment", 1, 1003),
        (4, "T02", "OpenSwathTimSim_treatment_02", "treatment", 2, 1004),
    ]
    import csv
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["RunOrdinal", "RunId", "RunName", "Condition", "Replicate", "TimSimSampleSeed", "AbundanceSeed", "InputDirectory", "InputDatabase"])
        for ordinal, run_id, run_name, condition, replicate, sample_seed in rows:
            input_dir = tmp_path / "inputs" / run_name
            input_dir.mkdir(parents=True)
            (input_dir / "synthetic_data.db").touch()
            writer.writerow([ordinal, run_id, run_name, condition, replicate, sample_seed, sample_seed + 100, input_dir, input_dir / "synthetic_data.db"])

    rendered = tmp_path / "rendered"
    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "render_study_configs.py"),
            "--repo-root", str(ROOT),
            "--output-root", str(tmp_path / "out"),
            "--rendered-dir", str(rendered),
            "--reference", str(reference),
            "--fasta", str(fasta),
            "--study-manifest", str(manifest),
            "--blueprint-sample-seed", "999",
            "--n-proteins", "1",
            "--num-peptides-total", "1",
            "--num-sample-peptides", "1",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    configs = sorted(rendered.glob("*.toml"))
    assert len(configs) == 5
    text = "\n".join(path.read_text(encoding="utf-8") for path in configs)
    assert 'experiment_name = "OpenSwathTimSim_control_01"' in text
    assert 'experiment_name = "OpenSwathTimSim_treatment_02"' in text
    assert "add_real_data_noise = false" in text
    assert "superimpose_on_reference = false" in text


def test_multirun_study_defaults_and_dynamic_opendia_contract() -> None:
    study = (ROOT / "scripts" / "generate_study.sh").read_text(encoding="utf-8")
    assert "CONTROL_RUNS=25" in study
    assert "TREATMENT_RUNS=25" in study
    assert "PRECURSORS=1000" in study
    assert "ENTRAPMENTS=1000" in study
    assert "SIMULATED_PEPTIDES=10000" in study
    assert "FASTA_PEPTIDES=20000" in study
    assert "select_reference_precursors.py" in study
    assert "prepare_study_databases.py" in study
    assert "build_study_realized_truth.py" in study

    opendia = (ROOT / "scripts" / "run_opendia.sh").read_text(encoding="utf-8")
    assert 'EXPECTED_RUNS="$(python3 - "$MANIFEST"' in opendia
    assert "expected exactly three fixture" not in opendia

    entrapment = (TOOLS / "benchmark_entrapment.py").read_text(encoding="utf-8")
    assert "RUN_ROLES =" not in entrapment
    assert 'run_roles = sorted(str(value) for value in run_best["run_role"].dropna().unique())' in entrapment


def test_benchmark_study_separates_design_realized_and_observed_truth(tmp_path: Path) -> None:
    import csv
    build = tmp_path / "study"
    opendia = tmp_path / "opendia"
    out = tmp_path / "results"
    build.mkdir(); opendia.mkdir()

    run_rows = [
        (1, "C01", "OpenSwathTimSim_control_01", "control", 1),
        (2, "C02", "OpenSwathTimSim_control_02", "control", 2),
        (3, "T01", "OpenSwathTimSim_treatment_01", "treatment", 1),
        (4, "T02", "OpenSwathTimSim_treatment_02", "treatment", 2),
    ]
    with (build / "OpenSwathTimSim.study_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["RunOrdinal", "RunId", "RunName", "Condition", "Replicate"])
        writer.writerows(run_rows)

    with (build / "OpenSwathTimSim.high_signal_selection.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["selection_rank", "precursor_key", "sequence", "protein_id", "charge"], delimiter="\t")
        writer.writeheader()
        writer.writerows([
            {"selection_rank": 1, "precursor_key": "PEPTIDEAK/2", "sequence": "PEPTIDEAK", "protein_id": "P1", "charge": 2},
            {"selection_rank": 2, "precursor_key": "PEPTIDEBK/2", "sequence": "PEPTIDEBK", "protein_id": "P2", "charge": 2},
        ])

    design_rows = [
        {"ProteinId": "P1", "TreatmentClass": "up", "DesignLog2FC": 1.0, "DesignFoldChange": 2.0},
        {"ProteinId": "P2", "TreatmentClass": "unchanged", "DesignLog2FC": 0.0, "DesignFoldChange": 1.0},
    ]
    with (build / "OpenSwathTimSim.study_design_truth.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(design_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(design_rows)

    realized_rows = []
    result_rows = []
    for _, _, run_name, condition, replicate in run_rows:
        for sequence, protein, design_fc, control_events, treatment_events, base_intensity in [
            ("PEPTIDEAK", "P1", 1.0, 100.0, 200.0, 1000.0),
            ("PEPTIDEBK", "P2", 0.0, 100.0, 100.0, 800.0),
        ]:
            events = treatment_events if condition == "treatment" else control_events
            intensity = base_intensity * (2.0 if sequence == "PEPTIDEAK" and condition == "treatment" else 1.0)
            realized_rows.append({
                "selection_rank": 1 if sequence == "PEPTIDEAK" else 2,
                "precursor_key": f"{sequence}/2",
                "TransitionGroupId": f"TIMSIM_{sequence}_2",
                "PeptideSequence": sequence,
                "PrecursorCharge": 2,
                "ProteinId": protein,
                "RunId": run_name,
                "RunName": run_name,
                "Condition": condition,
                "Replicate": replicate,
                "TreatmentClass": "up" if protein == "P1" else "unchanged",
                "DesignLog2FC": design_fc,
                "BaselineEvents": 100.0,
                "RealizedInputEvents": events,
                "TotalLog2FactorVsBlueprint": 1.0 if events == 200 else 0.0,
                "PrecursorMz": 500.0,
                "AssayRT": 10.0,
                "AssayIM": 1.0,
                "PeptideEvents": events,
                "IonRelativeAbundance": 1.0,
                "FrameAbundanceSum": 1.0,
                "ScanAbundanceSum": 1.0,
                "RealizedEventProxy": events,
                "RealizedRTApex": 10.0,
                "RealizedRTCentroid": 10.0,
                "RealizedIMApexSQLite": 1.0,
                "RealizedIMCentroidSQLite": 1.0,
                "FramePoints": 1,
                "ScanPoints": 1,
            })
            result_rows.append({
                "run_name": run_name,
                "Sequence": sequence,
                "Charge": 2,
                "RT": 10.0,
                "EXP_IM": 1.0,
                "Intensity": intensity,
                "decoy": 0,
                "m_score": 0.001,
                "d_score": 5.0,
            })
    with (build / "OpenSwathTimSim.realized_truth.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(realized_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(realized_rows)
    with (opendia / "OpenDIA.results.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(result_rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(result_rows)

    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "benchmark_study.py"),
            "--build-dir", str(build),
            "--opendia-dir", str(opendia),
            "--out-dir", str(out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    report = (out / "study_benchmark_report.md").read_text(encoding="utf-8")
    effects = (out / "peptide_effects.tsv").read_text(encoding="utf-8")
    assert "design vs realized" in report
    assert "realized vs OpenDIA" in report
    assert "PEPTIDEAK" in effects
    assert "1.0" in effects


def test_study_benchmark_runner_uses_project_python_and_supports_postprocess_only() -> None:
    runner = (ROOT / "scripts" / "run_study_benchmark_lanes.sh").read_text(encoding="utf-8")
    assert '"$ROOT/scripts/benchmark_study_lanes.sh" "$BUILD_DIR" "$BENCH_ROOT"' in runner
    assert "\npython " not in runner

    helper = (ROOT / "scripts" / "benchmark_study_lanes.sh").read_text(encoding="utf-8")
    assert 'PYTHON="$FIXTURE_VENV/bin/python"' in helper
    assert '"$PYTHON" "$ROOT/tools/benchmark_study.py"' in helper
    assert '"$PYTHON" "$ROOT/tools/benchmark_entrapment.py"' in helper


def test_quantification_diagnostics_reports_heavy_tail_and_run_signal_fidelity(tmp_path: Path) -> None:
    import pandas as pd

    study = tmp_path / "study"
    study.mkdir()
    runs = [
        ("C01", "control"), ("C02", "control"),
        ("T01", "treatment"), ("T02", "treatment"),
    ]
    peptide_rows = []
    measurement_rows = []
    for idx, (sequence, realized_fc, observed_fc, protein) in enumerate([
        ("PEPAK", 1.0, 1.0, "P1"),
        ("PEPBK", 0.0, 0.0, "P2"),
        ("PEPCK", -1.0, 0.5, "P3"),
    ]):
        peptide_rows.append({
            "PeptideSequence": sequence,
            "PrecursorCharge": 2,
            "ProteinId": protein,
            "TreatmentClass": "unchanged" if realized_fc == 0 else ("up" if realized_fc > 0 else "down"),
            "DesignLog2FC": realized_fc,
            "RealizedConditionLog2FC": realized_fc,
            "ObservedLog2FC": observed_fc,
            "ControlDetected": 2,
            "TreatmentDetected": 2,
            "ControlRuns": 2,
            "TreatmentRuns": 2,
            "WelchPValue": 0.1,
            "BH_QValue": 0.2,
            "CalledDifferential": False,
            "TruthDifferential": realized_fc != 0,
        })
        for run_name, condition in runs:
            truth = 100.0 * (2.0 ** realized_fc if condition == "treatment" else 1.0)
            observed = 1000.0 * (2.0 ** observed_fc if condition == "treatment" else 1.0)
            measurement_rows.append({
                "run_name": run_name,
                "condition": condition,
                "sequence": sequence,
                "charge": 2,
                "realized_input_events": truth * (idx + 1),
                "intensity": observed * (idx + 1),
                "rt_error_seconds": 0.1 + idx * 0.1,
                "im_error": 0.001 + idx * 0.001,
            })
    pd.DataFrame(peptide_rows).to_csv(study / "peptide_effects.tsv", sep="\t", index=False)
    pd.DataFrame(measurement_rows).to_csv(study / "precursor_measurements.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"ProteinId": "P1", "TreatmentClass": "up", "DesignLog2FC": 1.0, "RealizedConditionLog2FC": 1.0, "ObservedLog2FC": 1.0, "WelchPValue": 0.1, "ControlRuns": 2, "TreatmentRuns": 2, "SelectedPrecursors": 1, "BH_QValue": 0.2, "CalledDifferential": False, "TruthDifferential": True},
        {"ProteinId": "P2", "TreatmentClass": "unchanged", "DesignLog2FC": 0.0, "RealizedConditionLog2FC": 0.0, "ObservedLog2FC": 0.0, "WelchPValue": 0.1, "ControlRuns": 2, "TreatmentRuns": 2, "SelectedPrecursors": 2, "BH_QValue": 0.2, "CalledDifferential": False, "TruthDifferential": False},
        {"ProteinId": "P3", "TreatmentClass": "down", "DesignLog2FC": -1.0, "RealizedConditionLog2FC": -1.0, "ObservedLog2FC": 0.5, "WelchPValue": 0.1, "ControlRuns": 2, "TreatmentRuns": 2, "SelectedPrecursors": 3, "BH_QValue": 0.2, "CalledDifferential": False, "TruthDifferential": True},
    ]).to_csv(study / "protein_effects.tsv", sep="\t", index=False)

    subprocess.run(
        [sys.executable, str(TOOLS / "diagnose_study_quantification.py"), "--study-results-dir", str(study)],
        check=True,
        capture_output=True,
        text=True,
    )
    diagnostic = study / "quantification_diagnostics"
    assert (diagnostic / "quantification_diagnostics.md").is_file()
    assert (diagnostic / "run_signal_fidelity.tsv").is_file()
    outliers = pd.read_csv(diagnostic / "peptide_quant_outliers.tsv", sep="\t")
    assert outliers.iloc[0]["PeptideSequence"] == "PEPCK"
    summary = __import__("json").loads((diagnostic / "quantification_diagnostics.json").read_text())
    assert summary["peptide_residual"]["gt_1_0"] == 1


def test_raw_oracle_transition_integration_is_truth_centered() -> None:
    import numpy as np

    oracle = load_tool("extract_raw_signal_oracle")
    mz = np.array([499.9900, 500.0000, 500.0040, 600.0000], dtype=float)
    intensity = np.array([1.0, 10.0, 20.0, 30.0], dtype=float)
    mobility = np.array([0.90, 1.00, 1.02, 1.00], dtype=float)
    order = np.argsort(mz)
    signal, peaks = oracle.integrate_transition(
        mz[order],
        intensity[order],
        mobility[order],
        product_mz=500.0,
        target_im=1.0,
        ppm=20.0,
        im_half_window=0.03,
    )
    assert signal == 30.0
    assert peaks == 2


def test_raw_oracle_frame_adapter_accepts_dataframe_frames() -> None:
    import pandas as pd
    import numpy as np

    oracle = load_tool("extract_raw_signal_oracle")

    class FakeFrame:
        def df(self):
            return pd.DataFrame(
                {
                    "mz": [400.0, 500.0],
                    "intensity": [11.0, 22.0],
                    "inv_ion_mobility": [0.9, 1.1],
                }
            )

    mz, intensity, mobility = oracle.frame_arrays(object(), FakeFrame(), 1)
    assert np.allclose(mz, [400.0, 500.0])
    assert np.allclose(intensity, [11.0, 22.0])
    assert np.allclose(mobility, [0.9, 1.1])


def test_raw_oracle_benchmark_separates_raw_and_opendia_response(tmp_path: Path) -> None:
    import pandas as pd

    oracle_dir = tmp_path / "oracle"
    study_dir = tmp_path / "study"
    oracle_dir.mkdir(); study_dir.mkdir()
    rows = []
    observed_rows = []
    peptide_rows = []
    for idx, (sequence, realized_fc, raw_fc, observed_fc) in enumerate([
        ("PEPAK", 1.0, 1.0, 0.8),
        ("PEPBK", 0.0, 0.0, 0.0),
        ("PEPCK", -1.0, -1.0, -0.8),
    ]):
        for run_name, condition in [("C01", "control"), ("C02", "control"), ("T01", "treatment"), ("T02", "treatment")]:
            truth = 100.0 * (2.0 ** realized_fc if condition == "treatment" else 1.0) * (idx + 1)
            raw = 1000.0 * (2.0 ** raw_fc if condition == "treatment" else 1.0) * (idx + 1)
            observed = 5000.0 * (2.0 ** observed_fc if condition == "treatment" else 1.0) * (idx + 1)
            rows.append({
                "RunName": run_name,
                "PeptideSequence": sequence,
                "PrecursorCharge": 2,
                "ProteinId": f"P{idx+1}",
                "TransitionGroupId": f"TIMSIM_{sequence}_2",
                "PrecursorMz": 500.0 + idx,
                "RealizedRTApex": 10.0,
                "RealizedIMApex": 1.0,
                "RealizedInputEvents": truth,
                "RealizedEventProxy": truth,
                "RawOracleIntensity": raw,
                "MatchedTransitions": 8,
                "TotalTransitions": 8,
                "MatchedRawPeaks": 16,
                "FramesConsidered": 3,
                "WindowGroup": 1,
            })
            observed_rows.append({
                "run_name": run_name,
                "condition": condition,
                "sequence": sequence,
                "charge": 2,
                "intensity": observed,
                "rt_error_seconds": 0.1,
                "im_error": 0.001,
            })
        peptide_rows.append({
            "PeptideSequence": sequence,
            "PrecursorCharge": 2,
            "ProteinId": f"P{idx+1}",
            "TreatmentClass": "up" if realized_fc > 0 else ("down" if realized_fc < 0 else "unchanged"),
            "DesignLog2FC": realized_fc,
            "RealizedConditionLog2FC": realized_fc,
            "ObservedLog2FC": observed_fc,
            "ControlDetected": 2,
            "TreatmentDetected": 2,
            "ControlRuns": 2,
            "TreatmentRuns": 2,
            "WelchPValue": 0.01,
            "BH_QValue": 0.02,
            "CalledDifferential": realized_fc != 0,
            "TruthDifferential": realized_fc != 0,
        })
    pd.DataFrame(rows).to_csv(oracle_dir / "raw_signal_oracle_measurements.tsv", sep="\t", index=False)
    pd.DataFrame(observed_rows).to_csv(study_dir / "precursor_measurements.tsv", sep="\t", index=False)
    pd.DataFrame(peptide_rows).to_csv(study_dir / "peptide_effects.tsv", sep="\t", index=False)

    subprocess.run(
        [
            sys.executable,
            str(TOOLS / "benchmark_raw_signal_oracle.py"),
            "--oracle-dir", str(oracle_dir),
            "--study-results-dir", str(study_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    summary = __import__("json").loads((oracle_dir / "raw_signal_oracle_summary.json").read_text())
    assert abs(summary["condition_effects"]["realized_vs_raw"]["slope"] - 1.0) < 1e-9
    assert abs(summary["condition_effects"]["raw_vs_opendia"]["slope"] - 0.8) < 1e-9
    report = (oracle_dir / "raw_signal_oracle_report.md").read_text(encoding="utf-8")
    assert "Realized input effect → raw oracle effect" in report
    assert "Raw oracle effect → OpenDIA effect" in report


def test_raw_oracle_runner_is_postprocess_only_and_uses_project_python() -> None:
    runner = (ROOT / "scripts" / "run_raw_signal_oracle.sh").read_text(encoding="utf-8")
    assert 'PYTHON="$FIXTURE_VENV/bin/python"' in runner
    assert "extract_raw_signal_oracle.py" in runner
    assert "benchmark_raw_signal_oracle.py" in runner
    assert "run_opendia.sh" not in runner
    assert "generate_study.sh" not in runner


def test_raw_oracle_sweep_grid_has_expected_144_configurations() -> None:
    sweep = load_tool("sweep_raw_signal_oracle")
    rows = sweep.configuration_rows(
        sweep.DEFAULT_RT_WINDOWS,
        sweep.DEFAULT_IM_WINDOWS,
        sweep.DEFAULT_PPM_WINDOWS,
        sweep.DEFAULT_SUBSETS,
    )
    assert len(rows) == 4 * 3 * 3 * 4 == 144
    assert any(row["ConfigId"] == "rt6_im0.03_ppm25_top8" for row in rows)
    assert any(row["ConfigId"] == "rt1_im0.01_ppm10_top4" for row in rows)


def test_raw_oracle_sweep_transition_grid_is_nested() -> None:
    import numpy as np

    sweep = load_tool("sweep_raw_signal_oracle")
    mz = np.array([499.996, 500.000, 500.008, 500.020], dtype=float)
    intensity = np.array([5.0, 10.0, 20.0, 40.0], dtype=float)
    mobility = np.array([1.000, 1.005, 1.018, 1.040], dtype=float)
    grid = sweep.integrate_transition_grid(
        mz,
        intensity,
        mobility,
        product_mz=500.0,
        target_im=1.0,
        ppm_windows=(10.0, 15.0, 25.0),
        im_windows=(0.01, 0.02, 0.03),
    )
    assert grid.shape == (3, 3)
    # Every wider mass/mobility window must contain at least as much signal.
    assert np.all(np.diff(grid, axis=0) >= 0)
    assert np.all(np.diff(grid, axis=1) >= 0)
    assert grid[0, 0] == 15.0
    assert grid[2, 1] == 35.0


def test_raw_oracle_specificity_collision_metric_excludes_same_precursor() -> None:
    import pandas as pd

    sweep = load_tool("sweep_raw_signal_oracle")
    background = pd.DataFrame(
        [
            {"PrecursorKey": "TARGET/2", "RT": 100.0, "IM": 1.0, "ProductMz": 500.0, "PredictedIntensity": 100.0},
            {"PrecursorKey": "COMP1/2", "RT": 101.0, "IM": 1.005, "ProductMz": 500.004, "PredictedIntensity": 20.0},
            {"PrecursorKey": "COMP2/2", "RT": 120.0, "IM": 1.005, "ProductMz": 500.003, "PredictedIntensity": 30.0},
        ]
    )
    count, competitors, intensity, nearest_ppm = sweep.collision_metrics_for_transition(
        500.0,
        "TARGET/2",
        100.0,
        1.0,
        1,
        background,
        max_rt=6.0,
        max_im=0.03,
        max_ppm=25.0,
    )
    assert count == 1
    assert competitors == 1
    assert intensity == 20.0
    assert 7.9 < nearest_ppm < 8.1


def test_raw_oracle_sweep_runner_is_postprocess_only_and_uses_project_python() -> None:
    runner = (ROOT / "scripts" / "run_raw_signal_oracle_sweep.sh").read_text(encoding="utf-8")
    assert 'PYTHON="$FIXTURE_VENV/bin/python"' in runner
    assert "sweep_raw_signal_oracle.py" in runner
    assert "run_opendia.sh" not in runner
    assert "generate_study.sh" not in runner
    assert "1,2,4,6" in runner
    assert "0.01,0.02,0.03" in runner
    assert "10,15,25" in runner
    assert "8,6,4,3" in runner
