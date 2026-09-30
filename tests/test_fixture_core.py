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
