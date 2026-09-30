#!/usr/bin/env python3
"""Validate the structural and truth contract of a generated multi-run study."""
from __future__ import annotations

import argparse
import json
import sqlite3
from collections import Counter
from pathlib import Path

import pandas as pd


def validate_protein_selection_contract(selected: pd.DataFrame, manifest: dict) -> None:
    selection = manifest.get("selection", {})
    if selection.get("mode") != "protein_balanced":
        return
    target_proteins = int(selection.get("target_proteins", 0))
    precursors_per_protein = int(selection.get("precursors_per_protein", 0))
    expected_targets = target_proteins * precursors_per_protein
    if target_proteins < 1 or precursors_per_protein < 1:
        raise SystemExit("Protein-balanced selection manifest is missing positive target_proteins/precursors_per_protein")
    if int(manifest.get("requested_precursors", 0)) != expected_targets:
        raise SystemExit(
            "Protein-balanced manifest is inconsistent: requested_precursors != "
            "target_proteins * precursors_per_protein"
        )
    counts = selected.groupby("protein_id")["precursor_key"].nunique()
    if len(counts) != target_proteins:
        raise SystemExit(
            f"Protein-balanced selection contains {len(counts)} proteins; expected {target_proteins}"
        )
    bad = counts[counts.ne(precursors_per_protein)]
    if not bad.empty:
        examples = ", ".join(f"{protein}:{count}" for protein, count in bad.head(10).items())
        raise SystemExit(
            f"Protein-balanced selection requires exactly {precursors_per_protein} precursors/protein; "
            f"violations: {examples}"
        )
    required_specificity = {
        "specificity_collision_free_fragments",
        "specificity_all_competitor_precursors_sum",
        "specificity_all_collision_count_sum",
        "specificity_top_fragments",
        "specificity_top_competitor_precursors_sum",
        "specificity_top_collision_count_sum",
        "protein_selection_rank",
        "protein_precursor_rank",
    }
    missing = required_specificity - set(selected.columns)
    if missing:
        raise SystemExit(
            "Protein-balanced selection is missing simulator-only specificity/protein-rank fields: "
            + ", ".join(sorted(missing))
        )
    if selected["protein_selection_rank"].nunique() != target_proteins:
        raise SystemExit("protein_selection_rank is not unique per selected protein")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study_root", type=Path)
    args = parser.parse_args()
    root = args.study_root.resolve()
    manifest_path = root / "fixture_manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("fixture_kind") != "openms_timsim_multirun_study_v1":
        raise SystemExit(f"Unexpected fixture_kind: {manifest.get('fixture_kind')}")

    run_table = pd.read_csv(root / "OpenSwathTimSim.study_manifest.tsv", sep="\t")
    selected = pd.read_csv(root / "OpenSwathTimSim.high_signal_selection.tsv", sep="\t")
    realized = pd.read_csv(root / "OpenSwathTimSim.realized_truth.tsv", sep="\t")
    design = pd.read_csv(root / "OpenSwathTimSim.study_design_truth.tsv", sep="\t")
    abundance = pd.read_csv(root / "OpenSwathTimSim.study_abundance_truth.tsv", sep="\t")
    target_library = pd.read_csv(root / "OpenSwathTimSim.target_only.transitions.tsv", sep="\t")
    combined_library = pd.read_csv(root / "OpenSwathTimSim.transitions.tsv", sep="\t")
    entrapment_truth = pd.read_csv(root / "OpenSwathTimSim.entrapment_truth.tsv", sep="\t")

    expected_runs = int(manifest["study"]["control_runs"]) + int(manifest["study"]["treatment_runs"])
    expected_targets = int(manifest["requested_precursors"])
    if len(run_table) != expected_runs:
        raise SystemExit(f"Study manifest has {len(run_table)} runs; expected {expected_runs}")
    if run_table["RunName"].nunique() != expected_runs or run_table["RunId"].nunique() != expected_runs:
        raise SystemExit("Study RunName/RunId values are not unique")
    counts = run_table["Condition"].value_counts().to_dict()
    if int(counts.get("control", 0)) != int(manifest["study"]["control_runs"]):
        raise SystemExit("Control run count does not match fixture manifest")
    if int(counts.get("treatment", 0)) != int(manifest["study"]["treatment_runs"]):
        raise SystemExit("Treatment run count does not match fixture manifest")
    if len(selected) != expected_targets or selected["precursor_key"].nunique() != expected_targets:
        raise SystemExit(f"Selected precursor table does not contain exactly {expected_targets} unique targets")

    validate_protein_selection_contract(selected, manifest)

    expected_entrapments = int(manifest.get("entrapment", {}).get("precursors", 0))
    transitions_per = int(manifest["transitions_per_precursor"])
    target_counts = Counter(target_library["TransitionGroupId"].astype(str))
    combined_counts = Counter(combined_library["TransitionGroupId"].astype(str))
    if len(target_counts) != expected_targets or any(value != transitions_per for value in target_counts.values()):
        raise SystemExit("Target-only transition library does not match target/transition counts in manifest")
    if len(combined_counts) != expected_targets + expected_entrapments or any(value != transitions_per for value in combined_counts.values()):
        raise SystemExit("Combined transition library does not match target + entrapment counts in manifest")
    if len(entrapment_truth) != expected_entrapments:
        raise SystemExit(f"Entrapment truth has {len(entrapment_truth)} rows; expected {expected_entrapments}")
    if "Decoy" in combined_library and not combined_library["Decoy"].astype(str).str.lower().isin({"0", "false"}).all():
        raise SystemExit("Input target/entrapment library must contain Decoy=0 only")
    ent_sequences = set(entrapment_truth["EntrapmentSequence"].astype(str).str.upper())
    blueprint_db = root / str(manifest["blueprint_run"]) / "synthetic_data.db"
    with sqlite3.connect(blueprint_db) as connection:
        peptide_cols = [row[1] for row in connection.execute('PRAGMA table_info("peptides")')]
        sequence_col = next((name for name in ("sequence", "peptide") if name in peptide_cols), None)
        if sequence_col is None:
            raise SystemExit("Could not locate peptide sequence column in blueprint DB")
        simulated = {str(row[0]).upper() for row in connection.execute(f'SELECT "{sequence_col}" FROM peptides')}
    overlap = ent_sequences & simulated
    if overlap:
        raise SystemExit(f"Entrapment sequences occur in blueprint TimSim population; examples: {sorted(overlap)[:10]}")
    if len(realized) != expected_targets * expected_runs:
        raise SystemExit(
            f"Realized truth has {len(realized)} rows; expected {expected_targets * expected_runs}"
        )
    per_run = realized.groupby("RunName")["precursor_key"].nunique()
    if len(per_run) != expected_runs or not per_run.eq(expected_targets).all():
        raise SystemExit("Each study run must have realized truth for every selected precursor")

    # Assay coordinates must remain invariant for each selected precursor.
    for column in ("PrecursorMz", "AssayRT", "AssayIM"):
        if realized.groupby("precursor_key")[column].nunique(dropna=False).max() != 1:
            raise SystemExit(f"{column} changed across study runs")

    proteins = design["ProteinId"].nunique()
    if proteins == 0:
        raise SystemExit("Study design truth contains no proteins")
    if not set(design["TreatmentClass"]) <= {"up", "down", "unchanged"}:
        raise SystemExit("Study design truth contains an unknown TreatmentClass")
    if not (design.loc[design["TreatmentClass"] == "unchanged", "DesignLog2FC"].abs() < 1e-12).all():
        raise SystemExit("Unchanged proteins must have DesignLog2FC=0")
    if len(abundance) != expected_runs * abundance["PeptideId"].nunique():
        raise SystemExit("Study abundance truth is not rectangular across runs/peptides")

    for row in run_table.itertuples(index=False):
        run_dir = root / row.RunName
        tdf = run_dir / f"{row.RunName}.d"
        db = run_dir / "synthetic_data.db"
        link = root / "tdfs" / f"{row.RunName}.d"
        if not (tdf / "analysis.tdf").is_file() or not (tdf / "analysis.tdf_bin").is_file():
            raise SystemExit(f"Missing generated TDF for {row.RunName}")
        if not db.is_file():
            raise SystemExit(f"Missing synthetic_data.db for {row.RunName}")
        if not link.is_symlink() or link.resolve() != tdf.resolve():
            raise SystemExit(f"Invalid TDF symlink for {row.RunName}: {link}")

    print(
        f"Study contract OK: {expected_runs} runs "
        f"({counts.get('control', 0)} control + {counts.get('treatment', 0)} treatment), "
        f"{expected_targets} frozen targets, {proteins} design proteins"
    )
    selection = manifest.get("selection", {})
    if selection.get("mode") == "protein_balanced":
        print(
            "Production target contract OK: "
            f"{selection['target_proteins']} proteins x "
            f"{selection['precursors_per_protein']} precursors/protein"
        )
    print("Design truth, realized abundance truth, and per-run realized precursor truth are internally consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
