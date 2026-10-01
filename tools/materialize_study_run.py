#!/usr/bin/env python3
"""Materialize one planned TimSim run from an immutable blueprint database."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import random
import shutil
import sqlite3
from pathlib import Path
from typing import Iterable

import numpy as np


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [row[1] for row in connection.execute(f"PRAGMA table_info({quote_identifier(table)})")]


def first_present(options: Iterable[str], available: Iterable[str], label: str) -> str:
    available_set = set(available)
    for option in options:
        if option in available_set:
            return option
    raise RuntimeError(f"Could not identify {label}; available columns: {sorted(available_set)}")


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="", encoding="utf-8")
    return path.open(mode, newline="", encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "r") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_run(manifest: Path, ordinal: int) -> dict[str, str]:
    rows = read_tsv(manifest)
    matches = [row for row in rows if int(row["RunOrdinal"]) == ordinal]
    if len(matches) != 1:
        raise SystemExit(f"Expected exactly one run with ordinal {ordinal}; found {len(matches)}")
    return matches[0]


def load_design(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    return {row["ProteinId"]: row for row in rows}


def load_selected_sequences(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    selected: dict[str, dict[str, str]] = {}
    for row in rows:
        sequence = row["sequence"].strip().upper()
        if sequence in selected:
            raise SystemExit(f"Selected precursor table contains duplicate peptide sequence: {sequence}")
        selected[sequence] = row
    return selected


def load_peptides(connection: sqlite3.Connection) -> tuple[list[tuple[int, str, str, float]], dict[str, str]]:
    available = columns(connection, "peptides")
    id_col = first_present(("peptide_id", "id"), available, "peptide ID")
    sequence_col = first_present(("sequence", "peptide"), available, "peptide sequence")
    protein_col = first_present(("protein", "protein_id", "protein_name"), available, "protein")
    events_col = first_present(("events", "total_events", "abundance"), available, "peptide abundance")
    query = (
        f"SELECT {quote_identifier(id_col)}, {quote_identifier(sequence_col)}, "
        f"{quote_identifier(protein_col)}, {quote_identifier(events_col)} "
        f"FROM {quote_identifier('peptides')}"
    )
    rows = [
        (int(pid), str(sequence).strip().upper(), str(protein), float(events))
        for pid, sequence, protein, events in connection.execute(query)
    ]
    rows.sort(key=lambda row: (row[2], row[1], row[0]))
    return rows, {"id": id_col, "events": events_col}


def update_database(path: Path, columns_map: dict[str, str], values: list[tuple[int, int]]) -> None:
    with sqlite3.connect(path) as connection:
        available = set(columns(connection, "peptides"))
        abundance_columns = [columns_map["events"]]
        for optional in ("events", "total_events", "abundance"):
            if optional in available and optional not in abundance_columns:
                abundance_columns.append(optional)
        assignments = ", ".join(f"{quote_identifier(column)} = ?" for column in abundance_columns)
        query = (
            f"UPDATE {quote_identifier('peptides')} SET {assignments} "
            f"WHERE {quote_identifier(columns_map['id'])} = ?"
        )
        payload = [tuple([events] * len(abundance_columns) + [peptide_id]) for peptide_id, events in values]
        connection.executemany(query, payload)
        connection.commit()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blueprint-db", type=Path, required=True)
    parser.add_argument("--study-manifest", type=Path, required=True)
    parser.add_argument("--design-truth", type=Path, required=True)
    parser.add_argument("--selected-precursors", type=Path, required=True)
    parser.add_argument("--run-ordinal", type=int, required=True)
    parser.add_argument("--output-db", type=Path, required=True)
    parser.add_argument("--selected-input-truth-out", type=Path, required=True)
    parser.add_argument("--stats-out", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.blueprint_db, args.study_manifest, args.design_truth, args.selected_precursors):
        if not path.is_file():
            raise SystemExit(f"Required input does not exist: {path}")

    run = load_run(args.study_manifest, args.run_ordinal)
    design = load_design(args.design_truth)
    selected = load_selected_sequences(args.selected_precursors)

    args.output_db.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.blueprint_db, args.output_db)
    with sqlite3.connect(args.blueprint_db) as connection:
        peptides, columns_map = load_peptides(connection)

    abundance_seed = int(run["AbundanceSeed"])
    rng = random.Random(abundance_seed)
    np_rng = np.random.default_rng(abundance_seed + 104729)
    run_global = rng.gauss(0.0, float(run["RunLog2SD"]))
    cell_size = rng.gauss(0.0, float(run["CellSizeLog2SD"]))
    input_offset = float(run["GlobalInputLog2Offset"])
    profile = run["AbundanceProfile"]
    event_sampling = run["EventSampling"]
    allow_zero = str(run["AllowZeroEvents"]).strip().lower() in {"1", "true", "yes"}

    proteins = sorted({row[2] for row in peptides})
    protein_noise = {protein: rng.gauss(0.0, float(run["ProteinLog2SD"])) for protein in proteins}
    missing_design = sorted(set(proteins) - set(design))
    if missing_design:
        raise SystemExit(f"Design truth is missing {len(missing_design)} blueprint proteins")

    updates: list[tuple[int, int]] = []
    selected_truth: list[dict[str, object]] = []
    positive_peptides = 0
    zero_peptides = 0
    total_baseline = 0
    total_mean = 0.0
    total_realized = 0

    for peptide_id, sequence, protein, events in peptides:
        baseline_events = max(1, int(round(events)))
        peptide_noise = rng.gauss(0.0, float(run["PeptideLog2SD"]))
        condition_effect = float(design[protein]["DesignLog2FC"]) if run["Condition"] == "treatment" else 0.0
        total_log2 = input_offset + cell_size + run_global + protein_noise[protein] + peptide_noise + condition_effect
        mean_events = max(0.0, baseline_events * (2.0 ** total_log2))

        if event_sampling == "poisson":
            realized_events = int(np_rng.poisson(mean_events))
        elif event_sampling == "deterministic":
            realized_events = int(round(mean_events))
        else:
            raise SystemExit(f"Unsupported event sampling mode: {event_sampling}")
        if not allow_zero:
            realized_events = max(1, realized_events)

        updates.append((peptide_id, realized_events))
        total_baseline += baseline_events
        total_mean += mean_events
        total_realized += realized_events
        if realized_events > 0:
            positive_peptides += 1
        else:
            zero_peptides += 1

        selected_row = selected.get(sequence)
        if selected_row is not None:
            selected_truth.append({
                "selection_rank": selected_row["selection_rank"],
                "precursor_key": selected_row["precursor_key"],
                "PeptideSequence": sequence,
                "PrecursorCharge": selected_row["charge"],
                "ProteinId": protein,
                "RunOrdinal": run["RunOrdinal"],
                "RunId": run["RunId"],
                "RunName": run["RunName"],
                "Condition": run["Condition"],
                "Replicate": run["Replicate"],
                "InputLevel": run["InputLevel"],
                "AbundanceProfile": profile,
                "TreatmentClass": design[protein]["TreatmentClass"],
                "DesignLog2FC": design[protein]["DesignLog2FC"],
                "BaselineEvents": baseline_events,
                "GlobalInputLog2Offset": f"{input_offset:.8f}",
                "CellSizeRandomLog2Effect": f"{cell_size:.8f}",
                "RunGlobalLog2Effect": f"{run_global:.8f}",
                "ProteinRandomLog2Effect": f"{protein_noise[protein]:.8f}",
                "PeptideRandomLog2Effect": f"{peptide_noise:.8f}",
                "ConditionLog2Effect": f"{condition_effect:.8f}",
                "TotalLog2FactorVsBlueprint": f"{total_log2:.8f}",
                "ExpectedInputEvents": f"{mean_events:.8f}",
                "RealizedInputEvents": realized_events,
                "BiologicalPresentInRun": int(realized_events > 0),
            })

    update_database(args.output_db, columns_map, updates)
    selected_truth.sort(key=lambda row: int(row["selection_rank"]))
    if len(selected_truth) != len(selected):
        raise SystemExit(
            f"Only {len(selected_truth)} of {len(selected)} selected peptide sequences were found in the blueprint DB"
        )
    write_tsv(args.selected_input_truth_out, selected_truth)

    stats = {
        "run_ordinal": int(run["RunOrdinal"]),
        "run_name": run["RunName"],
        "condition": run["Condition"],
        "replicate": int(run["Replicate"]),
        "profile": profile,
        "input_level": run["InputLevel"],
        "input_log2_offset": input_offset,
        "cell_size_random_log2_effect": cell_size,
        "run_global_log2_effect": run_global,
        "event_sampling": event_sampling,
        "allow_zero_events": allow_zero,
        "blueprint_peptides": len(peptides),
        "positive_peptides": positive_peptides,
        "zero_event_peptides": zero_peptides,
        "positive_fraction": positive_peptides / max(len(peptides), 1),
        "total_baseline_events": total_baseline,
        "total_expected_events": total_mean,
        "total_realized_events": total_realized,
        "selected_targets": len(selected_truth),
        "selected_present_targets": sum(int(row["BiologicalPresentInRun"]) for row in selected_truth),
    }
    args.stats_out.parent.mkdir(parents=True, exist_ok=True)
    args.stats_out.write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")

    print(
        f"Materialized {run['RunName']}: profile={profile}, peptides={len(peptides)}, "
        f"positive={positive_peptides}, selected_present={stats['selected_present_targets']}/{len(selected_truth)}"
    )
    print(f"Source DB: {args.output_db}")
    print(f"Selected input truth: {args.selected_input_truth_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
