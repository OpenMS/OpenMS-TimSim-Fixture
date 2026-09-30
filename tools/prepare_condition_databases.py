#!/usr/bin/env python3
"""Create comparable replicate/treatment inputs from a TimSim baseline DB.

Run 2 is a control replicate with modest peptide-level abundance variation.
Run 3 adds coherent protein-level treatment effects plus modest sample noise.
RT, m/z, charge states, ion mobility, and peptide identities are left unchanged.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import shutil
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


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


@dataclass(frozen=True)
class PeptideRow:
    peptide_id: int
    sequence: str
    protein: str
    events: float


def load_peptides(connection: sqlite3.Connection) -> tuple[list[PeptideRow], dict[str, str]]:
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
        PeptideRow(int(pid), str(seq), str(protein), float(events))
        for pid, seq, protein, events in connection.execute(query)
    ]
    rows.sort(key=lambda row: (row.sequence, row.peptide_id))
    return rows, {
        "id": id_col,
        "sequence": sequence_col,
        "protein": protein_col,
        "events": events_col,
    }


def assign_treatment_classes(
    proteins: list[str],
    up_fraction: float,
    down_fraction: float,
    seed: int,
) -> dict[str, str]:
    shuffled = list(sorted(set(proteins)))
    random.Random(seed).shuffle(shuffled)
    count = len(shuffled)
    up_count = int(round(count * up_fraction))
    down_count = int(round(count * down_fraction))
    if count >= 3 and up_fraction > 0:
        up_count = max(1, up_count)
    if count >= 3 and down_fraction > 0:
        down_count = max(1, down_count)
    if up_count + down_count > count:
        raise ValueError("Treatment up/down fractions select more proteins than exist")

    classes = {protein: "unchanged" for protein in shuffled}
    for protein in shuffled[:up_count]:
        classes[protein] = "up"
    for protein in shuffled[up_count:up_count + down_count]:
        classes[protein] = "down"
    return classes


def update_database(
    path: Path,
    columns_map: dict[str, str],
    new_events: dict[int, int],
) -> None:
    with sqlite3.connect(path) as connection:
        available = set(columns(connection, "peptides"))
        abundance_columns = [columns_map["events"]]
        # Keep TimSim's parallel abundance column consistent when present.
        for optional in ("events", "total_events", "abundance"):
            if optional in available and optional not in abundance_columns:
                abundance_columns.append(optional)

        assignments = ", ".join(f"{quote_identifier(column)} = ?" for column in abundance_columns)
        query = (
            f"UPDATE {quote_identifier('peptides')} SET {assignments} "
            f"WHERE {quote_identifier(columns_map['id'])} = ?"
        )
        payload = [
            tuple([value] * len(abundance_columns) + [peptide_id])
            for peptide_id, value in new_events.items()
        ]
        connection.executemany(query, payload)
        connection.commit()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-db", type=Path, required=True)
    parser.add_argument("--replicate-db", type=Path, required=True)
    parser.add_argument("--treatment-db", type=Path, required=True)
    parser.add_argument("--truth-out", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--replicate-log2-sd", type=float, default=0.10)
    parser.add_argument("--treatment-log2-sd", type=float, default=0.08)
    parser.add_argument("--treatment-up-fraction", type=float, default=0.15)
    parser.add_argument("--treatment-down-fraction", type=float, default=0.15)
    parser.add_argument("--treatment-up-fold", type=float, default=2.0)
    parser.add_argument("--treatment-down-fold", type=float, default=0.5)
    args = parser.parse_args()

    if not args.baseline_db.is_file():
        raise SystemExit(f"Baseline DB does not exist: {args.baseline_db}")
    for value, label in (
        (args.replicate_log2_sd, "replicate-log2-sd"),
        (args.treatment_log2_sd, "treatment-log2-sd"),
    ):
        if value < 0:
            parser.error(f"--{label} must be non-negative")
    if not 0 <= args.treatment_up_fraction <= 1:
        parser.error("--treatment-up-fraction must be between 0 and 1")
    if not 0 <= args.treatment_down_fraction <= 1:
        parser.error("--treatment-down-fraction must be between 0 and 1")
    if args.treatment_up_fraction + args.treatment_down_fraction > 1:
        parser.error("treatment up/down fractions must sum to <= 1")
    if args.treatment_up_fold <= 0 or args.treatment_down_fold <= 0:
        parser.error("treatment fold changes must be positive")

    with sqlite3.connect(args.baseline_db) as connection:
        peptide_rows, columns_map = load_peptides(connection)

    if not peptide_rows:
        raise SystemExit("Baseline peptides table is empty")

    args.replicate_db.parent.mkdir(parents=True, exist_ok=True)
    args.treatment_db.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.baseline_db, args.replicate_db)
    shutil.copy2(args.baseline_db, args.treatment_db)

    treatment_classes = assign_treatment_classes(
        [row.protein for row in peptide_rows],
        args.treatment_up_fraction,
        args.treatment_down_fraction,
        args.seed + 1,
    )

    replicate_rng = random.Random(args.seed + 2)
    treatment_rng = random.Random(args.seed + 3)
    replicate_events: dict[int, int] = {}
    treatment_events: dict[int, int] = {}
    truth_rows: list[dict[str, object]] = []

    for row in peptide_rows:
        replicate_log2 = replicate_rng.gauss(0.0, args.replicate_log2_sd)
        replicate_factor = 2.0 ** replicate_log2

        treatment_class = treatment_classes[row.protein]
        biological_factor = {
            "up": args.treatment_up_fold,
            "down": args.treatment_down_fold,
            "unchanged": 1.0,
        }[treatment_class]
        treatment_noise_log2 = treatment_rng.gauss(0.0, args.treatment_log2_sd)
        treatment_factor = biological_factor * (2.0 ** treatment_noise_log2)

        baseline_events = max(1, int(round(row.events)))
        replicate_value = max(1, int(round(baseline_events * replicate_factor)))
        treatment_value = max(1, int(round(baseline_events * treatment_factor)))
        replicate_events[row.peptide_id] = replicate_value
        treatment_events[row.peptide_id] = treatment_value

        truth_rows.append({
            "PeptideId": row.peptide_id,
            "PeptideSequence": row.sequence,
            "ProteinId": row.protein,
            "TreatmentClass": treatment_class,
            "BaselineEvents": baseline_events,
            "ControlReplicateFactor": f"{replicate_factor:.8f}",
            "ControlReplicateEvents": replicate_value,
            "TreatmentBiologicalFactor": f"{biological_factor:.8f}",
            "TreatmentNoiseFactor": f"{2.0 ** treatment_noise_log2:.8f}",
            "TreatmentTotalFactor": f"{treatment_factor:.8f}",
            "TreatmentEvents": treatment_value,
        })

    update_database(args.replicate_db, columns_map, replicate_events)
    update_database(args.treatment_db, columns_map, treatment_events)

    args.truth_out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(truth_rows[0])
    with args.truth_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(truth_rows)

    counts = {name: list(treatment_classes.values()).count(name) for name in ("up", "down", "unchanged")}
    print(f"Wrote control-replicate source DB: {args.replicate_db}")
    print(f"Wrote treatment source DB:       {args.treatment_db}")
    print(f"Wrote condition truth:           {args.truth_out}")
    print(
        "Treatment proteins: "
        f"{counts['up']} up, {counts['down']} down, {counts['unchanged']} unchanged"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
