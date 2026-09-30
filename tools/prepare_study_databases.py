#!/usr/bin/env python3
"""Create condition-specific TimSim source databases for a multi-run study.

The blueprint TimSim database defines the immutable molecular universe and baseline
peptide abundance. This tool creates independent control/treatment replicate source
databases by changing peptide abundance only. Sequence, RT, precursor m/z, charge,
ion mobility and all molecular identities remain unchanged.

Two quantitative truth layers are written:

* study_design_truth.tsv: condition-level protein effects fixed before any run.
* study_abundance_truth.tsv: run/peptide realized TimSim input abundance after
  deterministic biological replicate variation.
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


@dataclass(frozen=True)
class ProteinDesign:
    protein: str
    treatment_class: str
    design_log2fc: float


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
        PeptideRow(int(pid), str(sequence), str(protein), float(events))
        for pid, sequence, protein, events in connection.execute(query)
    ]
    rows.sort(key=lambda row: (row.protein, row.sequence, row.peptide_id))
    return rows, {"id": id_col, "sequence": sequence_col, "protein": protein_col, "events": events_col}


def update_database(path: Path, columns_map: dict[str, str], events_by_peptide: dict[int, int]) -> None:
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
        payload = [
            tuple([value] * len(abundance_columns) + [peptide_id])
            for peptide_id, value in events_by_peptide.items()
        ]
        connection.executemany(query, payload)
        connection.commit()


def assign_protein_design(
    proteins: list[str],
    *,
    up_fraction: float,
    down_fraction: float,
    min_abs_log2fc: float,
    max_abs_log2fc: float,
    seed: int,
) -> dict[str, ProteinDesign]:
    unique = sorted(set(proteins))
    rng = random.Random(seed)
    shuffled = unique.copy()
    rng.shuffle(shuffled)
    up_count = int(round(len(shuffled) * up_fraction))
    down_count = int(round(len(shuffled) * down_fraction))
    if len(shuffled) >= 3 and up_fraction > 0:
        up_count = max(1, up_count)
    if len(shuffled) >= 3 and down_fraction > 0:
        down_count = max(1, down_count)
    if up_count + down_count > len(shuffled):
        raise ValueError("Treatment up/down fractions select more proteins than exist")

    up = set(shuffled[:up_count])
    down = set(shuffled[up_count:up_count + down_count])
    result: dict[str, ProteinDesign] = {}
    for protein in unique:
        if protein in up:
            magnitude = rng.uniform(min_abs_log2fc, max_abs_log2fc)
            result[protein] = ProteinDesign(protein, "up", magnitude)
        elif protein in down:
            magnitude = rng.uniform(min_abs_log2fc, max_abs_log2fc)
            result[protein] = ProteinDesign(protein, "down", -magnitude)
        else:
            result[protein] = ProteinDesign(protein, "unchanged", 0.0)
    return result


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"Refusing to write empty TSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blueprint-db", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--study-manifest-out", type=Path, required=True)
    parser.add_argument("--design-truth-out", type=Path, required=True)
    parser.add_argument("--abundance-truth-out", type=Path, required=True)
    parser.add_argument("--control-runs", type=int, default=25)
    parser.add_argument("--treatment-runs", type=int, default=25)
    parser.add_argument("--study-seed", type=int, default=2026093001)
    parser.add_argument("--sample-seed-base", type=int, default=2026094000)
    parser.add_argument("--run-log2-sd", type=float, default=0.05)
    parser.add_argument("--protein-log2-sd", type=float, default=0.15)
    parser.add_argument("--peptide-log2-sd", type=float, default=0.08)
    parser.add_argument("--treatment-up-fraction", type=float, default=0.15)
    parser.add_argument("--treatment-down-fraction", type=float, default=0.15)
    parser.add_argument("--treatment-min-abs-log2fc", type=float, default=0.5)
    parser.add_argument("--treatment-max-abs-log2fc", type=float, default=2.0)
    args = parser.parse_args()

    if not args.blueprint_db.is_file():
        raise SystemExit(f"Blueprint DB does not exist: {args.blueprint_db}")
    if args.control_runs < 1 or args.treatment_runs < 1:
        parser.error("--control-runs and --treatment-runs must both be positive")
    if args.study_seed < 0 or args.sample_seed_base < 0:
        parser.error("seeds must be non-negative")
    for value, label in (
        (args.run_log2_sd, "run-log2-sd"),
        (args.protein_log2_sd, "protein-log2-sd"),
        (args.peptide_log2_sd, "peptide-log2-sd"),
    ):
        if value < 0:
            parser.error(f"--{label} must be non-negative")
    if not 0 <= args.treatment_up_fraction <= 1 or not 0 <= args.treatment_down_fraction <= 1:
        parser.error("treatment fractions must be in [0, 1]")
    if args.treatment_up_fraction + args.treatment_down_fraction > 1:
        parser.error("treatment up/down fractions must sum to <= 1")
    if not 0 < args.treatment_min_abs_log2fc <= args.treatment_max_abs_log2fc:
        parser.error("treatment effect bounds must satisfy 0 < min <= max")

    with sqlite3.connect(args.blueprint_db) as connection:
        peptides, columns_map = load_peptides(connection)
    if not peptides:
        raise SystemExit("Blueprint peptides table is empty")

    design = assign_protein_design(
        [row.protein for row in peptides],
        up_fraction=args.treatment_up_fraction,
        down_fraction=args.treatment_down_fraction,
        min_abs_log2fc=args.treatment_min_abs_log2fc,
        max_abs_log2fc=args.treatment_max_abs_log2fc,
        seed=args.study_seed + 1,
    )

    design_rows = [
        {
            "ProteinId": item.protein,
            "TreatmentClass": item.treatment_class,
            "DesignLog2FC": f"{item.design_log2fc:.8f}",
            "DesignFoldChange": f"{2.0 ** item.design_log2fc:.8f}",
        }
        for item in sorted(design.values(), key=lambda value: value.protein)
    ]
    write_tsv(args.design_truth_out, design_rows)

    run_specs: list[tuple[str, int]] = [
        *(('control', index) for index in range(1, args.control_runs + 1)),
        *(('treatment', index) for index in range(1, args.treatment_runs + 1)),
    ]
    manifest_rows: list[dict[str, object]] = []
    abundance_rows: list[dict[str, object]] = []
    args.output_root.mkdir(parents=True, exist_ok=True)

    for ordinal, (condition, replicate) in enumerate(run_specs, start=1):
        prefix = "C" if condition == "control" else "T"
        run_id = f"{prefix}{replicate:02d}"
        run_name = f"OpenSwathTimSim_{condition}_{replicate:02d}"
        sample_seed = args.sample_seed_base + ordinal
        abundance_seed = args.study_seed + 1000 + ordinal
        if sample_seed > 4294967295 or abundance_seed > 4294967295:
            raise SystemExit("Derived run seed exceeds uint32 range")
        rng = random.Random(abundance_seed)
        run_global = rng.gauss(0.0, args.run_log2_sd)
        protein_noise = {
            protein: rng.gauss(0.0, args.protein_log2_sd)
            for protein in sorted(design)
        }

        input_dir = args.output_root / run_name
        input_dir.mkdir(parents=True, exist_ok=True)
        input_db = input_dir / "synthetic_data.db"
        shutil.copy2(args.blueprint_db, input_db)
        new_events: dict[int, int] = {}

        for peptide in peptides:
            peptide_noise = rng.gauss(0.0, args.peptide_log2_sd)
            condition_effect = design[peptide.protein].design_log2fc if condition == "treatment" else 0.0
            total_log2 = run_global + protein_noise[peptide.protein] + peptide_noise + condition_effect
            baseline_events = max(1, int(round(peptide.events)))
            realized_events = max(1, int(round(baseline_events * (2.0 ** total_log2))))
            new_events[peptide.peptide_id] = realized_events
            abundance_rows.append({
                "RunId": run_id,
                "RunName": run_name,
                "Condition": condition,
                "Replicate": replicate,
                "PeptideId": peptide.peptide_id,
                "PeptideSequence": peptide.sequence,
                "ProteinId": peptide.protein,
                "TreatmentClass": design[peptide.protein].treatment_class,
                "DesignLog2FC": f"{design[peptide.protein].design_log2fc:.8f}",
                "BaselineEvents": baseline_events,
                "RunGlobalLog2Effect": f"{run_global:.8f}",
                "ProteinRandomLog2Effect": f"{protein_noise[peptide.protein]:.8f}",
                "PeptideRandomLog2Effect": f"{peptide_noise:.8f}",
                "ConditionLog2Effect": f"{condition_effect:.8f}",
                "TotalLog2FactorVsBlueprint": f"{total_log2:.8f}",
                "RealizedInputEvents": realized_events,
            })

        update_database(input_db, columns_map, new_events)
        manifest_rows.append({
            "RunOrdinal": ordinal,
            "RunId": run_id,
            "RunName": run_name,
            "Condition": condition,
            "Replicate": replicate,
            "TimSimSampleSeed": sample_seed,
            "AbundanceSeed": abundance_seed,
            "InputDirectory": str(input_dir),
            "InputDatabase": str(input_db),
        })

    write_tsv(args.study_manifest_out, manifest_rows)
    write_tsv(args.abundance_truth_out, abundance_rows)

    class_counts = {key: 0 for key in ("up", "down", "unchanged")}
    for item in design.values():
        class_counts[item.treatment_class] += 1
    print(
        f"Study design: {args.control_runs} control + {args.treatment_runs} treatment runs; "
        f"{len(peptides)} peptides from {len(design)} proteins"
    )
    print(
        "Treatment proteins: "
        f"{class_counts['up']} up, {class_counts['down']} down, {class_counts['unchanged']} unchanged"
    )
    print(f"Study manifest:  {args.study_manifest_out}")
    print(f"Design truth:    {args.design_truth_out}")
    print(f"Abundance truth: {args.abundance_truth_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
