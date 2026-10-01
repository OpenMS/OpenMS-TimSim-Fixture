#!/usr/bin/env python3
"""Plan a large TimSim study without copying one source database per run.

The plan records deterministic run seeds and abundance-model parameters. Individual
runs are materialized from the immutable blueprint DB only when a worker executes.
This keeps persistent storage proportional to final raw data instead of run count.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
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
class ProteinDesign:
    protein: str
    treatment_class: str
    design_log2fc: float


def load_proteins(path: Path) -> tuple[list[str], int]:
    with sqlite3.connect(path) as connection:
        available = columns(connection, "peptides")
        protein_col = first_present(("protein", "protein_id", "protein_name"), available, "protein")
        proteins = [
            str(row[0])
            for row in connection.execute(
                f"SELECT DISTINCT {quote_identifier(protein_col)} FROM {quote_identifier('peptides')}"
            )
        ]
        peptide_count = int(connection.execute("SELECT COUNT(*) FROM peptides").fetchone()[0])
    return sorted(set(proteins)), peptide_count


def assign_protein_design(
    proteins: list[str],
    *,
    up_fraction: float,
    down_fraction: float,
    min_abs_log2fc: float,
    max_abs_log2fc: float,
    seed: int,
) -> dict[str, ProteinDesign]:
    rng = random.Random(seed)
    shuffled = proteins.copy()
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
    design: dict[str, ProteinDesign] = {}
    for protein in proteins:
        if protein in up:
            effect = rng.uniform(min_abs_log2fc, max_abs_log2fc)
            design[protein] = ProteinDesign(protein, "up", effect)
        elif protein in down:
            effect = -rng.uniform(min_abs_log2fc, max_abs_log2fc)
            design[protein] = ProteinDesign(protein, "down", effect)
        else:
            design[protein] = ProteinDesign(protein, "unchanged", 0.0)
    return design


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def parse_levels(text: str) -> list[float]:
    if not text.strip():
        return [0.0]
    values = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("input offset levels must contain at least one number")
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blueprint-db", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--selected-precursors", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--gradient-length", type=float, default=30.0)
    parser.add_argument("--control-runs", type=int, default=50)
    parser.add_argument("--treatment-runs", type=int, default=50)
    parser.add_argument("--profile", choices=("bulk", "single_cell"), default="bulk")
    parser.add_argument("--study-seed", type=int, default=2026100101)
    parser.add_argument("--sample-seed-base", type=int, default=2026101000)
    parser.add_argument("--run-log2-sd", type=float, default=0.05)
    parser.add_argument("--protein-log2-sd", type=float, default=0.15)
    parser.add_argument("--peptide-log2-sd", type=float, default=0.08)
    parser.add_argument("--cell-size-log2-sd", type=float, default=0.0)
    parser.add_argument("--input-log2-offset", type=float, default=0.0)
    parser.add_argument(
        "--input-log2-offset-levels",
        default="",
        help="Comma-separated offsets assigned evenly within each condition (single-cell calibration).",
    )
    parser.add_argument("--event-sampling", choices=("deterministic", "poisson"), default="deterministic")
    parser.add_argument("--allow-zero-events", action="store_true")
    parser.add_argument("--treatment-up-fraction", type=float, default=0.15)
    parser.add_argument("--treatment-down-fraction", type=float, default=0.15)
    parser.add_argument("--treatment-min-abs-log2fc", type=float, default=0.5)
    parser.add_argument("--treatment-max-abs-log2fc", type=float, default=2.0)
    args = parser.parse_args()

    if not args.blueprint_db.is_file():
        raise SystemExit(f"Blueprint DB does not exist: {args.blueprint_db}")
    if not args.selected_precursors.is_file():
        raise SystemExit(f"Selected precursor table does not exist: {args.selected_precursors}")
    if not args.fasta.is_file():
        raise SystemExit(f"Synthetic FASTA does not exist: {args.fasta}")
    if args.control_runs < 1 or args.treatment_runs < 1:
        parser.error("control/treatment run counts must both be positive")
    if args.study_seed < 0 or args.sample_seed_base < 0:
        parser.error("seeds must be non-negative")
    for value, label in (
        (args.run_log2_sd, "run-log2-sd"),
        (args.protein_log2_sd, "protein-log2-sd"),
        (args.peptide_log2_sd, "peptide-log2-sd"),
        (args.cell_size_log2_sd, "cell-size-log2-sd"),
    ):
        if value < 0:
            parser.error(f"--{label} must be non-negative")
    if args.profile == "single_cell" and args.event_sampling != "poisson":
        parser.error("single_cell profile requires --event-sampling poisson")
    if args.profile == "single_cell" and not args.allow_zero_events:
        parser.error("single_cell profile requires --allow-zero-events")
    if args.treatment_up_fraction + args.treatment_down_fraction > 1:
        parser.error("treatment fractions must sum to <= 1")
    if not 0 < args.treatment_min_abs_log2fc <= args.treatment_max_abs_log2fc:
        parser.error("treatment effect bounds must satisfy 0 < min <= max")

    levels = parse_levels(args.input_log2_offset_levels)
    proteins, peptide_count = load_proteins(args.blueprint_db)
    design = assign_protein_design(
        proteins,
        up_fraction=args.treatment_up_fraction,
        down_fraction=args.treatment_down_fraction,
        min_abs_log2fc=args.treatment_min_abs_log2fc,
        max_abs_log2fc=args.treatment_max_abs_log2fc,
        seed=args.study_seed + 1,
    )

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    design_path = output_root / "OpenSwathTimSim.study_design_truth.tsv"
    manifest_path = output_root / "OpenSwathTimSim.study_manifest.tsv"
    plan_path = output_root / "OpenSwathTimSim.materialized_study_plan.json"

    design_rows = [
        {
            "ProteinId": item.protein,
            "TreatmentClass": item.treatment_class,
            "DesignLog2FC": f"{item.design_log2fc:.8f}",
            "DesignFoldChange": f"{2.0 ** item.design_log2fc:.8f}",
        }
        for item in sorted(design.values(), key=lambda value: value.protein)
    ]
    write_tsv(design_path, design_rows)

    run_specs = [
        *(('control', index) for index in range(1, args.control_runs + 1)),
        *(('treatment', index) for index in range(1, args.treatment_runs + 1)),
    ]
    manifest_rows: list[dict[str, object]] = []
    for ordinal, (condition, replicate) in enumerate(run_specs, start=1):
        prefix = "C" if condition == "control" else "T"
        run_id = f"{prefix}{replicate:03d}"
        run_name = f"OpenSwathTimSim_{condition}_{replicate:03d}"
        sample_seed = args.sample_seed_base + ordinal
        abundance_seed = args.study_seed + 1000 + ordinal
        if sample_seed > 4294967295 or abundance_seed > 4294967295:
            raise SystemExit("Derived run seed exceeds uint32 range")
        level_index = (replicate - 1) % len(levels)
        level_offset = levels[level_index] + args.input_log2_offset
        manifest_rows.append({
            "RunOrdinal": ordinal,
            "RunId": run_id,
            "RunName": run_name,
            "Condition": condition,
            "Replicate": replicate,
            "TimSimSampleSeed": sample_seed,
            "AbundanceSeed": abundance_seed,
            "AbundanceProfile": args.profile,
            "InputLevel": f"L{level_index + 1}",
            "GlobalInputLog2Offset": f"{level_offset:.8f}",
            "RunLog2SD": f"{args.run_log2_sd:.8f}",
            "ProteinLog2SD": f"{args.protein_log2_sd:.8f}",
            "PeptideLog2SD": f"{args.peptide_log2_sd:.8f}",
            "CellSizeLog2SD": f"{args.cell_size_log2_sd:.8f}",
            "EventSampling": args.event_sampling,
            "AllowZeroEvents": int(args.allow_zero_events),
        })
    write_tsv(manifest_path, manifest_rows)

    selected_count = sum(1 for _ in csv.DictReader(args.selected_precursors.open(newline="", encoding="utf-8"), delimiter="\t"))
    class_counts = {key: 0 for key in ("up", "down", "unchanged")}
    for item in design.values():
        class_counts[item.treatment_class] += 1

    plan = {
        "plan_version": "materialized_study_v1",
        "materialization": "per_run_scratch",
        "blueprint_db": str(args.blueprint_db.resolve()),
        "selected_precursors": str(args.selected_precursors.resolve()),
        "fasta": str(args.fasta.resolve()),
        "gradient_length_seconds": args.gradient_length,
        "profile": args.profile,
        "control_runs": args.control_runs,
        "treatment_runs": args.treatment_runs,
        "total_runs": len(manifest_rows),
        "blueprint_peptides": peptide_count,
        "blueprint_proteins": len(proteins),
        "selected_precursors_count": selected_count,
        "input_log2_offset_levels": levels,
        "base_input_log2_offset": args.input_log2_offset,
        "event_sampling": args.event_sampling,
        "allow_zero_events": args.allow_zero_events,
        "run_log2_sd": args.run_log2_sd,
        "protein_log2_sd": args.protein_log2_sd,
        "peptide_log2_sd": args.peptide_log2_sd,
        "cell_size_log2_sd": args.cell_size_log2_sd,
        "treatment": {
            "up_fraction": args.treatment_up_fraction,
            "down_fraction": args.treatment_down_fraction,
            "min_abs_log2fc": args.treatment_min_abs_log2fc,
            "max_abs_log2fc": args.treatment_max_abs_log2fc,
            "protein_class_counts": class_counts,
        },
        "artifacts": {
            "manifest": manifest_path.name,
            "design_truth": design_path.name,
            "run_input_truth_dir": "run_input_truth",
            "run_realized_truth_dir": "run_realized_truth",
            "run_stats_dir": "run_stats",
        },
    }
    plan_path.write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")

    print(
        f"Materialized study plan: {args.control_runs} control + {args.treatment_runs} treatment "
        f"runs; profile={args.profile}; selected targets={selected_count}"
    )
    print(f"Blueprint: {peptide_count} peptides / {len(proteins)} proteins")
    print(f"Manifest: {manifest_path}")
    print(f"Design truth: {design_path}")
    print(f"Plan: {plan_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
