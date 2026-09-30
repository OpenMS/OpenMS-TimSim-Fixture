#!/usr/bin/env python3
"""Validate that the three TimSim runs form one comparable experiment."""
from __future__ import annotations

import argparse
import sqlite3
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


def load_run(path: Path) -> tuple[dict[str, tuple[float, float]], dict[tuple[str, int], tuple[float, float]]]:
    with sqlite3.connect(path) as connection:
        peptide_cols = columns(connection, "peptides")
        ion_cols = columns(connection, "ions")

        peptide_id_col = first_present(("peptide_id", "id"), peptide_cols, "peptide ID")
        sequence_col = first_present(("sequence", "peptide"), peptide_cols, "sequence")
        events_col = first_present(("events", "total_events", "abundance"), peptide_cols, "events")
        rt_col = first_present(
            ("retention_time_gru_predictor", "retention_time", "rt", "predicted_rt"),
            peptide_cols,
            "retention time",
        )
        peptides_by_id: dict[int, str] = {}
        peptides: dict[str, tuple[float, float]] = {}
        query = (
            f"SELECT {quote_identifier(peptide_id_col)}, {quote_identifier(sequence_col)}, "
            f"{quote_identifier(events_col)}, {quote_identifier(rt_col)} FROM peptides"
        )
        for peptide_id, sequence, events, rt in connection.execute(query):
            sequence = str(sequence)
            peptides_by_id[int(peptide_id)] = sequence
            peptides[sequence] = (float(events), float(rt))

        ion_peptide_col = first_present(("peptide_id",), ion_cols, "ion peptide ID")
        charge_col = first_present(("charge", "precursor_charge"), ion_cols, "charge")
        mz_col = first_present(("mz", "precursor_mz"), ion_cols, "precursor m/z")
        im_col = first_present(
            ("inv_mobility_gru_predictor", "inverse_mobility", "ion_mobility", "mobility"),
            ion_cols,
            "ion mobility",
        )
        ions: dict[tuple[str, int], tuple[float, float]] = {}
        query = (
            f"SELECT {quote_identifier(ion_peptide_col)}, {quote_identifier(charge_col)}, "
            f"{quote_identifier(mz_col)}, {quote_identifier(im_col)} FROM ions"
        )
        for peptide_id, charge, mz, mobility in connection.execute(query):
            sequence = peptides_by_id[int(peptide_id)]
            ions[(sequence, int(charge))] = (float(mz), float(mobility))

    return peptides, ions


def max_abs(values: list[float]) -> float:
    return max((abs(value) for value in values), default=0.0)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("db", type=Path, nargs=3, metavar=("RUN1_DB", "RUN2_DB", "RUN3_DB"))
    args = parser.parse_args()

    loaded = [load_run(path) for path in args.db]
    peptide_maps = [item[0] for item in loaded]
    ion_maps = [item[1] for item in loaded]

    baseline_sequences = set(peptide_maps[0])
    baseline_ions = set(ion_maps[0])
    for index in (1, 2):
        if set(peptide_maps[index]) != baseline_sequences:
            missing = baseline_sequences - set(peptide_maps[index])
            extra = set(peptide_maps[index]) - baseline_sequences
            raise SystemExit(
                f"Run {index + 1} peptide identities differ: {len(missing)} missing, {len(extra)} extra"
            )
        if set(ion_maps[index]) != baseline_ions:
            missing = baseline_ions - set(ion_maps[index])
            extra = set(ion_maps[index]) - baseline_ions
            raise SystemExit(
                f"Run {index + 1} precursor identities differ: {len(missing)} missing, {len(extra)} extra"
            )

    for index in (1, 2):
        rt_delta = [
            peptide_maps[index][sequence][1] - peptide_maps[0][sequence][1]
            for sequence in baseline_sequences
        ]
        mz_delta = [
            ion_maps[index][key][0] - ion_maps[0][key][0]
            for key in baseline_ions
        ]
        im_delta = [
            ion_maps[index][key][1] - ion_maps[0][key][1]
            for key in baseline_ions
        ]
        if max_abs(rt_delta) > 1e-8:
            raise SystemExit(f"Run {index + 1} changed RT values; max absolute delta={max_abs(rt_delta)}")
        if max_abs(mz_delta) > 1e-8:
            raise SystemExit(f"Run {index + 1} changed precursor m/z; max absolute delta={max_abs(mz_delta)}")
        if max_abs(im_delta) > 1e-8:
            raise SystemExit(f"Run {index + 1} changed ion mobility; max absolute delta={max_abs(im_delta)}")

    replicate_changed = sum(
        peptide_maps[1][sequence][0] != peptide_maps[0][sequence][0]
        for sequence in baseline_sequences
    )
    treatment_changed = sum(
        peptide_maps[2][sequence][0] != peptide_maps[0][sequence][0]
        for sequence in baseline_sequences
    )
    if replicate_changed == 0:
        raise SystemExit("Run 2 has no abundance changes")
    if treatment_changed == 0:
        raise SystemExit("Run 3 has no abundance changes")

    print(
        f"Comparable experiment OK: {len(baseline_sequences)} peptides, "
        f"{len(baseline_ions)} precursor ions"
    )
    print(f"RT/mz/IM identities are unchanged across all runs")
    print(f"Run 2 abundance changes: {replicate_changed}/{len(baseline_sequences)} peptides")
    print(f"Run 3 abundance changes: {treatment_changed}/{len(baseline_sequences)} peptides")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
