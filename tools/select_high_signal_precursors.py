#!/usr/bin/env python3
"""Select a simulator-defined high-signal TimSim precursor correctness set.

Selection is intentionally independent of OpenDIA identification, scores, selected
features, Percolator, or q-values. The only inputs are TimSim synthetic_data.db
files plus the simulator/library geometry encoded in the baseline database.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Any, Iterable

from export_openswath_tsv import (
    Candidate,
    columns,
    first_present,
    optional_present,
    parse_array,
    read_candidates,
)

RUNS = (
    ("baseline", "OpenSwathTimSim_01_baseline"),
    ("control_replicate", "OpenSwathTimSim_02_biological_replicate"),
    ("treatment", "OpenSwathTimSim_03_treatment"),
)


@dataclass(frozen=True)
class RunTruth:
    run_role: str
    run_name: str
    peptide_id: int
    ion_id: int
    sequence: str
    protein: str
    charge: int
    precursor_mz: float
    assay_rt: float
    assay_im: float
    peptide_events: float
    ion_relative_abundance: float
    frame_abundance_sum: float
    scan_abundance_sum: float
    realized_event_proxy: float
    realized_rt_apex: float
    realized_rt_centroid: float
    realized_im_apex: float
    realized_im_centroid: float
    frame_points: int
    scan_points: int
    acquisition_rt_start: float
    acquisition_rt_end: float


@dataclass(frozen=True)
class WindowRecord:
    window_index: int
    window_group: int
    mz_lower: float
    mz_upper: float
    im_lower: float
    im_upper: float

    @property
    def label(self) -> str:
        return f"map_{self.window_index:02d}"


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def weighted_coordinates(
    occurrences_raw: Any,
    abundance_raw: Any,
    coordinate_by_index: dict[int, float],
) -> tuple[float, float, float, int]:
    occurrences = [int(round(value)) for value in parse_array(occurrences_raw)]
    abundances = [float(value) for value in parse_array(abundance_raw)]
    if len(occurrences) != len(abundances):
        raise RuntimeError(
            f"Occurrence/abundance length mismatch: {len(occurrences)} != {len(abundances)}"
        )

    points = [
        (coordinate_by_index[index], weight)
        for index, weight in zip(occurrences, abundances)
        if index in coordinate_by_index and math.isfinite(weight) and weight >= 0.0
    ]
    abundance_sum = sum(weight for _, weight in points)
    if not points or abundance_sum <= 0.0:
        return abundance_sum, math.nan, math.nan, len(points)

    apex_coordinate = max(points, key=lambda item: item[1])[0]
    centroid = sum(coordinate * weight for coordinate, weight in points) / abundance_sum
    return abundance_sum, float(apex_coordinate), float(centroid), len(points)


def load_run_truth(path: Path, run_role: str, run_name: str) -> dict[tuple[str, int], RunTruth]:
    if not path.is_file():
        raise FileNotFoundError(f"TimSim database does not exist: {path}")

    with sqlite3.connect(path) as connection:
        peptide_cols = columns(connection, "peptides")
        ion_cols = columns(connection, "ions")
        frame_cols = columns(connection, "frames")
        scan_cols = columns(connection, "scans")

        peptide_id_col = first_present(("peptide_id", "id"), peptide_cols, label="peptide ID")
        sequence_col = first_present(("sequence", "peptide"), peptide_cols, label="sequence")
        protein_col = optional_present(("protein", "protein_id", "protein_name"), peptide_cols)
        events_col = first_present(("events", "total_events", "abundance"), peptide_cols, label="peptide abundance")
        rt_col = first_present(
            ("retention_time_gru_predictor", "retention_time", "rt", "predicted_rt"),
            peptide_cols,
            label="assay RT",
        )
        frame_occurrence_col = first_present(("frame_occurrence",), peptide_cols, label="frame occurrence")
        frame_abundance_col = first_present(("frame_abundance",), peptide_cols, label="frame abundance")
        decoy_col = optional_present(("decoy", "is_decoy"), peptide_cols)

        frame_id_col = first_present(("frame_id", "id"), frame_cols, label="frame ID")
        frame_time_col = first_present(("time", "rt"), frame_cols, label="frame time")
        frame_coordinates = {
            int(frame_id): float(time)
            for frame_id, time in connection.execute(
                f"SELECT {quote_identifier(frame_id_col)}, {quote_identifier(frame_time_col)} FROM frames"
            )
        }
        if not frame_coordinates:
            raise RuntimeError(f"No frame coordinates found in {path}")
        acquisition_start = min(frame_coordinates.values())
        acquisition_end = max(frame_coordinates.values())

        scan_id_col = first_present(("scan", "scan_id", "id"), scan_cols, label="scan ID")
        scan_im_col = first_present(
            ("mobility", "inv_mobility", "one_over_k0", "inverse_mobility"),
            scan_cols,
            label="scan mobility",
        )
        scan_coordinates = {
            int(scan): float(mobility)
            for scan, mobility in connection.execute(
                f"SELECT {quote_identifier(scan_id_col)}, {quote_identifier(scan_im_col)} FROM scans"
            )
        }

        peptide_select = [
            peptide_id_col,
            sequence_col,
            events_col,
            rt_col,
            frame_occurrence_col,
            frame_abundance_col,
        ]
        if protein_col:
            peptide_select.append(protein_col)
        if decoy_col:
            peptide_select.append(decoy_col)

        peptides: dict[int, dict[str, Any]] = {}
        query = "SELECT " + ", ".join(quote_identifier(value) for value in peptide_select) + " FROM peptides"
        for row in connection.execute(query):
            data = dict(zip(peptide_select, row))
            if decoy_col and bool(data.get(decoy_col)):
                continue
            peptide_id = int(data[peptide_id_col])
            sequence = str(data[sequence_col]).strip().upper()
            frame_sum, rt_apex, rt_centroid, frame_points = weighted_coordinates(
                data[frame_occurrence_col], data[frame_abundance_col], frame_coordinates
            )
            peptides[peptide_id] = {
                "sequence": sequence,
                "protein": str(data.get(protein_col) or f"TIMSIM_PROTEIN_{peptide_id}"),
                "events": float(data[events_col]),
                "assay_rt": float(data[rt_col]),
                "frame_sum": frame_sum,
                "rt_apex": rt_apex,
                "rt_centroid": rt_centroid,
                "frame_points": frame_points,
            }

        ion_id_col = first_present(("ion_id", "id"), ion_cols, label="ion ID")
        ion_peptide_col = first_present(("peptide_id",), ion_cols, label="ion peptide ID")
        charge_col = first_present(("charge", "precursor_charge"), ion_cols, label="charge")
        mz_col = first_present(("mz", "precursor_mz"), ion_cols, label="precursor m/z")
        abundance_col = first_present(("relative_abundance",), ion_cols, label="ion relative abundance")
        im_col = first_present(
            ("inv_mobility_gru_predictor", "inverse_mobility", "ion_mobility", "mobility"),
            ion_cols,
            label="assay ion mobility",
        )
        scan_occurrence_col = first_present(("scan_occurrence",), ion_cols, label="scan occurrence")
        scan_abundance_col = first_present(("scan_abundance",), ion_cols, label="scan abundance")

        query = (
            "SELECT "
            + ", ".join(
                quote_identifier(value)
                for value in (
                    ion_id_col,
                    ion_peptide_col,
                    charge_col,
                    mz_col,
                    abundance_col,
                    im_col,
                    scan_occurrence_col,
                    scan_abundance_col,
                )
            )
            + " FROM ions"
        )

        result: dict[tuple[str, int], RunTruth] = {}
        for ion_id, peptide_id, charge, mz, ion_fraction, assay_im, scan_occurrence, scan_abundance in connection.execute(query):
            peptide_id = int(peptide_id)
            if peptide_id not in peptides:
                continue
            peptide = peptides[peptide_id]
            scan_sum, im_apex, im_centroid, scan_points = weighted_coordinates(
                scan_occurrence, scan_abundance, scan_coordinates
            )
            realized_proxy = (
                peptide["events"]
                * peptide["frame_sum"]
                * float(ion_fraction)
                * scan_sum
            )
            key = (peptide["sequence"], int(charge))
            if key in result:
                raise RuntimeError(f"Duplicate precursor key {key} in {path}")
            result[key] = RunTruth(
                run_role=run_role,
                run_name=run_name,
                peptide_id=peptide_id,
                ion_id=int(ion_id),
                sequence=peptide["sequence"],
                protein=peptide["protein"],
                charge=int(charge),
                precursor_mz=float(mz),
                assay_rt=peptide["assay_rt"],
                assay_im=float(assay_im),
                peptide_events=peptide["events"],
                ion_relative_abundance=float(ion_fraction),
                frame_abundance_sum=peptide["frame_sum"],
                scan_abundance_sum=scan_sum,
                realized_event_proxy=realized_proxy,
                realized_rt_apex=peptide["rt_apex"],
                realized_rt_centroid=peptide["rt_centroid"],
                realized_im_apex=im_apex,
                realized_im_centroid=im_centroid,
                frame_points=peptide["frame_points"],
                scan_points=scan_points,
                acquisition_rt_start=acquisition_start,
                acquisition_rt_end=acquisition_end,
            )
        return result


def load_windows_with_identity(connection: sqlite3.Connection) -> list[WindowRecord]:
    window_cols = columns(connection, "dia_ms_ms_windows")
    scan_cols = columns(connection, "scans")
    group_col = optional_present(("window_group", "WindowGroup"), window_cols)
    start_col = first_present(("scan_start", "ScanNumBegin"), window_cols, label="window scan start")
    end_col = first_present(("scan_end", "ScanNumEnd"), window_cols, label="window scan end")
    center_col = first_present(("isolation_mz", "IsolationMz"), window_cols, label="window isolation m/z")
    width_col = first_present(("isolation_width", "IsolationWidth"), window_cols, label="window isolation width")
    scan_id_col = first_present(("scan", "scan_id", "Scan"), scan_cols, label="scan index")
    mobility_col = first_present(
        ("mobility", "inv_mobility", "one_over_k0", "inverse_mobility"),
        scan_cols,
        label="scan mobility",
    )

    scan_to_im = {
        int(scan): float(mobility)
        for scan, mobility in connection.execute(
            f"SELECT {quote_identifier(scan_id_col)}, {quote_identifier(mobility_col)} FROM scans"
        )
    }
    available_scans = sorted(scan_to_im)
    select = [start_col, end_col, center_col, width_col]
    if group_col:
        select.append(group_col)
    query = "SELECT " + ", ".join(quote_identifier(value) for value in select) + " FROM dia_ms_ms_windows"

    result: list[WindowRecord] = []
    for window_index, row in enumerate(connection.execute(query), start=1):
        data = dict(zip(select, row))
        start = int(data[start_col])
        end = int(data[end_col])
        start_key = min(available_scans, key=lambda value: abs(value - start))
        end_key = min(available_scans, key=lambda value: abs(value - max(start, end - 1)))
        im_a, im_b = scan_to_im[start_key], scan_to_im[end_key]
        center = float(data[center_col])
        width = float(data[width_col])
        result.append(
            WindowRecord(
                window_index=window_index,
                window_group=int(data[group_col]) if group_col else window_index,
                mz_lower=center - width / 2.0,
                mz_upper=center + width / 2.0,
                im_lower=min(im_a, im_b),
                im_upper=max(im_a, im_b),
            )
        )
    return result


def assign_window(
    candidate: Candidate,
    windows: list[WindowRecord],
    mz_margin: float,
    im_margin: float,
) -> WindowRecord | None:
    matches: list[tuple[float, WindowRecord]] = []
    for window in windows:
        if not (
            window.mz_lower + mz_margin < candidate.precursor_mz < window.mz_upper - mz_margin
            and window.im_lower + im_margin < candidate.precursor_im < window.im_upper - im_margin
        ):
            continue
        mz_room = min(
            candidate.precursor_mz - (window.mz_lower + mz_margin),
            (window.mz_upper - mz_margin) - candidate.precursor_mz,
        )
        im_room = min(
            candidate.precursor_im - (window.im_lower + im_margin),
            (window.im_upper - im_margin) - candidate.precursor_im,
        )
        mz_scale = max((window.mz_upper - window.mz_lower) - 2.0 * mz_margin, 1e-12)
        im_scale = max((window.im_upper - window.im_lower) - 2.0 * im_margin, 1e-12)
        score = min(mz_room / mz_scale, im_room / im_scale)
        matches.append((score, window))
    if not matches:
        return None
    return max(matches, key=lambda item: (item[0], -item[1].window_index))[1]


def assign_rank_bins(rows: list[dict[str, Any]], value_key: str, output_key: str, bins: int) -> None:
    if bins < 1:
        raise ValueError("Number of bins must be positive")
    ordered = sorted(range(len(rows)), key=lambda index: (float(rows[index][value_key]), rows[index]["sequence"]))
    total = len(ordered)
    for rank, index in enumerate(ordered):
        rows[index][output_key] = min(bins - 1, rank * bins // max(total, 1))


def proportional_targets(values: list[Any], total_selected: int) -> dict[Any, int]:
    counts = Counter(values)
    total_available = sum(counts.values())
    if total_available == 0:
        return {}
    raw = {key: total_selected * count / total_available for key, count in counts.items()}
    targets = {key: min(counts[key], int(math.floor(value))) for key, value in raw.items()}
    remaining = total_selected - sum(targets.values())
    order = sorted(
        counts,
        key=lambda key: (raw[key] - math.floor(raw[key]), counts[key], str(key)),
        reverse=True,
    )
    for key in order:
        if remaining <= 0:
            break
        if targets[key] < counts[key]:
            targets[key] += 1
            remaining -= 1
    return targets


def stratified_select(rows: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    if len(rows) < count:
        raise ValueError(f"Only {len(rows)} eligible unique precursors are available for {count} selections")

    dimensions = {
        "rt_bin": 3.0,
        "mz_bin": 2.0,
        "im_bin": 2.0,
        "swath_window_index": 2.5,
        "charge": 1.5,
    }
    targets = {
        dimension: proportional_targets([row[dimension] for row in rows], count)
        for dimension in dimensions
    }
    selected_counts = {dimension: Counter() for dimension in dimensions}
    selected_proteins: Counter[str] = Counter()

    log_signal = [math.log10(max(float(row["min_realized_event_proxy"]), 1.0)) for row in rows]
    signal_min = min(log_signal)
    signal_max = max(log_signal)
    signal_span = max(signal_max - signal_min, 1e-12)

    remaining = list(rows)
    selected: list[dict[str, Any]] = []
    for selection_rank in range(1, count + 1):
        best_index = -1
        best_score: tuple[float, float, str] | None = None
        for index, row in enumerate(remaining):
            score = 0.0
            for dimension, weight in dimensions.items():
                category = row[dimension]
                target = max(targets[dimension].get(category, 0), 1)
                deficit = (target - selected_counts[dimension][category]) / target
                score += weight * deficit

            # Encourage broad protein coverage without making it a hard constraint.
            if selected_proteins[row["protein_id"]] == 0:
                score += 0.75
            else:
                score += 0.15 / (1.0 + selected_proteins[row["protein_id"]])

            normalized_signal = (
                math.log10(max(float(row["min_realized_event_proxy"]), 1.0)) - signal_min
            ) / signal_span
            score += 0.10 * normalized_signal

            # Deterministic tie-breakers prefer higher simulator signal, then key order.
            tie = (score, float(row["min_realized_event_proxy"]), row["precursor_key"])
            if best_score is None or tie > best_score:
                best_score = tie
                best_index = index

        chosen = remaining.pop(best_index)
        chosen["selection_rank"] = selection_rank
        selected.append(chosen)
        for dimension in dimensions:
            selected_counts[dimension][chosen[dimension]] += 1
        selected_proteins[chosen["protein_id"]] += 1

    return selected


def describe(values: Iterable[float]) -> dict[str, float]:
    data = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not data:
        return {key: math.nan for key in ("min", "q25", "median", "q75", "max")}

    def quantile(fraction: float) -> float:
        if len(data) == 1:
            return data[0]
        position = fraction * (len(data) - 1)
        lower = int(math.floor(position))
        upper = int(math.ceil(position))
        if lower == upper:
            return data[lower]
        weight = position - lower
        return data[lower] * (1.0 - weight) + data[upper] * weight

    return {
        "min": data[0],
        "q25": quantile(0.25),
        "median": quantile(0.50),
        "q75": quantile(0.75),
        "max": data[-1],
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        if fieldnames:
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t").writeheader()
        return
    fields = fieldnames or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-db", type=Path, required=True)
    parser.add_argument("--control-db", type=Path, required=True)
    parser.add_argument("--treatment-db", type=Path, required=True)
    parser.add_argument("--selected-out", type=Path, required=True)
    parser.add_argument("--truth-out", type=Path, required=True)
    parser.add_argument("--candidate-qc-out", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    parser.add_argument("--qc-report", type=Path, required=True)
    parser.add_argument("--precursors", type=int, default=500)
    parser.add_argument("--min-realized-event-proxy", type=float, default=50000.0)
    parser.add_argument("--min-frame-abundance-sum", type=float, default=0.90)
    parser.add_argument("--min-scan-abundance-sum", type=float, default=0.95)
    parser.add_argument("--min-ion-relative-abundance", type=float, default=0.25)
    parser.add_argument("--rt-edge-margin-fraction", type=float, default=0.05)
    parser.add_argument(
        "--rt-edge-min-seconds",
        type=float,
        default=6.0,
        help=(
            "Minimum absolute realized-RT-apex margin from each acquisition edge. "
            "The effective margin is max(acquisition_duration * --rt-edge-margin-fraction, "
            "--rt-edge-min-seconds). Default: 6.0 s."
        ),
    )
    parser.add_argument("--minimum-fragments", type=int, default=8)
    parser.add_argument("--min-product-mz", type=float, default=350.0)
    parser.add_argument("--max-product-mz", type=float, default=2000.0)
    parser.add_argument("--mz-window-margin", type=float, default=1.0)
    parser.add_argument("--im-window-margin", type=float, default=0.005)
    parser.add_argument("--rt-bins", type=int, default=5)
    parser.add_argument("--mz-bins", type=int, default=5)
    parser.add_argument("--im-bins", type=int, default=5)
    args = parser.parse_args()

    if args.precursors < 1:
        parser.error("--precursors must be positive")
    if not 0.0 <= args.rt_edge_margin_fraction < 0.5:
        parser.error("--rt-edge-margin-fraction must be in [0, 0.5)")
    if args.rt_edge_min_seconds < 0.0:
        parser.error("--rt-edge-min-seconds must be non-negative")

    db_paths = {
        "baseline": args.baseline_db,
        "control_replicate": args.control_db,
        "treatment": args.treatment_db,
    }
    run_truth = {
        role: load_run_truth(db_paths[role], role, run_name)
        for role, run_name in RUNS
    }

    with sqlite3.connect(args.baseline_db) as connection:
        windows = load_windows_with_identity(connection)
        fragment_candidates = read_candidates(
            connection,
            min_product_mz=args.min_product_mz,
            max_product_mz=args.max_product_mz,
            minimum_fragments=args.minimum_fragments,
        )
        total_ions = int(connection.execute("SELECT COUNT(*) FROM ions").fetchone()[0])

    candidate_rows: list[dict[str, Any]] = []
    stage_counts = Counter()
    for candidate in fragment_candidates:
        stage_counts["sufficient_fragments"] += 1
        key = (candidate.sequence, candidate.precursor_charge)
        precursor_key = f"{candidate.sequence}/{candidate.precursor_charge}"
        window = assign_window(candidate, windows, args.mz_window_margin, args.im_window_margin)
        geometry_ok = window is not None
        if geometry_ok:
            stage_counts["safe_dia_geometry"] += 1

        per_run = [run_truth[role].get(key) for role, _ in RUNS]
        all_runs_present = all(value is not None for value in per_run)
        if geometry_ok and all_runs_present:
            stage_counts["present_all_runs"] += 1

        valid_runs = [value for value in per_run if value is not None]
        min_proxy = min((value.realized_event_proxy for value in valid_runs), default=math.nan)
        min_frame = min((value.frame_abundance_sum for value in valid_runs), default=math.nan)
        min_scan = min((value.scan_abundance_sum for value in valid_runs), default=math.nan)
        min_ion = min((value.ion_relative_abundance for value in valid_runs), default=math.nan)

        signal_ok = all_runs_present and min_proxy >= args.min_realized_event_proxy
        frame_ok = all_runs_present and min_frame >= args.min_frame_abundance_sum
        scan_ok = all_runs_present and min_scan >= args.min_scan_abundance_sum
        ion_ok = all_runs_present and min_ion >= args.min_ion_relative_abundance
        rt_edge_ok = all_runs_present
        rt_edge_margin = math.nan
        rt_edge_margins: list[float] = []
        if all_runs_present:
            for value in valid_runs:
                duration = value.acquisition_rt_end - value.acquisition_rt_start
                fractional_margin = duration * args.rt_edge_margin_fraction
                margin = max(fractional_margin, args.rt_edge_min_seconds)
                rt_edge_margins.append(margin)
                if not (
                    math.isfinite(value.realized_rt_apex)
                    and value.acquisition_rt_start + margin <= value.realized_rt_apex <= value.acquisition_rt_end - margin
                ):
                    rt_edge_ok = False
            if rt_edge_margins:
                # Candidate-level provenance records the most stringent margin
                # used across the three runs. Truth rows below record run-specific
                # effective margins.
                rt_edge_margin = max(rt_edge_margins)

        cumulative = geometry_ok and all_runs_present
        for name, passed in (
            ("realized_signal", signal_ok),
            ("frame_mass", frame_ok),
            ("scan_mass", scan_ok),
            ("ion_fraction", ion_ok),
            ("rt_edge", rt_edge_ok),
        ):
            cumulative = cumulative and passed
            if cumulative:
                stage_counts[name] += 1

        baseline = run_truth["baseline"].get(key)
        row = {
            "precursor_key": precursor_key,
            "sequence": candidate.sequence,
            "protein_id": candidate.protein,
            "charge": candidate.precursor_charge,
            "precursor_mz": candidate.precursor_mz,
            "assay_rt": candidate.rt_seconds,
            "assay_im": candidate.precursor_im,
            "usable_fragments": len(candidate.fragments),
            "total_fragment_intensity": candidate.total_fragment_intensity,
            "swath_window_index": window.window_index if window else "",
            "swath_window_group": window.window_group if window else "",
            "swath_window_label": window.label if window else "",
            "swath_mz_lower": window.mz_lower if window else math.nan,
            "swath_mz_upper": window.mz_upper if window else math.nan,
            "swath_im_lower": window.im_lower if window else math.nan,
            "swath_im_upper": window.im_upper if window else math.nan,
            "min_realized_event_proxy": min_proxy,
            "min_frame_abundance_sum": min_frame,
            "min_scan_abundance_sum": min_scan,
            "min_ion_relative_abundance": min_ion,
            "realized_rt_apex_baseline": baseline.realized_rt_apex if baseline else math.nan,
            "realized_im_apex_baseline": baseline.realized_im_apex if baseline else math.nan,
            "rt_edge_margin_seconds": rt_edge_margin,
            "safe_dia_geometry": int(geometry_ok),
            "present_all_runs": int(all_runs_present),
            "signal_ok": int(signal_ok),
            "frame_ok": int(frame_ok),
            "scan_ok": int(scan_ok),
            "ion_fraction_ok": int(ion_ok),
            "rt_edge_ok": int(rt_edge_ok),
            "eligible": int(geometry_ok and all_runs_present and signal_ok and frame_ok and scan_ok and ion_ok and rt_edge_ok),
        }
        candidate_rows.append(row)

    args.candidate_qc_out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(args.candidate_qc_out, candidate_rows)

    eligible_ions = [row for row in candidate_rows if row["eligible"]]
    # Preserve one precursor charge per peptide, prioritizing the charge with the
    # largest minimum realized signal across all three runs.
    best_by_sequence: dict[str, dict[str, Any]] = {}
    for row in eligible_ions:
        previous = best_by_sequence.get(row["sequence"])
        ranking = (
            float(row["min_realized_event_proxy"]),
            float(row["min_ion_relative_abundance"]),
            float(row["total_fragment_intensity"]),
            -int(row["charge"]),
        )
        previous_ranking = None if previous is None else (
            float(previous["min_realized_event_proxy"]),
            float(previous["min_ion_relative_abundance"]),
            float(previous["total_fragment_intensity"]),
            -int(previous["charge"]),
        )
        if previous is None or ranking > previous_ranking:
            best_by_sequence[row["sequence"]] = row.copy()

    eligible = list(best_by_sequence.values())
    if eligible:
        assign_rank_bins(eligible, "assay_rt", "rt_bin", args.rt_bins)
        assign_rank_bins(eligible, "precursor_mz", "mz_bin", args.mz_bins)
        assign_rank_bins(eligible, "assay_im", "im_bin", args.im_bins)

    summary: dict[str, Any] = {
        "selection_is_opendia_independent": True,
        "requested_precursors": args.precursors,
        "total_ions_in_baseline_db": total_ions,
        "fragment_eligible_ions": len(fragment_candidates),
        "safe_geometry_ions": stage_counts["safe_dia_geometry"],
        "eligible_ions_after_all_filters": len(eligible_ions),
        "eligible_unique_peptides_after_best_charge": len(eligible),
        "criteria": {
            "min_realized_event_proxy_all_runs": args.min_realized_event_proxy,
            "min_frame_abundance_sum_all_runs": args.min_frame_abundance_sum,
            "min_scan_abundance_sum_all_runs": args.min_scan_abundance_sum,
            "min_ion_relative_abundance_all_runs": args.min_ion_relative_abundance,
            "rt_edge_margin_fraction_of_acquisition": args.rt_edge_margin_fraction,
            "rt_edge_min_seconds": args.rt_edge_min_seconds,
            "rt_edge_policy": "max(fraction_of_acquisition, absolute_seconds)",
            "minimum_usable_fragments": args.minimum_fragments,
            "mz_window_margin_th": args.mz_window_margin,
            "im_window_margin_1_over_k0": args.im_window_margin,
        },
        "cumulative_filter_counts": dict(stage_counts),
        "selected_precursors": 0,
        "status": "insufficient_candidates" if len(eligible) < args.precursors else "ready",
    }

    if len(eligible) < args.precursors:
        args.qc_json.parent.mkdir(parents=True, exist_ok=True)
        args.qc_json.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
        args.qc_report.parent.mkdir(parents=True, exist_ok=True)
        args.qc_report.write_text(
            "# High-signal correctness fixture selection QC\n\n"
            f"Selection failed without weakening criteria: **{len(eligible)}** unique peptide precursors qualified; "
            f"**{args.precursors}** are required. Increase the TimSim candidate population and regenerate all three runs.\n\n"
            f"Candidate-level QC: `{args.candidate_qc_out.name}`\n",
            encoding="utf-8",
        )
        raise SystemExit(
            f"Only {len(eligible)} unique high-signal precursor groups satisfy all simulator-only criteria; "
            f"{args.precursors} are required. Increase --simulated-peptides/--fasta-peptides and regenerate; "
            "do not relax the correctness thresholds merely to reach the target count."
        )

    selected = stratified_select(eligible, args.precursors)
    selected.sort(key=lambda row: int(row["selection_rank"]))
    selected_fields = [
        "selection_rank", "precursor_key", "sequence", "protein_id", "charge",
        "precursor_mz", "assay_rt", "assay_im", "usable_fragments",
        "total_fragment_intensity", "min_realized_event_proxy",
        "min_frame_abundance_sum", "min_scan_abundance_sum",
        "min_ion_relative_abundance", "realized_rt_apex_baseline",
        "realized_im_apex_baseline", "rt_edge_margin_seconds",
        "swath_window_index", "swath_window_group", "swath_window_label",
        "swath_mz_lower", "swath_mz_upper", "swath_im_lower", "swath_im_upper",
        "rt_bin", "mz_bin", "im_bin",
    ]
    write_tsv(args.selected_out, selected, selected_fields)

    truth_rows: list[dict[str, Any]] = []
    selected_by_key = {(row["sequence"], int(row["charge"])): row for row in selected}
    for role, run_name in RUNS:
        for key, selected_row in selected_by_key.items():
            run = run_truth[role][key]
            truth_rows.append({
                "selection_rank": selected_row["selection_rank"],
                "precursor_key": selected_row["precursor_key"],
                "TransitionGroupId": f"TIMSIM_{run.sequence}_{run.charge}",
                "PeptideSequence": run.sequence,
                "PrecursorCharge": run.charge,
                "ProteinId": run.protein,
                "PrecursorMz": run.precursor_mz,
                "AssayRT": run.assay_rt,
                "AssayIM": run.assay_im,
                "RunRole": run.run_role,
                "RunName": run.run_name,
                "PeptideId": run.peptide_id,
                "IonId": run.ion_id,
                "PeptideEvents": run.peptide_events,
                "IonRelativeAbundance": run.ion_relative_abundance,
                "FrameAbundanceSum": run.frame_abundance_sum,
                "ScanAbundanceSum": run.scan_abundance_sum,
                "RealizedEventProxy": run.realized_event_proxy,
                "RealizedRTApex": run.realized_rt_apex,
                "RealizedRTCentroid": run.realized_rt_centroid,
                "RealizedIMApexSQLite": run.realized_im_apex,
                "RealizedIMCentroidSQLite": run.realized_im_centroid,
                "FramePoints": run.frame_points,
                "ScanPoints": run.scan_points,
                "AcquisitionRTStart": run.acquisition_rt_start,
                "AcquisitionRTEnd": run.acquisition_rt_end,
                "RTEdgeMarginSeconds": max(
                    (run.acquisition_rt_end - run.acquisition_rt_start) * args.rt_edge_margin_fraction,
                    args.rt_edge_min_seconds,
                ),
                "SwathWindowIndex": selected_row["swath_window_index"],
                "SwathWindowGroup": selected_row["swath_window_group"],
                "SwathMzLower": selected_row["swath_mz_lower"],
                "SwathMzUpper": selected_row["swath_mz_upper"],
                "SwathIMLowerSQLite": selected_row["swath_im_lower"],
                "SwathIMUpperSQLite": selected_row["swath_im_upper"],
                "UsableFragments": selected_row["usable_fragments"],
            })
    truth_rows.sort(key=lambda row: (int(row["selection_rank"]), str(row["RunRole"])))
    write_tsv(args.truth_out, truth_rows)

    summary["selected_precursors"] = len(selected)
    summary["selected_proteins"] = len({row["protein_id"] for row in selected})
    summary["selected_charge_counts"] = dict(sorted(Counter(str(row["charge"]) for row in selected).items()))
    summary["selected_swath_map_counts"] = dict(sorted(Counter(str(row["swath_window_index"]) for row in selected).items(), key=lambda item: int(item[0])))
    summary["selected_distributions"] = {
        "min_realized_event_proxy": describe(row["min_realized_event_proxy"] for row in selected),
        "min_frame_abundance_sum": describe(row["min_frame_abundance_sum"] for row in selected),
        "min_scan_abundance_sum": describe(row["min_scan_abundance_sum"] for row in selected),
        "min_ion_relative_abundance": describe(row["min_ion_relative_abundance"] for row in selected),
        "assay_rt": describe(row["assay_rt"] for row in selected),
        "precursor_mz": describe(row["precursor_mz"] for row in selected),
        "assay_im": describe(row["assay_im"] for row in selected),
    }
    summary["status"] = "selected"
    args.qc_json.parent.mkdir(parents=True, exist_ok=True)
    args.qc_json.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")

    report_lines = [
        "# High-signal correctness fixture selection QC",
        "",
        "This selection used only TimSim truth, simulated fragment/library properties, and DIA acquisition geometry. ",
        "No OpenDIA identification, q-value, d-score, selected feature, or Percolator result was used.",
        "",
        f"- Selected precursor groups: **{len(selected)}/{args.precursors}**",
        f"- Eligible unique peptide precursors before stratification: **{len(eligible)}**",
        f"- Selected proteins: **{summary['selected_proteins']}**",
        f"- Minimum realized event proxy threshold: **{args.min_realized_event_proxy:g}** in every run",
        f"- Minimum frame abundance mass: **{args.min_frame_abundance_sum:g}** in every run",
        f"- Minimum scan abundance mass: **{args.min_scan_abundance_sum:g}** in every run",
        f"- Minimum selected charge fraction: **{args.min_ion_relative_abundance:g}** in every run",
        f"- RT edge margin: **max({100.0 * args.rt_edge_margin_fraction:g}% of acquisition duration, {args.rt_edge_min_seconds:g} s)** at each edge",
        "",
        "## Selected numeric distributions",
        "",
        "| Metric | Min | Q25 | Median | Q75 | Max |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for metric, stats in summary["selected_distributions"].items():
        report_lines.append(
            f"| {metric} | {stats['min']:.6g} | {stats['q25']:.6g} | {stats['median']:.6g} | {stats['q75']:.6g} | {stats['max']:.6g} |"
        )
    report_lines += [
        "",
        "## Charge distribution",
        "",
        "| Charge | Selected |",
        "|---:|---:|",
    ]
    for charge, count in summary["selected_charge_counts"].items():
        report_lines.append(f"| {charge} | {count} |")
    report_lines += [
        "",
        "## diaPASEF SWATH-map distribution",
        "",
        "| Map | Selected |",
        "|---:|---:|",
    ]
    for map_id, count in summary["selected_swath_map_counts"].items():
        report_lines.append(f"| {map_id} | {count} |")
    report_lines += [
        "",
        "The detailed long-form truth table contains assay coordinates and per-run simulator realization components. ",
        "SQLite-derived IM apex/centroid values are retained for provenance; final IM validation should prefer the raw TDF/XIPM signal coordinates from the oracle audit.",
        "",
    ]
    args.qc_report.parent.mkdir(parents=True, exist_ok=True)
    args.qc_report.write_text("\n".join(report_lines), encoding="utf-8")

    print(f"Selected exactly {len(selected)} high-signal precursor groups")
    print(f"Eligible unique peptide precursors before stratification: {len(eligible)}")
    print(f"Wrote selected precursors: {args.selected_out}")
    print(f"Wrote realized truth:      {args.truth_out}")
    print(f"Wrote selection QC:        {args.qc_report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
