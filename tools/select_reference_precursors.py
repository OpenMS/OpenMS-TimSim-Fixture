#!/usr/bin/env python3
"""Freeze a high-signal precursor set from a TimSim blueprint simulation only.

This is the canonical selector for multi-run studies. It deliberately does not
inspect condition-specific runs so strong treatment effects and biological
missingness cannot influence which targets enter the benchmark.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
import sys
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import select_high_signal_precursors as shared
from export_openswath_tsv import read_candidates


def _signal_sort_key(row: dict[str, Any]) -> tuple[float, float, float, int, str]:
    return (
        -float(row["min_realized_event_proxy"]),
        -float(row["min_ion_relative_abundance"]),
        -float(row["total_fragment_intensity"]),
        int(row["charge"]),
        str(row["precursor_key"]),
    )


def _build_fragment_collision_index(
    fragment_candidates: list[Any],
    windows: list[Any],
    mz_window_margin: float,
    im_window_margin: float,
) -> dict[int, tuple[list[float], list[tuple[float, float, float, str, float]]]]:
    by_group: dict[int, list[tuple[float, float, float, str, float]]] = defaultdict(list)
    for candidate in fragment_candidates:
        window = shared.assign_window(candidate, windows, mz_window_margin, im_window_margin)
        if window is None:
            continue
        key = f"{candidate.sequence}/{candidate.precursor_charge}"
        for fragment in candidate.fragments:
            by_group[int(window.window_group)].append(
                (
                    float(fragment["mz"]),
                    float(candidate.rt_seconds),
                    float(candidate.precursor_im),
                    key,
                    float(fragment["intensity"]),
                )
            )
    index: dict[int, tuple[list[float], list[tuple[float, float, float, str, float]]]] = {}
    for group, rows in by_group.items():
        rows.sort(key=lambda item: item[0])
        index[group] = ([item[0] for item in rows], rows)
    return index


def _specificity_metrics(
    candidate: Any,
    window_group: int,
    collision_index: dict[int, tuple[list[float], list[tuple[float, float, float, str, float]]]],
    *,
    max_rt_seconds: float,
    max_im_1_over_k0: float,
    max_fragment_ppm: float,
    top_fragments: int,
) -> dict[str, float | int]:
    target_key = f"{candidate.sequence}/{candidate.precursor_charge}"
    mz_values, background = collision_index.get(window_group, ([], []))
    fragment_rows: list[tuple[int, int, float, float, float]] = []
    for fragment in candidate.fragments:
        product_mz = float(fragment["mz"])
        delta = product_mz * max_fragment_ppm * 1e-6
        left = bisect_left(mz_values, product_mz - delta)
        right = bisect_right(mz_values, product_mz + delta)
        competitor_keys: set[str] = set()
        collision_count = 0
        collision_intensity = 0.0
        for mz, rt, im, precursor_key, predicted_intensity in background[left:right]:
            if precursor_key == target_key:
                continue
            if abs(rt - float(candidate.rt_seconds)) > max_rt_seconds:
                continue
            if abs(im - float(candidate.precursor_im)) > max_im_1_over_k0:
                continue
            collision_count += 1
            competitor_keys.add(precursor_key)
            collision_intensity += predicted_intensity
        fragment_rows.append(
            (
                len(competitor_keys),
                collision_count,
                collision_intensity,
                -float(fragment["intensity"]),
                product_mz,
            )
        )

    fragment_rows.sort()
    top = fragment_rows[: min(top_fragments, len(fragment_rows))]
    return {
        "specificity_collision_free_fragments": sum(row[0] == 0 for row in fragment_rows),
        "specificity_all_competitor_precursors_sum": sum(row[0] for row in fragment_rows),
        "specificity_all_collision_count_sum": sum(row[1] for row in fragment_rows),
        "specificity_all_collision_predicted_intensity": sum(row[2] for row in fragment_rows),
        "specificity_top_fragments": len(top),
        "specificity_top_competitor_precursors_sum": sum(row[0] for row in top),
        "specificity_top_collision_count_sum": sum(row[1] for row in top),
        "specificity_top_collision_predicted_intensity": sum(row[2] for row in top),
    }


def _annotate_collision_specificity(
    eligible: list[dict[str, Any]],
    fragment_candidates: list[Any],
    windows: list[Any],
    *,
    mz_window_margin: float,
    im_window_margin: float,
    max_rt_seconds: float,
    max_im_1_over_k0: float,
    max_fragment_ppm: float,
    top_fragments: int,
) -> None:
    candidate_by_key = {
        (str(candidate.sequence), int(candidate.precursor_charge)): candidate
        for candidate in fragment_candidates
    }
    collision_index = _build_fragment_collision_index(
        fragment_candidates, windows, mz_window_margin, im_window_margin
    )
    for row in eligible:
        key = (str(row["sequence"]), int(row["charge"]))
        candidate = candidate_by_key.get(key)
        if candidate is None:
            raise RuntimeError(f"Missing blueprint fragment candidate for {row['precursor_key']}")
        row.update(
            _specificity_metrics(
                candidate,
                int(row["swath_window_group"]),
                collision_index,
                max_rt_seconds=max_rt_seconds,
                max_im_1_over_k0=max_im_1_over_k0,
                max_fragment_ppm=max_fragment_ppm,
                top_fragments=top_fragments,
            )
        )



def _scale_variable_select(
    eligible: list[dict[str, Any]],
    count: int,
) -> list[dict[str, Any]]:
    """Select a deterministic abundance/geometry-stratified scale library.

    This selector is intentionally different from the correctness-tier high-signal
    selector. Every input row is already structurally valid and present in the
    blueprint. The realized-abundance marginal is preserved first; RT/mz/IM
    coverage is then preserved within each abundance quantile. A stable hash is
    used only to choose among otherwise equivalent candidates.
    """
    if count > len(eligible):
        raise ValueError(f"Cannot select {count} scale precursors from {len(eligible)} candidates")

    def geometry_stratum(row: dict[str, Any]) -> tuple[int, int, int]:
        return int(row["rt_bin"]), int(row["mz_bin"]), int(row["im_bin"])

    def stable_key(row: dict[str, Any]) -> tuple[int, str]:
        digest = hashlib.sha256(
            f"scale_variable_v1\0{row['precursor_key']}".encode("utf-8")
        ).digest()
        return int.from_bytes(digest[:8], "big"), str(row["precursor_key"])

    by_abundance: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in eligible:
        by_abundance[int(row["abundance_bin"])].append(row)
    abundance_targets = shared.proportional_targets(
        [int(row["abundance_bin"]) for row in eligible], count
    )

    selected: list[dict[str, Any]] = []
    global_remainder: list[dict[str, Any]] = []
    for abundance_bin in sorted(by_abundance):
        rows = by_abundance[abundance_bin]
        wanted = min(int(abundance_targets.get(abundance_bin, 0)), len(rows))
        by_geometry: dict[tuple[int, int, int], list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            by_geometry[geometry_stratum(row)].append(row)
        geometry_targets = shared.proportional_targets(
            [geometry_stratum(row) for row in rows], wanted
        ) if wanted else {}
        abundance_selected: list[dict[str, Any]] = []
        abundance_remainder: list[dict[str, Any]] = []
        for category in sorted(by_geometry):
            category_rows = sorted(by_geometry[category], key=stable_key)
            take = min(int(geometry_targets.get(category, 0)), len(category_rows))
            abundance_selected.extend(category_rows[:take])
            abundance_remainder.extend(category_rows[take:])
        if len(abundance_selected) < wanted:
            extra = sorted(abundance_remainder, key=stable_key)[: wanted - len(abundance_selected)]
            abundance_selected.extend(extra)
            used = {str(row["precursor_key"]) for row in extra}
            abundance_remainder = [
                row for row in abundance_remainder if str(row["precursor_key"]) not in used
            ]
        selected.extend(abundance_selected[:wanted])
        global_remainder.extend(abundance_remainder)

    if len(selected) < count:
        selected.extend(sorted(global_remainder, key=stable_key)[: count - len(selected)])
    elif len(selected) > count:
        selected = sorted(
            selected,
            key=lambda row: (int(row["abundance_bin"]), geometry_stratum(row), stable_key(row)),
        )[:count]

    selected = sorted(
        selected,
        key=lambda row: (int(row["abundance_bin"]), geometry_stratum(row), stable_key(row)),
    )
    result: list[dict[str, Any]] = []
    for rank, row in enumerate(selected, start=1):
        chosen = row.copy()
        chosen["selection_rank"] = rank
        result.append(chosen)
    return result

def _protein_precursor_shortlist(
    rows: list[dict[str, Any]],
    precursors_per_protein: int,
    shortlist_multiplier: int,
) -> list[dict[str, Any]]:
    signal_order = sorted(rows, key=_signal_sort_key)
    shortlist_size = min(
        len(signal_order),
        max(precursors_per_protein, precursors_per_protein * shortlist_multiplier),
    )
    shortlist = signal_order[:shortlist_size]
    specificity_order = sorted(
        shortlist,
        key=lambda row: (
            int(row["specificity_top_competitor_precursors_sum"]),
            int(row["specificity_all_competitor_precursors_sum"]),
            int(row["specificity_top_collision_count_sum"]),
            int(row["specificity_all_collision_count_sum"]),
            float(row["specificity_top_collision_predicted_intensity"]),
            float(row["specificity_all_collision_predicted_intensity"]),
            -int(row["specificity_collision_free_fragments"]),
            *_signal_sort_key(row),
        ),
    )
    chosen = [row.copy() for row in specificity_order[:precursors_per_protein]]
    for rank, row in enumerate(chosen, start=1):
        row["protein_precursor_rank"] = rank
    return chosen


def _balanced_protein_select(
    eligible: list[dict[str, Any]],
    *,
    target_proteins: int,
    precursors_per_protein: int,
    shortlist_multiplier: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    by_protein: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in eligible:
        by_protein[str(row["protein_id"])].append(row)

    eligible_per_protein = Counter(len(rows) for rows in by_protein.values())
    qualified = {
        protein: _protein_precursor_shortlist(rows, precursors_per_protein, shortlist_multiplier)
        for protein, rows in by_protein.items()
        if len(rows) >= precursors_per_protein
    }
    if len(qualified) < target_proteins:
        raise ValueError(
            f"Only {len(qualified)} proteins have at least {precursors_per_protein} eligible precursor groups; "
            f"{target_proteins} proteins are required"
        )

    dimensions = {
        "rt_bin": 3.0,
        "mz_bin": 2.0,
        "im_bin": 2.0,
        "swath_window_index": 2.5,
        "charge": 1.5,
    }
    all_unit_rows = [row for rows in qualified.values() for row in rows]
    total_precursors = target_proteins * precursors_per_protein
    targets = {
        dimension: shared.proportional_targets(
            [row[dimension] for row in all_unit_rows], total_precursors
        )
        for dimension in dimensions
    }
    selected_counts = {dimension: Counter() for dimension in dimensions}

    weakest_signal = {
        protein: min(math.log10(max(float(row["min_realized_event_proxy"]), 1.0)) for row in rows)
        for protein, rows in qualified.items()
    }
    signal_min = min(weakest_signal.values())
    signal_max = max(weakest_signal.values())
    signal_span = max(signal_max - signal_min, 1e-12)
    specificity_penalty = {
        protein: sum(int(row["specificity_top_competitor_precursors_sum"]) for row in rows)
        for protein, rows in qualified.items()
    }

    remaining = dict(qualified)
    selected_units: list[tuple[str, list[dict[str, Any]]]] = []
    for _ in range(target_proteins):
        best_protein: str | None = None
        best_tie: tuple[float, float, float] | None = None
        for protein, rows in remaining.items():
            coverage_score = 0.0
            for row in rows:
                for dimension, weight in dimensions.items():
                    category = row[dimension]
                    target = max(targets[dimension].get(category, 0), 1)
                    deficit = (target - selected_counts[dimension][category]) / target
                    coverage_score += weight * deficit
            coverage_score /= max(len(rows), 1)
            signal_quality = (weakest_signal[protein] - signal_min) / signal_span
            specificity_quality = 1.0 / (1.0 + specificity_penalty[protein])
            score = coverage_score + 0.15 * signal_quality + 0.10 * specificity_quality
            tie = (score, signal_quality, specificity_quality)
            if (
                best_tie is None
                or tie > best_tie
                or (tie == best_tie and protein < str(best_protein))
            ):
                best_protein = protein
                best_tie = tie
        assert best_protein is not None
        rows = remaining.pop(best_protein)
        selected_units.append((best_protein, rows))
        for row in rows:
            for dimension in dimensions:
                selected_counts[dimension][row[dimension]] += 1

    selected: list[dict[str, Any]] = []
    for protein_rank, (protein, rows) in enumerate(selected_units, start=1):
        for row in sorted(rows, key=lambda item: int(item["protein_precursor_rank"])):
            chosen = row.copy()
            chosen["protein_selection_rank"] = protein_rank
            chosen["selection_rank"] = len(selected) + 1
            selected.append(chosen)

    metadata = {
        "proteins_with_any_eligible_precursor": len(by_protein),
        "proteins_with_required_precursors": len(qualified),
        "eligible_precursors_per_protein_distribution": {
            str(count): proteins for count, proteins in sorted(eligible_per_protein.items())
        },
    }
    return selected, metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blueprint-db", type=Path, required=True)
    parser.add_argument("--selected-out", type=Path, required=True)
    parser.add_argument("--reference-truth-out", type=Path, required=True)
    parser.add_argument("--candidate-qc-out", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    parser.add_argument("--qc-report", type=Path, required=True)
    parser.add_argument("--precursors", type=int, default=1000)
    parser.add_argument(
        "--selection-mode",
        choices=("global_stratified", "protein_balanced", "scale_variable"),
        default="global_stratified",
    )
    parser.add_argument("--target-proteins", type=int, default=250)
    parser.add_argument("--precursors-per-protein", type=int, default=4)
    parser.add_argument("--signal-shortlist-multiplier", type=int, default=2)
    parser.add_argument("--specificity-max-rt-seconds", type=float, default=6.0)
    parser.add_argument("--specificity-max-im", type=float, default=0.03)
    parser.add_argument("--specificity-max-fragment-ppm", type=float, default=25.0)
    parser.add_argument("--specificity-top-fragments", type=int, default=4)
    parser.add_argument("--min-realized-event-proxy", type=float, default=50000.0)
    parser.add_argument("--min-frame-abundance-sum", type=float, default=0.90)
    parser.add_argument("--min-scan-abundance-sum", type=float, default=0.95)
    parser.add_argument("--min-ion-relative-abundance", type=float, default=0.25)
    parser.add_argument("--rt-edge-margin-fraction", type=float, default=0.05)
    parser.add_argument("--rt-edge-min-seconds", type=float, default=6.0)
    parser.add_argument("--minimum-fragments", type=int, default=8)
    parser.add_argument("--min-product-mz", type=float, default=350.0)
    parser.add_argument("--max-product-mz", type=float, default=2000.0)
    parser.add_argument("--mz-window-margin", type=float, default=1.0)
    parser.add_argument("--im-window-margin", type=float, default=0.005)
    parser.add_argument("--rt-bins", type=int, default=5)
    parser.add_argument("--scale-abundance-bins", type=int, default=10)
    parser.add_argument("--mz-bins", type=int, default=5)
    parser.add_argument("--im-bins", type=int, default=5)
    args = parser.parse_args()

    if args.precursors < 1:
        parser.error("--precursors must be positive")
    if args.target_proteins < 1:
        parser.error("--target-proteins must be positive")
    if args.precursors_per_protein < 1:
        parser.error("--precursors-per-protein must be positive")
    if args.signal_shortlist_multiplier < 1:
        parser.error("--signal-shortlist-multiplier must be positive")
    if args.scale_abundance_bins < 2:
        parser.error("--scale-abundance-bins must be >= 2")
    if args.specificity_top_fragments < 1:
        parser.error("--specificity-top-fragments must be positive")
    if any(value <= 0 for value in (args.specificity_max_rt_seconds, args.specificity_max_im, args.specificity_max_fragment_ppm)):
        parser.error("specificity geometry limits must be positive")
    if args.selection_mode == "protein_balanced":
        expected = args.target_proteins * args.precursors_per_protein
        if args.precursors != expected:
            parser.error(
                "protein_balanced selection requires --precursors == "
                "--target-proteins * --precursors-per-protein "
                f"({args.precursors} != {expected})"
            )
    if not args.blueprint_db.is_file():
        raise SystemExit(f"Blueprint DB does not exist: {args.blueprint_db}")
    if not 0 <= args.rt_edge_margin_fraction < 0.5:
        parser.error("--rt-edge-margin-fraction must be in [0, 0.5)")
    if args.rt_edge_min_seconds < 0:
        parser.error("--rt-edge-min-seconds must be non-negative")

    reference_truth = shared.load_run_truth(args.blueprint_db, "blueprint", "OpenSwathTimSim_blueprint")
    with sqlite3.connect(args.blueprint_db) as connection:
        windows = shared.load_windows_with_identity(connection)
        fragment_candidates = read_candidates(
            connection,
            min_product_mz=args.min_product_mz,
            max_product_mz=args.max_product_mz,
            minimum_fragments=args.minimum_fragments,
        )
        total_ions = int(connection.execute("SELECT COUNT(*) FROM ions").fetchone()[0])

    candidate_rows: list[dict[str, Any]] = []
    stages = Counter()
    for candidate in fragment_candidates:
        stages["sufficient_fragments"] += 1
        key = (candidate.sequence, candidate.precursor_charge)
        truth = reference_truth.get(key)
        window = shared.assign_window(candidate, windows, args.mz_window_margin, args.im_window_margin)
        geometry_ok = window is not None
        if geometry_ok:
            stages["safe_dia_geometry"] += 1
        present = truth is not None
        if geometry_ok and present:
            stages["present_in_blueprint"] += 1

        proxy = truth.realized_event_proxy if truth else math.nan
        frame_mass = truth.frame_abundance_sum if truth else math.nan
        scan_mass = truth.scan_abundance_sum if truth else math.nan
        ion_fraction = truth.ion_relative_abundance if truth else math.nan
        signal_ok = present and proxy >= args.min_realized_event_proxy
        frame_ok = present and frame_mass >= args.min_frame_abundance_sum
        scan_ok = present and scan_mass >= args.min_scan_abundance_sum
        ion_ok = present and ion_fraction >= args.min_ion_relative_abundance
        margin = math.nan
        rt_edge_ok = False
        if truth:
            duration = truth.acquisition_rt_end - truth.acquisition_rt_start
            margin = max(duration * args.rt_edge_margin_fraction, args.rt_edge_min_seconds)
            rt_edge_ok = (
                math.isfinite(truth.realized_rt_apex)
                and truth.acquisition_rt_start + margin <= truth.realized_rt_apex <= truth.acquisition_rt_end - margin
            )

        strict_high_signal_eligible = (
            geometry_ok and present and signal_ok and frame_ok and scan_ok and ion_ok and rt_edge_ok
        )
        scale_structural_eligible = (
            geometry_ok
            and present
            and math.isfinite(proxy) and proxy > 0
            and math.isfinite(frame_mass) and frame_mass > 0
            and math.isfinite(scan_mass) and scan_mass > 0
            and math.isfinite(ion_fraction) and ion_fraction > 0
            and rt_edge_ok
        )
        eligible = (
            scale_structural_eligible
            if args.selection_mode == "scale_variable"
            else strict_high_signal_eligible
        )
        if strict_high_signal_eligible:
            stages["strict_high_signal_eligible"] += 1
        if scale_structural_eligible:
            stages["scale_structural_eligible"] += 1
        if eligible:
            stages["eligible"] += 1
        candidate_rows.append({
            "precursor_key": f"{candidate.sequence}/{candidate.precursor_charge}",
            "sequence": candidate.sequence,
            "protein_id": candidate.protein,
            "charge": candidate.precursor_charge,
            "precursor_mz": candidate.precursor_mz,
            "assay_rt": candidate.rt_seconds,
            "assay_im": candidate.precursor_im,
            "usable_fragments": len(candidate.fragments),
            "total_fragment_intensity": candidate.total_fragment_intensity,
            # Keep these field names aligned with the existing stratified selector.
            "min_realized_event_proxy": proxy,
            "min_frame_abundance_sum": frame_mass,
            "min_scan_abundance_sum": scan_mass,
            "min_ion_relative_abundance": ion_fraction,
            "realized_rt_apex_baseline": truth.realized_rt_apex if truth else math.nan,
            "realized_im_apex_baseline": truth.realized_im_apex if truth else math.nan,
            "rt_edge_margin_seconds": margin,
            "swath_window_index": window.window_index if window else "",
            "swath_window_group": window.window_group if window else "",
            "swath_window_label": window.label if window else "",
            "swath_mz_lower": window.mz_lower if window else math.nan,
            "swath_mz_upper": window.mz_upper if window else math.nan,
            "swath_im_lower": window.im_lower if window else math.nan,
            "swath_im_upper": window.im_upper if window else math.nan,
            "safe_dia_geometry": int(geometry_ok),
            "present_in_blueprint": int(present),
            "signal_ok": int(signal_ok),
            "frame_ok": int(frame_ok),
            "scan_ok": int(scan_ok),
            "ion_fraction_ok": int(ion_ok),
            "rt_edge_ok": int(rt_edge_ok),
            "strict_high_signal_eligible": int(strict_high_signal_eligible),
            "scale_structural_eligible": int(scale_structural_eligible),
            "eligible": int(eligible),
        })

    shared.write_tsv(args.candidate_qc_out, candidate_rows)
    strict_high_signal_unique = {
        str(row["sequence"]) for row in candidate_rows if row["strict_high_signal_eligible"]
    }
    scale_structural_unique = {
        str(row["sequence"]) for row in candidate_rows if row["scale_structural_eligible"]
    }
    eligible_ions = [row for row in candidate_rows if row["eligible"]]
    best_by_sequence: dict[str, dict[str, Any]] = {}
    for row in eligible_ions:
        ranking = (
            float(row["min_realized_event_proxy"]),
            float(row["min_ion_relative_abundance"]),
            float(row["total_fragment_intensity"]),
            -int(row["charge"]),
        )
        previous = best_by_sequence.get(row["sequence"])
        if previous is None:
            best_by_sequence[row["sequence"]] = row.copy()
            continue
        previous_ranking = (
            float(previous["min_realized_event_proxy"]),
            float(previous["min_ion_relative_abundance"]),
            float(previous["total_fragment_intensity"]),
            -int(previous["charge"]),
        )
        if ranking > previous_ranking:
            best_by_sequence[row["sequence"]] = row.copy()

    eligible = list(best_by_sequence.values())
    if eligible:
        shared.assign_rank_bins(eligible, "assay_rt", "rt_bin", args.rt_bins)
        shared.assign_rank_bins(eligible, "precursor_mz", "mz_bin", args.mz_bins)
        shared.assign_rank_bins(eligible, "assay_im", "im_bin", args.im_bins)
        for row in eligible:
            row["log10_realized_event_proxy"] = math.log10(
                max(float(row["min_realized_event_proxy"]), 1e-300)
            )
        if args.selection_mode == "scale_variable":
            shared.assign_rank_bins(
                eligible,
                "log10_realized_event_proxy",
                "abundance_bin",
                args.scale_abundance_bins,
            )

    summary = {
        "selection_is_opendia_independent": True,
        "selection_scope": "blueprint_only_before_condition_effects",
        "selection_mode": args.selection_mode,
        "requested_precursors": args.precursors,
        "requested_target_proteins": args.target_proteins if args.selection_mode == "protein_balanced" else None,
        "precursors_per_protein": args.precursors_per_protein if args.selection_mode == "protein_balanced" else None,
        "total_ions_in_blueprint_db": total_ions,
        "fragment_eligible_ions": len(fragment_candidates),
        "eligible_ions_after_all_filters": len(eligible_ions),
        "eligible_unique_peptides_after_best_charge": len(eligible),
        "strict_high_signal_unique_peptides": len(strict_high_signal_unique),
        "scale_structural_unique_peptides": len(scale_structural_unique),
        "criteria": {
            "min_realized_event_proxy_blueprint": args.min_realized_event_proxy,
            "min_frame_abundance_sum_blueprint": args.min_frame_abundance_sum,
            "min_scan_abundance_sum_blueprint": args.min_scan_abundance_sum,
            "min_ion_relative_abundance_blueprint": args.min_ion_relative_abundance,
            "rt_edge_margin_fraction_of_acquisition": args.rt_edge_margin_fraction,
            "rt_edge_min_seconds": args.rt_edge_min_seconds,
            "minimum_usable_fragments": args.minimum_fragments,
            "mz_window_margin_th": args.mz_window_margin,
            "im_window_margin_1_over_k0": args.im_window_margin,
        },
        "cumulative_filter_counts": dict(stages),
        "selected_precursors": 0,
        "status": "ready",
    }

    failure_message: str | None = None
    if len(eligible) < args.precursors:
        summary["status"] = "insufficient_candidates"
        if args.selection_mode == "scale_variable":
            failure_message = (
                f"Only {len(eligible)} unique structurally valid precursor groups satisfy the scale blueprint contract; "
                f"{args.precursors} are required. Increase the candidate population; do not relax fragment/geometry safety criteria merely to reach the target count."
            )
        else:
            failure_message = (
                f"Only {len(eligible)} unique high-signal precursor groups satisfy blueprint criteria; "
                f"{args.precursors} are required. Increase candidate population; do not relax thresholds merely to reach the target count."
            )

    selected: list[dict[str, Any]]
    if failure_message is None and args.selection_mode == "protein_balanced":
        eligible_counts = Counter(str(row["protein_id"]) for row in eligible)
        qualified_proteins = sum(count >= args.precursors_per_protein for count in eligible_counts.values())
        summary["proteins_with_any_eligible_precursor"] = len(eligible_counts)
        summary["proteins_with_required_precursors"] = qualified_proteins
        summary["specificity"] = {
            "source": "complete_TimSim_blueprint_fragment_geometry",
            "uses_opendia": False,
            "uses_observed_raw_intensity": False,
            "max_rt_seconds": args.specificity_max_rt_seconds,
            "max_im_1_over_k0": args.specificity_max_im,
            "max_fragment_ppm": args.specificity_max_fragment_ppm,
            "top_fragments": args.specificity_top_fragments,
            "signal_shortlist_multiplier": args.signal_shortlist_multiplier,
        }
        if qualified_proteins < args.target_proteins:
            summary["status"] = "insufficient_protein_redundancy"
            failure_message = (
                f"Only {qualified_proteins} proteins have at least {args.precursors_per_protein} eligible precursor groups; "
                f"{args.target_proteins} proteins are required for the production contract. "
                "Increase the synthetic candidate population; do not relax the high-signal criteria."
            )
        else:
            _annotate_collision_specificity(
                eligible,
                fragment_candidates,
                windows,
                mz_window_margin=args.mz_window_margin,
                im_window_margin=args.im_window_margin,
                max_rt_seconds=args.specificity_max_rt_seconds,
                max_im_1_over_k0=args.specificity_max_im,
                max_fragment_ppm=args.specificity_max_fragment_ppm,
                top_fragments=args.specificity_top_fragments,
            )
            selected, protein_metadata = _balanced_protein_select(
                eligible,
                target_proteins=args.target_proteins,
                precursors_per_protein=args.precursors_per_protein,
                shortlist_multiplier=args.signal_shortlist_multiplier,
            )
            summary.update(protein_metadata)
    elif failure_message is None and args.selection_mode == "scale_variable":
        selected = _scale_variable_select(eligible, args.precursors)
    elif failure_message is None:
        selected = shared.stratified_select(eligible, args.precursors)
    else:
        selected = []

    if failure_message is not None:
        args.qc_json.parent.mkdir(parents=True, exist_ok=True)
        args.qc_json.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
        args.qc_report.parent.mkdir(parents=True, exist_ok=True)
        args.qc_report.write_text(
            "# Blueprint precursor-selection QC\n\n"
            f"Selection failed without weakening criteria: {failure_message}\n",
            encoding="utf-8",
        )
        raise SystemExit(failure_message)

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
        "rt_bin", "mz_bin", "im_bin", "abundance_bin", "log10_realized_event_proxy",
        "protein_selection_rank", "protein_precursor_rank",
        "specificity_collision_free_fragments",
        "specificity_all_competitor_precursors_sum",
        "specificity_all_collision_count_sum",
        "specificity_all_collision_predicted_intensity",
        "specificity_top_fragments",
        "specificity_top_competitor_precursors_sum",
        "specificity_top_collision_count_sum",
        "specificity_top_collision_predicted_intensity",
    ]
    shared.write_tsv(args.selected_out, selected, selected_fields)

    reference_rows: list[dict[str, Any]] = []
    for selected_row in selected:
        key = (str(selected_row["sequence"]), int(selected_row["charge"]))
        run = reference_truth[key]
        reference_rows.append({
            "selection_rank": selected_row["selection_rank"],
            "precursor_key": selected_row["precursor_key"],
            "TransitionGroupId": f"TIMSIM_{run.sequence}_{run.charge}",
            "PeptideSequence": run.sequence,
            "PrecursorCharge": run.charge,
            "ProteinId": run.protein,
            "PrecursorMz": run.precursor_mz,
            "AssayRT": run.assay_rt,
            "AssayIM": run.assay_im,
            "PeptideEvents": run.peptide_events,
            "IonRelativeAbundance": run.ion_relative_abundance,
            "FrameAbundanceSum": run.frame_abundance_sum,
            "ScanAbundanceSum": run.scan_abundance_sum,
            "RealizedEventProxy": run.realized_event_proxy,
            "RealizedRTApex": run.realized_rt_apex,
            "RealizedRTCentroid": run.realized_rt_centroid,
            "RealizedIMApexSQLite": run.realized_im_apex,
            "RealizedIMCentroidSQLite": run.realized_im_centroid,
            "RTEdgeMarginSeconds": max(
                (run.acquisition_rt_end - run.acquisition_rt_start) * args.rt_edge_margin_fraction,
                args.rt_edge_min_seconds,
            ),
            "SwathWindowIndex": selected_row["swath_window_index"],
            "SwathWindowGroup": selected_row["swath_window_group"],
        })
    shared.write_tsv(args.reference_truth_out, reference_rows)

    summary["selected_precursors"] = len(selected)
    summary["selected_proteins"] = len({row["protein_id"] for row in selected})
    summary["selected_precursors_per_protein"] = dict(
        sorted(Counter(Counter(str(row["protein_id"]) for row in selected).values()).items())
    )
    summary["selected_charge_counts"] = dict(sorted(Counter(str(row["charge"]) for row in selected).items()))
    summary["selected_distributions"] = {
        "reference_realized_event_proxy": shared.describe(row["min_realized_event_proxy"] for row in selected),
        "reference_frame_abundance_sum": shared.describe(row["min_frame_abundance_sum"] for row in selected),
        "reference_scan_abundance_sum": shared.describe(row["min_scan_abundance_sum"] for row in selected),
        "reference_ion_relative_abundance": shared.describe(row["min_ion_relative_abundance"] for row in selected),
        "assay_rt": shared.describe(row["assay_rt"] for row in selected),
        "precursor_mz": shared.describe(row["precursor_mz"] for row in selected),
        "assay_im": shared.describe(row["assay_im"] for row in selected),
    }
    summary["status"] = "selected"
    args.qc_json.parent.mkdir(parents=True, exist_ok=True)
    args.qc_json.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")

    lines = [
        "# Blueprint precursor-selection QC",
        "",
        "Targets were frozen from the TimSim blueprint before control/treatment replicate effects were created.",
        "OpenDIA results and condition-specific realized abundance were not used for selection.",
        "",
        f"- Selection mode: **{args.selection_mode}**",
        f"- Selected precursor groups: **{len(selected)}/{args.precursors}**",
        f"- Eligible unique precursor groups before selection: **{len(eligible)}**",
        f"- Selected proteins: **{summary['selected_proteins']}**",
        *(
            [
                f"- Production composition: **{args.target_proteins} proteins × {args.precursors_per_protein} precursors/protein**",
                f"- Proteins with at least {args.precursors_per_protein} eligible precursors: **{summary['proteins_with_required_precursors']}**",
                "- Fragment specificity: **simulator-only blueprint collision geometry**",
            ]
            if args.selection_mode == "protein_balanced"
            else []
        ),
        *(
            [
                "- Scale contract: **structurally valid, blueprint-present precursors sampled across realized abundance**",
                f"- Abundance quantile bins: **{args.scale_abundance_bins}**",
                f"- Strict correctness-tier high-signal candidates in the same blueprint: **{len(strict_high_signal_unique)}**",
            ]
            if args.selection_mode == "scale_variable"
            else [f"- Minimum blueprint realized event proxy: **{args.min_realized_event_proxy:g}**"]
        ),
        "",
        "Condition-specific missingness is intentionally allowed after selection and is part of the benchmark.",
        "",
    ]
    args.qc_report.parent.mkdir(parents=True, exist_ok=True)
    args.qc_report.write_text("\n".join(lines), encoding="utf-8")
    if args.selection_mode == "scale_variable":
        print(f"Selected exactly {len(selected)} blueprint scale-variable precursor groups")
    else:
        print(f"Selected exactly {len(selected)} blueprint high-signal precursor groups")
    if args.selection_mode == "protein_balanced":
        print(
            f"Production composition: {args.target_proteins} proteins x "
            f"{args.precursors_per_protein} precursors/protein"
        )
    print(f"Eligible unique peptide precursors before selection: {len(eligible)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
