#!/usr/bin/env python3
"""Freeze a high-signal precursor set from a TimSim blueprint simulation only.

This is the canonical selector for multi-run studies. It deliberately does not
inspect condition-specific runs so strong treatment effects and biological
missingness cannot influence which targets enter the benchmark.
"""
from __future__ import annotations

import argparse
import json
import math
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any

import select_high_signal_precursors as shared
from export_openswath_tsv import read_candidates


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blueprint-db", type=Path, required=True)
    parser.add_argument("--selected-out", type=Path, required=True)
    parser.add_argument("--reference-truth-out", type=Path, required=True)
    parser.add_argument("--candidate-qc-out", type=Path, required=True)
    parser.add_argument("--qc-json", type=Path, required=True)
    parser.add_argument("--qc-report", type=Path, required=True)
    parser.add_argument("--precursors", type=int, default=1000)
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
    parser.add_argument("--mz-bins", type=int, default=5)
    parser.add_argument("--im-bins", type=int, default=5)
    args = parser.parse_args()

    if args.precursors < 1:
        parser.error("--precursors must be positive")
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

        eligible = geometry_ok and present and signal_ok and frame_ok and scan_ok and ion_ok and rt_edge_ok
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
            "eligible": int(eligible),
        })

    shared.write_tsv(args.candidate_qc_out, candidate_rows)
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

    summary = {
        "selection_is_opendia_independent": True,
        "selection_scope": "blueprint_only_before_condition_effects",
        "requested_precursors": args.precursors,
        "total_ions_in_blueprint_db": total_ions,
        "fragment_eligible_ions": len(fragment_candidates),
        "eligible_ions_after_all_filters": len(eligible_ions),
        "eligible_unique_peptides_after_best_charge": len(eligible),
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
        "status": "insufficient_candidates" if len(eligible) < args.precursors else "ready",
    }
    if len(eligible) < args.precursors:
        args.qc_json.parent.mkdir(parents=True, exist_ok=True)
        args.qc_json.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
        args.qc_report.parent.mkdir(parents=True, exist_ok=True)
        args.qc_report.write_text(
            "# Blueprint precursor-selection QC\n\n"
            f"Selection failed without weakening criteria: **{len(eligible)}** unique precursor groups qualified; "
            f"**{args.precursors}** are required. Increase the synthetic candidate population.\n",
            encoding="utf-8",
        )
        raise SystemExit(
            f"Only {len(eligible)} unique high-signal precursor groups satisfy blueprint criteria; "
            f"{args.precursors} are required. Increase candidate population; do not relax thresholds merely to reach the target count."
        )

    selected = shared.stratified_select(eligible, args.precursors)
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
        f"- Selected precursor groups: **{len(selected)}/{args.precursors}**",
        f"- Eligible unique precursor groups before stratification: **{len(eligible)}**",
        f"- Selected proteins: **{summary['selected_proteins']}**",
        f"- Minimum blueprint realized event proxy: **{args.min_realized_event_proxy:g}**",
        "",
        "Condition-specific missingness is intentionally allowed after selection and is part of the benchmark.",
        "",
    ]
    args.qc_report.parent.mkdir(parents=True, exist_ok=True)
    args.qc_report.write_text("\n".join(lines), encoding="utf-8")
    print(f"Selected exactly {len(selected)} blueprint high-signal precursor groups")
    print(f"Eligible unique peptide precursors before stratification: {len(eligible)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
