#!/usr/bin/env python3
"""Summarize correctness-fixture recovery/localization and enumerate every failure."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

RUN_ROLES = ("baseline", "control_replicate", "treatment")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--benchmark-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--oracle-audit", type=Path)
    parser.add_argument("--rt-threshold", type=float, default=1.0)
    parser.add_argument("--im-threshold", type=float, default=0.02)
    parser.add_argument("--require-perfect-recovery", action="store_true")
    args = parser.parse_args()

    selected = pd.read_csv(args.build_dir / "OpenSwathTimSim.high_signal_selection.tsv", sep="\t")
    bench = pd.read_csv(args.benchmark_dir / "precursor_benchmark.tsv", sep="\t", low_memory=False)
    expected = len(selected)
    if expected != 500:
        print(f"WARNING: correctness fixture contains {expected} targets rather than the production target of 500")

    work = bench.copy()
    work["detected"] = work["detected"].astype(str).str.lower().isin({"true", "1", "yes"})
    work["abs_rt_error"] = pd.to_numeric(work["rt_error_seconds"], errors="coerce").abs()
    work["abs_sqlite_im_error"] = pd.to_numeric(work["im_error"], errors="coerce").abs()

    if args.oracle_audit and args.oracle_audit.is_file():
        oracle = pd.read_csv(args.oracle_audit, sep="\t", low_memory=False)
        keep = [
            "run_name", "precursor_key", "raw_ms2_rt_centroid", "raw_ms2_im_centroid",
            "transitions_with_oracle_signal", "raw_ms2_oracle_intensity",
        ]
        available = [column for column in keep if column in oracle.columns]
        work = work.merge(oracle[available], on=["run_name", "precursor_key"], how="left", validate="one_to_one")
        if "raw_ms2_im_centroid" in work:
            work["tdf_im_error"] = pd.to_numeric(work["observed_im"], errors="coerce") - pd.to_numeric(work["raw_ms2_im_centroid"], errors="coerce")
            work["abs_tdf_im_error"] = work["tdf_im_error"].abs()
    else:
        work["abs_tdf_im_error"] = np.nan

    run_rows = []
    for role in RUN_ROLES:
        subset = work[work["run_role"] == role].copy()
        detected = subset[subset["detected"]]
        run_rows.append({
            "run_role": role,
            "targets": len(subset),
            "identified": int(subset["detected"].sum()),
            "recovery": float(subset["detected"].mean()),
            "rt_within_0_5s": int((detected["abs_rt_error"] <= 0.5).sum()),
            "rt_within_1_0s": int((detected["abs_rt_error"] <= 1.0).sum()),
            "rt_within_2_0s": int((detected["abs_rt_error"] <= 2.0).sum()),
            "identified_and_rt_correct": int((subset["detected"] & (subset["abs_rt_error"] <= args.rt_threshold)).sum()),
            "median_abs_rt_error": float(detected["abs_rt_error"].median()) if len(detected) else math.nan,
            "median_abs_sqlite_im_error": float(detected["abs_sqlite_im_error"].median()) if len(detected) else math.nan,
            "tdf_im_reference_available": int(detected["abs_tdf_im_error"].notna().sum()),
            "median_abs_tdf_im_error": float(detected["abs_tdf_im_error"].median()) if detected["abs_tdf_im_error"].notna().any() else math.nan,
            "identified_and_tdf_im_correct": int((subset["detected"] & (subset["abs_tdf_im_error"] <= args.im_threshold)).sum()) if subset["abs_tdf_im_error"].notna().any() else 0,
        })
    run_metrics = pd.DataFrame(run_rows)

    complete = work.pivot_table(index="precursor_key", columns="run_role", values="detected", aggfunc="max").reindex(columns=RUN_ROLES).fillna(False)
    complete_all = int(complete.all(axis=1).sum())

    failure_rows = []
    for row in work.itertuples(index=False):
        reasons = []
        if not bool(row.detected):
            reasons.append("not_identified")
        else:
            if math.isfinite(float(row.abs_rt_error)) and float(row.abs_rt_error) > args.rt_threshold:
                reasons.append(f"rt_error_gt_{args.rt_threshold:g}s")
            if hasattr(row, "abs_tdf_im_error") and math.isfinite(float(row.abs_tdf_im_error)) and float(row.abs_tdf_im_error) > args.im_threshold:
                reasons.append(f"tdf_im_error_gt_{args.im_threshold:g}")
        if reasons:
            failure_rows.append({
                "run_role": row.run_role,
                "run_name": row.run_name,
                "precursor_key": row.precursor_key,
                "sequence": row.sequence,
                "charge": row.charge,
                "failure_reasons": ";".join(reasons),
                "q_value": row.q_value,
                "intensity": row.intensity,
                "observed_rt": row.observed_rt,
                "realized_rt": row.evaluation_rt,
                "rt_error_seconds": row.rt_error_seconds,
                "observed_im": row.observed_im,
                "sqlite_realized_im": row.evaluation_im,
                "sqlite_im_error": row.im_error,
                "tdf_raw_ms2_im_centroid": getattr(row, "raw_ms2_im_centroid", math.nan),
                "tdf_im_error": getattr(row, "tdf_im_error", math.nan),
                "realized_event_proxy": row.realized_event_proxy,
                "frame_abundance_sum": row.frame_abundance_sum,
                "scan_abundance_sum": row.scan_abundance_sum,
                "ion_relative_abundance": row.ion_relative_abundance,
            })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    run_metrics.to_csv(args.out_dir / "correctness_run_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(failure_rows).to_csv(args.out_dir / "correctness_failures.tsv", sep="\t", index=False)

    summary = {
        "expected_targets": expected,
        "expected_run_target_observations": expected * len(RUN_ROLES),
        "complete_all_three_runs": complete_all,
        "perfect_recovery": complete_all == expected and all(int(row["identified"]) == expected for row in run_rows),
        "rt_correctness_threshold_seconds": args.rt_threshold,
        "tdf_im_correctness_threshold": args.im_threshold,
        "run_metrics": run_rows,
        "failure_rows": len(failure_rows),
        "oracle_tdf_im_reference_used": bool(args.oracle_audit and args.oracle_audit.is_file()),
    }
    (args.out_dir / "correctness_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")

    report = [
        "# OpenDIA high-signal correctness benchmark",
        "",
        f"Frozen simulator-only target set: **{expected} precursor groups x 3 runs**.",
        f"Complete across all three runs: **{complete_all}/{expected}**.",
        "",
        "| Run | Identified | Recovery | RT <=0.5 s | RT <=1.0 s | RT <=2.0 s |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in run_rows:
        report.append(
            f"| {row['run_role']} | {row['identified']}/{row['targets']} | {100*row['recovery']:.2f}% | "
            f"{row['rt_within_0_5s']} | {row['rt_within_1_0s']} | {row['rt_within_2_0s']} |"
        )
    report += [
        "",
        f"Every missed target or target exceeding the configured RT/IM localization thresholds is listed in `correctness_failures.tsv`; failures are never silently replaced.",
    ]
    if summary["oracle_tdf_im_reference_used"]:
        report.append("Raw TDF/XIPM MS2 mobility centroids from the oracle audit were used as the preferred IM localization reference where available.")
    else:
        report.append("No oracle audit was supplied, so IM values in this pass are SQLite-derived only. Run the oracle and repeat this postprocessor for TDF/XIPM-based IM validation.")
    (args.out_dir / "correctness_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print(run_metrics.to_string(index=False))
    print(f"Complete all runs: {complete_all}/{expected}")
    print(f"Failure rows: {len(failure_rows)}")
    if args.require_perfect_recovery and not summary["perfect_recovery"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
