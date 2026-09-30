#!/usr/bin/env python3
"""Build long-form realized precursor truth across all study runs."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import select_high_signal_precursors as shared


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"Refusing to write empty TSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", type=Path, required=True)
    parser.add_argument("--study-manifest", type=Path, required=True)
    parser.add_argument("--selected-precursors", type=Path, required=True)
    parser.add_argument("--abundance-truth", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    manifest_rows = read_tsv(args.study_manifest)
    selected_rows = read_tsv(args.selected_precursors)
    abundance_rows = read_tsv(args.abundance_truth)
    if not manifest_rows or not selected_rows:
        raise SystemExit("Study manifest and selected-precursor table must be non-empty")

    abundance = {
        (row["RunName"], row["PeptideSequence"].strip().upper()): row
        for row in abundance_rows
    }
    selected = {
        (row["sequence"].strip().upper(), int(row["charge"])): row
        for row in selected_rows
    }

    output: list[dict[str, Any]] = []
    for manifest in manifest_rows:
        run_name = manifest["RunName"]
        run_db = args.fixture_root / run_name / "synthetic_data.db"
        truth = shared.load_run_truth(run_db, manifest["Condition"], run_name)
        for key, selected_row in selected.items():
            if key not in truth:
                raise SystemExit(f"Selected precursor {key} missing from study run {run_name}")
            run = truth[key]
            abundance_row = abundance.get((run_name, run.sequence))
            if abundance_row is None:
                raise SystemExit(f"Missing abundance truth for {run_name} / {run.sequence}")
            output.append({
                "selection_rank": selected_row["selection_rank"],
                "precursor_key": selected_row["precursor_key"],
                "TransitionGroupId": f"TIMSIM_{run.sequence}_{run.charge}",
                "PeptideSequence": run.sequence,
                "PrecursorCharge": run.charge,
                "ProteinId": run.protein,
                "RunId": manifest["RunId"],
                "RunName": run_name,
                "Condition": manifest["Condition"],
                "Replicate": manifest["Replicate"],
                "TreatmentClass": abundance_row["TreatmentClass"],
                "DesignLog2FC": abundance_row["DesignLog2FC"],
                "BaselineEvents": abundance_row["BaselineEvents"],
                "RealizedInputEvents": abundance_row["RealizedInputEvents"],
                "TotalLog2FactorVsBlueprint": abundance_row["TotalLog2FactorVsBlueprint"],
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
                "FramePoints": run.frame_points,
                "ScanPoints": run.scan_points,
            })

    output.sort(key=lambda row: (int(row["selection_rank"]), int(next(m["RunOrdinal"] for m in manifest_rows if m["RunName"] == row["RunName"]))))
    write_tsv(args.out, output)
    print(
        f"Wrote realized truth: {args.out} "
        f"({len(selected)} precursors x {len(manifest_rows)} runs = {len(output)} rows)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
