#!/usr/bin/env python3
"""Extract selected-target realized truth from one completed materialized TimSim run."""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path
from typing import Any

import select_high_signal_precursors as shared


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="", encoding="utf-8")
    return path.open(mode, newline="", encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "r") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"Refusing to write empty TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-db", type=Path, required=True)
    parser.add_argument("--selected-precursors", type=Path, required=True)
    parser.add_argument("--selected-input-truth", type=Path, required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    for path in (args.run_db, args.selected_precursors, args.selected_input_truth):
        if not path.is_file():
            raise SystemExit(f"Required input does not exist: {path}")

    selected_rows = read_tsv(args.selected_precursors)
    input_rows = read_tsv(args.selected_input_truth)
    input_by_key = {
        (row["PeptideSequence"].strip().upper(), int(row["PrecursorCharge"])): row
        for row in input_rows
    }
    truth = shared.load_run_truth(args.run_db, args.condition, args.run_name)

    output: list[dict[str, Any]] = []
    observable_count = 0
    biological_count = 0
    for selected in sorted(selected_rows, key=lambda row: int(row["selection_rank"])):
        key = (selected["sequence"].strip().upper(), int(selected["charge"]))
        input_row = input_by_key.get(key)
        if input_row is None:
            raise SystemExit(f"Missing selected input truth for {key}")
        biological_present = int(float(input_row["RealizedInputEvents"]) > 0)
        biological_count += biological_present
        run = truth.get(key)
        observable = int(run is not None and math.isfinite(run.realized_event_proxy) and run.realized_event_proxy > 0)
        observable_count += observable

        row: dict[str, Any] = {
            "selection_rank": selected["selection_rank"],
            "precursor_key": selected["precursor_key"],
            "TransitionGroupId": f"TIMSIM_{key[0]}_{key[1]}",
            "PeptideSequence": key[0],
            "PrecursorCharge": key[1],
            "ProteinId": selected["protein_id"],
            "RunOrdinal": input_row["RunOrdinal"],
            "RunId": input_row["RunId"],
            "RunName": input_row["RunName"],
            "Condition": input_row["Condition"],
            "Replicate": input_row["Replicate"],
            "InputLevel": input_row["InputLevel"],
            "AbundanceProfile": input_row["AbundanceProfile"],
            "TreatmentClass": input_row["TreatmentClass"],
            "DesignLog2FC": input_row["DesignLog2FC"],
            "BaselineEvents": input_row["BaselineEvents"],
            "ExpectedInputEvents": input_row["ExpectedInputEvents"],
            "RealizedInputEvents": input_row["RealizedInputEvents"],
            "TotalLog2FactorVsBlueprint": input_row["TotalLog2FactorVsBlueprint"],
            "BiologicalPresentInRun": biological_present,
            "ObservableInSimulation": observable,
        }
        if run is None:
            row.update({
                "PrecursorMz": selected["precursor_mz"],
                "AssayRT": selected["assay_rt"],
                "AssayIM": selected["assay_im"],
                "PeptideEvents": 0,
                "IonRelativeAbundance": 0.0,
                "FrameAbundanceSum": 0.0,
                "ScanAbundanceSum": 0.0,
                "RealizedEventProxy": 0.0,
                "RealizedRTApex": "nan",
                "RealizedRTCentroid": "nan",
                "RealizedIMApexSQLite": "nan",
                "RealizedIMCentroidSQLite": "nan",
                "FramePoints": 0,
                "ScanPoints": 0,
            })
        else:
            row.update({
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
        output.append(row)

    write_tsv(args.out, output)
    print(
        f"Run truth {args.run_name}: biological_present={biological_count}/{len(output)}, "
        f"observable={observable_count}/{len(output)}"
    )
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
