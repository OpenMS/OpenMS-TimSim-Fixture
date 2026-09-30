#!/usr/bin/env python3
"""Validate the frozen simulator-only high-signal precursor selection."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

RUN_ROLES = {"baseline", "control_replicate", "treatment"}
FIXTURE_KIND = "openms_timsim_identification_quantification_entrapment"


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise SystemExit(f"Missing required TSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit(f"TSV has no rows: {path}")
    return rows


def f(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(f"Invalid numeric value for {key}: {row.get(key)!r}") from exc


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("build_dir", type=Path)
    args = parser.parse_args()
    root = args.build_dir.resolve()

    manifest_path = root / "fixture_manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("fixture_kind") != FIXTURE_KIND:
        raise SystemExit(
            "fixture_manifest.json has unexpected fixture_kind: "
            f"{manifest.get('fixture_kind')!r}; expected {FIXTURE_KIND!r}"
        )
    if manifest.get("selection_frozen_before_opendia") is not True:
        raise SystemExit("Manifest does not assert that selection was frozen before OpenDIA")

    requested = int(manifest["requested_precursors"])
    criteria = manifest["selection"]
    selected = read_tsv(root / criteria["selected_precursors_tsv"])
    truth = read_tsv(root / criteria["realized_truth_tsv"])

    if len(selected) != requested:
        raise SystemExit(f"Selected precursor count is {len(selected)}; expected exactly {requested}")
    precursor_keys = [row["precursor_key"] for row in selected]
    sequences = [row["sequence"] for row in selected]
    if len(set(precursor_keys)) != requested:
        raise SystemExit("Frozen selection contains duplicate precursor keys")
    if len(set(sequences)) != requested:
        raise SystemExit("Frozen selection contains more than one precursor charge per peptide sequence")

    expected_ranks = list(range(1, requested + 1))
    observed_ranks = sorted(int(row["selection_rank"]) for row in selected)
    if observed_ranks != expected_ranks:
        raise SystemExit("SelectionRank values are not exactly 1..N")

    truth_by_key: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in truth:
        truth_by_key[row["precursor_key"]].append(row)
    if set(truth_by_key) != set(precursor_keys):
        missing = sorted(set(precursor_keys) - set(truth_by_key))[:10]
        extra = sorted(set(truth_by_key) - set(precursor_keys))[:10]
        raise SystemExit(f"Truth/selection key mismatch; missing={missing}, extra={extra}")

    signal_floor = float(criteria["min_realized_event_proxy_all_runs"])
    frame_floor = float(criteria["min_frame_abundance_sum_all_runs"])
    scan_floor = float(criteria["min_scan_abundance_sum_all_runs"])
    ion_floor = float(criteria["min_ion_relative_abundance_all_runs"])
    edge_fraction = float(criteria["rt_edge_margin_fraction"])
    if "rt_edge_min_seconds" not in criteria:
        raise SystemExit(
            "Fixture predates the RT flank-support policy (missing selection.rt_edge_min_seconds); "
            "regenerate it with the current fixture before correctness/FDR validation"
        )
    edge_min_seconds = float(criteria["rt_edge_min_seconds"])

    failures: list[str] = []
    for key in precursor_keys:
        rows = truth_by_key[key]
        roles = {row["RunRole"] for row in rows}
        if len(rows) != 3 or roles != RUN_ROLES:
            failures.append(f"{key}: expected 3 run rows with roles {sorted(RUN_ROLES)}, observed {sorted(roles)}")
            continue
        for row in rows:
            if f(row, "RealizedEventProxy") < signal_floor:
                failures.append(f"{key}/{row['RunRole']}: realized signal below {signal_floor}")
            if f(row, "FrameAbundanceSum") < frame_floor:
                failures.append(f"{key}/{row['RunRole']}: frame mass below {frame_floor}")
            if f(row, "ScanAbundanceSum") < scan_floor:
                failures.append(f"{key}/{row['RunRole']}: scan mass below {scan_floor}")
            if f(row, "IonRelativeAbundance") < ion_floor:
                failures.append(f"{key}/{row['RunRole']}: ion fraction below {ion_floor}")
            start = f(row, "AcquisitionRTStart")
            end = f(row, "AcquisitionRTEnd")
            apex = f(row, "RealizedRTApex")
            margin = max((end - start) * edge_fraction, edge_min_seconds)
            if not (start + margin <= apex <= end - margin):
                failures.append(
                    f"{key}/{row['RunRole']}: RT apex {apex:.6g} outside interior [{start + margin:.6g}, {end - margin:.6g}]"
                )

    if failures:
        preview = "\n".join(f"  - {failure}" for failure in failures[:20])
        raise SystemExit(f"High-signal truth validation failed ({len(failures)} issues):\n{preview}")

    qc_json = root / criteria["qc_summary_json"]
    if not qc_json.is_file():
        raise SystemExit(f"Missing selection QC summary: {qc_json}")
    qc = json.loads(qc_json.read_text(encoding="utf-8"))
    if int(qc.get("selected_precursors", -1)) != requested or qc.get("status") != "selected":
        raise SystemExit("Selection QC summary does not confirm a complete selected set")
    if qc.get("selection_is_opendia_independent") is not True:
        raise SystemExit("Selection QC summary does not assert OpenDIA-independent selection")

    print(f"High-signal selection OK: {requested} precursor groups x 3 runs")
    print(
        "Thresholds: "
        f"realized_event_proxy >= {signal_floor:g}, frame mass >= {frame_floor:g}, "
        f"scan mass >= {scan_floor:g}, ion fraction >= {ion_floor:g}"
    )
    print(
        "RT edge exclusion: "
        f"max({100.0 * edge_fraction:g}% of acquisition duration, {edge_min_seconds:g} s) at each edge"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
