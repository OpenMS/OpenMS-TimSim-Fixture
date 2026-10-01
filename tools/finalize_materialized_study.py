#!/usr/bin/env python3
"""Finalize and validate a scratch-materialized large TimSim study."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="", encoding="utf-8")
    return path.open(mode, newline="", encoding="utf-8")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path, "r") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def median(values: list[float]) -> float:
    return float(statistics.median(values)) if values else math.nan


def summarize_truth(path: Path, expected_rows: int, expected_run_name: str) -> tuple[dict[str, object], int]:
    biological_targets = 0
    observable_targets = 0
    biological_proteins: set[str] = set()
    observable_proteins: set[str] = set()
    row_count = 0
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_count += 1
            if row["RunName"] != expected_run_name:
                raise SystemExit(f"Run identity mismatch in {path}: {row['RunName']} != {expected_run_name}")
            if int(row["BiologicalPresentInRun"]) == 1:
                biological_targets += 1
                biological_proteins.add(row["ProteinId"])
            if int(row["ObservableInSimulation"]) == 1:
                observable_targets += 1
                observable_proteins.add(row["ProteinId"])
    if row_count != expected_rows:
        raise SystemExit(f"Run truth row count mismatch for {path}: {row_count} != {expected_rows}")
    return (
        {
            "biological_selected_targets": biological_targets,
            "observable_selected_targets": observable_targets,
            "biological_selected_proteins": len(biological_proteins),
            "observable_selected_proteins": len(observable_proteins),
        },
        row_count,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-root", type=Path, required=True)
    args = parser.parse_args()

    root = args.study_root.resolve()
    manifest_path = root / "OpenSwathTimSim.study_manifest.tsv"
    selected_path = root / "OpenSwathTimSim.high_signal_selection.tsv"
    design_path = root / "OpenSwathTimSim.study_design_truth.tsv"
    plan_path = root / "OpenSwathTimSim.materialized_study_plan.json"
    for path in (manifest_path, selected_path, design_path, plan_path):
        if not path.is_file():
            raise SystemExit(f"Required study artifact missing: {path}")

    manifest = read_tsv(manifest_path)
    selected = read_tsv(selected_path)
    if not manifest or not selected:
        raise SystemExit("Study manifest and target selection must be non-empty")

    tdfs = root / "tdfs"
    if tdfs.exists():
        for item in list(tdfs.iterdir()):
            if item.is_symlink() or item.is_file():
                item.unlink()
            elif item.is_dir():
                import shutil
                shutil.rmtree(item)
    tdfs.mkdir(parents=True, exist_ok=True)

    run_summaries: list[dict[str, object]] = []
    truth_index_rows: list[dict[str, object]] = []
    total_truth_rows = 0
    for run in sorted(manifest, key=lambda row: int(row["RunOrdinal"])):
        ordinal = int(run["RunOrdinal"])
        truth_path = root / "run_realized_truth" / f"{ordinal:03d}_{run['RunId']}.tsv.gz"
        stats_path = root / "run_stats" / f"{ordinal:03d}_{run['RunId']}.json"
        raw_dir = root / run["RunName"] / f"{run['RunName']}.d"
        if not truth_path.is_file():
            raise SystemExit(f"Missing run truth: {truth_path}")
        if not stats_path.is_file():
            raise SystemExit(f"Missing run stats: {stats_path}")
        if not (raw_dir / "analysis.tdf").is_file() or not (raw_dir / "analysis.tdf_bin").is_file():
            raise SystemExit(f"Incomplete generated .d for {run['RunName']}: {raw_dir}")

        truth_summary, row_count = summarize_truth(truth_path, len(selected), run["RunName"])
        total_truth_rows += row_count
        link = tdfs / f"{run['RunName']}.d"
        link.symlink_to(Path("..") / run["RunName"] / f"{run['RunName']}.d")

        stats = json.loads(stats_path.read_text(encoding="utf-8"))
        stats.update(truth_summary)
        run_summaries.append(stats)
        truth_index_rows.append({
            "RunOrdinal": ordinal,
            "RunId": run["RunId"],
            "RunName": run["RunName"],
            "Condition": run["Condition"],
            "Replicate": run["Replicate"],
            "InputLevel": run["InputLevel"],
            "TruthFile": str(truth_path.relative_to(root)),
            "Rows": row_count,
            "SHA256": sha256_file(truth_path),
            "BiologicalTargets": truth_summary["biological_selected_targets"],
            "ObservableTargets": truth_summary["observable_selected_targets"],
            "ObservableProteins": truth_summary["observable_selected_proteins"],
        })

    truth_index_path = root / "OpenSwathTimSim.realized_truth_index.tsv"
    with truth_index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(truth_index_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(truth_index_rows)

    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for item in run_summaries:
        grouped[(str(item["profile"]), str(item["input_level"]))].append(item)
    level_summary: list[dict[str, object]] = []
    for (profile, level), items in sorted(grouped.items()):
        level_summary.append({
            "profile": profile,
            "input_level": level,
            "runs": len(items),
            "median_input_log2_offset": median([float(item["input_log2_offset"]) for item in items]),
            "median_positive_blueprint_peptides": median([float(item["positive_peptides"]) for item in items]),
            "median_selected_targets_present": median([float(item["selected_present_targets"]) for item in items]),
            "median_selected_targets_observable": median([float(item["observable_selected_targets"]) for item in items]),
            "median_selected_proteins_observable": median([float(item["observable_selected_proteins"]) for item in items]),
        })

    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    summary = {
        "summary_version": "materialized_study_summary_v1",
        "profile": plan["profile"],
        "runs": len(manifest),
        "selected_precursors": len(selected),
        "selected_proteins": len({row["protein_id"] for row in selected}),
        "partitioned_truth_rows": total_truth_rows,
        "truth_layout": "gzip_partition_per_run",
        "materialization": "per_run_scratch_compacted",
        "level_summary": level_summary,
    }
    summary_path = root / "materialized_study_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    report = [
        "# Materialized TimSim study",
        "",
        f"- Profile: **{plan['profile']}**",
        f"- Runs: **{len(manifest)}**",
        f"- Frozen assay targets: **{len(selected)}**",
        f"- Frozen assay proteins: **{summary['selected_proteins']}**",
        f"- Partitioned realized-truth rows: **{total_truth_rows}**",
        "- Per-run source DBs: **materialized in scratch and not persisted**",
        "- Realized truth: **one gzip TSV partition per run plus a SHA-256 index**",
        "- Run truth distinguishes biological presence from simulator-observable precursor signal.",
        "",
        "## Input-level summary",
        "",
        "| Profile | Level | Runs | log2 input offset | Present targets | Observable targets | Observable proteins |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for item in level_summary:
        report.append(
            f"| {item['profile']} | {item['input_level']} | {item['runs']} | "
            f"{item['median_input_log2_offset']:.3f} | {item['median_selected_targets_present']:.0f} | "
            f"{item['median_selected_targets_observable']:.0f} | {item['median_selected_proteins_observable']:.0f} |"
        )
    (root / "materialized_study_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    artifact_candidates = {
        "study_manifest": manifest_path,
        "study_design_truth": design_path,
        "materialized_plan": plan_path,
        "selected_precursors": selected_path,
        "realized_truth_index": truth_index_path,
    }
    for name in (
        "OpenSwathTimSim.target_only.transitions.tsv",
        "OpenSwathTimSim.transitions.tsv",
        "OpenSwathTimSim.entrapment_truth.tsv",
        "OpenSwathTimSim.reference_truth.tsv",
    ):
        path = root / name
        if path.is_file():
            artifact_candidates[name] = path

    fixture_manifest = {
        "fixture_kind": "openms_timsim_materialized_large_study_v1",
        "profile": plan["profile"],
        "selection_frozen_before_condition_effects": True,
        "selection_frozen_before_opendia": True,
        "run_materialization": "scratch_local_then_compact",
        "runs": [row["RunName"] for row in manifest],
        "requested_precursors": len(selected),
        "selected_proteins": summary["selected_proteins"],
        "partitioned_truth": {
            "index_tsv": truth_index_path.name,
            "directory": "run_realized_truth",
            "compression": "gzip",
            "rows": total_truth_rows,
        },
        "truth_semantics": {
            "BiologicalPresentInRun": "RealizedInputEvents > 0",
            "ObservableInSimulation": "selected precursor has positive realized TimSim precursor signal",
        },
        "plan": plan,
        "artifact_sha256": {name: sha256_file(path) for name, path in artifact_candidates.items()},
    }
    (root / "fixture_manifest.json").write_text(json.dumps(fixture_manifest, indent=2) + "\n", encoding="utf-8")

    print(
        f"Materialized study OK: {len(manifest)} runs, {len(selected)} targets, "
        f"{summary['selected_proteins']} proteins, {total_truth_rows} partitioned truth rows"
    )
    print(root / "materialized_study_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
