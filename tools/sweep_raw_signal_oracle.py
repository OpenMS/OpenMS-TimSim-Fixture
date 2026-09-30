#!/usr/bin/env python3
"""Sweep truth-centered raw-oracle windows and simulator-only transition specificity.

The sweep is diagnostic only. It scans each generated DIA-PASEF run once, accumulates
all declared RT/IM/mass-window combinations, and evaluates transition subsets ranked
by collision geometry in the complete TimSim blueprint fragment universe. OpenDIA
coordinates, scores, peak boundaries, or quantitative performance are never used to
rank transitions or choose extraction parameters.
"""
from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

import benchmark_raw_signal_oracle as benchmark
import extract_raw_signal_oracle as oracle
import select_high_signal_precursors as shared
from export_openswath_tsv import read_candidates

DEFAULT_RT_WINDOWS = (1.0, 2.0, 4.0, 6.0)
DEFAULT_IM_WINDOWS = (0.01, 0.02, 0.03)
DEFAULT_PPM_WINDOWS = (10.0, 15.0, 25.0)
DEFAULT_SUBSETS = (8, 6, 4, 3)


def parse_float_grid(text: str, label: str) -> tuple[float, ...]:
    values = tuple(sorted({float(token.strip()) for token in text.split(",") if token.strip()}))
    if not values or any(not math.isfinite(value) or value <= 0 for value in values):
        raise ValueError(f"{label} must contain positive finite comma-separated values")
    return values


def parse_int_grid(text: str, label: str) -> tuple[int, ...]:
    values = tuple(sorted({int(token.strip()) for token in text.split(",") if token.strip()}, reverse=True))
    if not values or any(value <= 0 for value in values):
        raise ValueError(f"{label} must contain positive comma-separated integers")
    return values


def load_library_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"PeptideSequence", "PrecursorCharge", "ProductMz", "TransitionGroupId"}
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(f"Transition library missing columns: {sorted(missing)}")
    if "Decoy" in table.columns:
        table = table[pd.to_numeric(table["Decoy"], errors="coerce").fillna(0).eq(0)].copy()
    if "QuantifyingTransition" in table.columns:
        quant = pd.to_numeric(table["QuantifyingTransition"], errors="coerce").fillna(0).ne(0)
        table = table[quant].copy()
    table["PeptideSequence"] = table["PeptideSequence"].astype(str).str.strip().str.upper()
    table["PrecursorCharge"] = pd.to_numeric(table["PrecursorCharge"], errors="raise").astype(int)
    table["ProductMz"] = pd.to_numeric(table["ProductMz"], errors="raise")
    if "LibraryIntensity" not in table.columns:
        table["LibraryIntensity"] = 1.0
    table["LibraryIntensity"] = pd.to_numeric(table["LibraryIntensity"], errors="coerce").fillna(0.0)
    table["TransitionOrdinal"] = table.groupby(
        ["PeptideSequence", "PrecursorCharge"], sort=False
    ).cumcount() + 1
    return table.reset_index(drop=True)


def load_reference_coordinates(study: Path) -> pd.DataFrame:
    path = study / "OpenSwathTimSim.reference_truth.tsv"
    if path.is_file():
        table = pd.read_csv(path, sep="\t", low_memory=False)
        im_col = "RealizedIMApexSQLite" if "RealizedIMApexSQLite" in table.columns else "AssayIM"
        rt_col = "RealizedRTApex" if "RealizedRTApex" in table.columns else "AssayRT"
        result = table[["PeptideSequence", "PrecursorCharge", "PrecursorMz", rt_col, im_col]].copy()
        result.columns = ["PeptideSequence", "PrecursorCharge", "PrecursorMz", "ReferenceRT", "ReferenceIM"]
    else:
        realized = pd.read_csv(study / "OpenSwathTimSim.realized_truth.tsv", sep="\t", low_memory=False)
        first_run = realized.sort_values("RunName").drop_duplicates(["PeptideSequence", "PrecursorCharge"])
        result = first_run[
            ["PeptideSequence", "PrecursorCharge", "PrecursorMz", "RealizedRTApex", "RealizedIMApexSQLite"]
        ].copy()
        result.columns = ["PeptideSequence", "PrecursorCharge", "PrecursorMz", "ReferenceRT", "ReferenceIM"]
    result["PeptideSequence"] = result["PeptideSequence"].astype(str).str.upper()
    result["PrecursorCharge"] = pd.to_numeric(result["PrecursorCharge"], errors="raise").astype(int)
    return result


def blueprint_db_path(study: Path) -> Path:
    manifest_path = study / "fixture_manifest.json"
    blueprint_name = "OpenSwathTimSim_blueprint"
    if manifest_path.is_file():
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        blueprint_name = str(payload.get("blueprint_run") or blueprint_name)
    path = study / blueprint_name / "synthetic_data.db"
    if not path.is_file():
        raise FileNotFoundError(f"Missing TimSim blueprint database: {path}")
    return path


def build_background_fragments(
    blueprint_db: Path,
    *,
    min_product_mz: float,
    max_product_mz: float,
) -> pd.DataFrame:
    with sqlite3.connect(blueprint_db) as connection:
        windows = shared.load_windows_with_identity(connection)
        candidates = read_candidates(
            connection,
            min_product_mz=min_product_mz,
            max_product_mz=max_product_mz,
            minimum_fragments=1,
        )
    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        window_group = oracle.assign_window_group(candidate.precursor_mz, candidate.precursor_im, windows)
        if window_group is None:
            continue
        precursor_key = f"{candidate.sequence}/{candidate.precursor_charge}"
        for fragment in candidate.fragments:
            rows.append(
                {
                    "PrecursorKey": precursor_key,
                    "WindowGroup": int(window_group),
                    "RT": float(candidate.rt_seconds),
                    "IM": float(candidate.precursor_im),
                    "ProductMz": float(fragment["mz"]),
                    "PredictedIntensity": float(fragment["intensity"]),
                }
            )
    if not rows:
        raise RuntimeError("No blueprint fragment background could be constructed")
    return pd.DataFrame(rows)


def collision_metrics_for_transition(
    product_mz: float,
    target_key: str,
    target_rt: float,
    target_im: float,
    target_group: int,
    group_background: pd.DataFrame,
    *,
    max_rt: float,
    max_im: float,
    max_ppm: float,
) -> tuple[int, int, float, float]:
    # group_background is pre-sorted once by ProductMz in rank_transition_specificity.
    sorted_mz = pd.to_numeric(group_background["ProductMz"], errors="coerce").to_numpy(dtype=float)
    delta = product_mz * max_ppm * 1e-6
    left = int(np.searchsorted(sorted_mz, product_mz - delta, side="left"))
    right = int(np.searchsorted(sorted_mz, product_mz + delta, side="right"))
    if right <= left:
        return 0, 0, 0.0, math.inf
    subset = group_background.iloc[left:right].copy()
    subset = subset[
        subset["PrecursorKey"].astype(str).ne(target_key)
        & (pd.to_numeric(subset["RT"], errors="coerce").sub(target_rt).abs() <= max_rt)
        & (pd.to_numeric(subset["IM"], errors="coerce").sub(target_im).abs() <= max_im)
    ]
    if subset.empty:
        return 0, 0, 0.0, math.inf
    ppm_error = (
        pd.to_numeric(subset["ProductMz"], errors="coerce").sub(product_mz).abs() / product_mz * 1e6
    )
    return (
        int(len(subset)),
        int(subset["PrecursorKey"].nunique()),
        float(pd.to_numeric(subset["PredictedIntensity"], errors="coerce").fillna(0.0).sum()),
        float(ppm_error.min()),
    )


def merge_library_reference_coordinates(
    library: pd.DataFrame, reference: pd.DataFrame
) -> pd.DataFrame:
    """Attach canonical reference precursor coordinates without pandas suffix collisions."""
    keys = ["PeptideSequence", "PrecursorCharge"]
    left = library.copy()
    if "PrecursorMz" in left.columns:
        left = left.rename(columns={"PrecursorMz": "LibraryPrecursorMz"})
    merged = left.merge(
        reference,
        on=keys,
        how="left",
        validate="many_to_one",
    )
    required = ["PrecursorMz", "ReferenceRT", "ReferenceIM"]
    missing = [column for column in required if column not in merged.columns]
    if missing:
        raise RuntimeError(
            "Reference-coordinate merge did not produce required columns: "
            + ", ".join(missing)
        )
    if merged[required].isna().any().any():
        raise RuntimeError("Could not map all target transitions to reference RT/IM coordinates")
    return merged


def rank_transition_specificity(
    study: Path,
    library: pd.DataFrame,
    *,
    max_rt: float,
    max_im: float,
    max_ppm: float,
) -> pd.DataFrame:
    reference = load_reference_coordinates(study)
    merged = merge_library_reference_coordinates(library, reference)

    blueprint = blueprint_db_path(study)
    min_mz = max(1.0, float(merged["ProductMz"].min()) - 5.0)
    max_mz = float(merged["ProductMz"].max()) + 5.0
    background = build_background_fragments(blueprint, min_product_mz=min_mz, max_product_mz=max_mz)

    with sqlite3.connect(blueprint) as connection:
        windows = shared.load_windows_with_identity(connection)
    merged["WindowGroup"] = [
        oracle.assign_window_group(float(mz), float(im), windows)
        for mz, im in zip(merged["PrecursorMz"], merged["ReferenceIM"])
    ]
    if merged["WindowGroup"].isna().any():
        raise RuntimeError("Could not assign all target transitions to a DIA window group")

    background_by_group = {
        int(group): table.sort_values("ProductMz", kind="mergesort").reset_index(drop=True)
        for group, table in background.groupby("WindowGroup", sort=False)
    }
    rows: list[dict[str, Any]] = []
    for row in merged.itertuples(index=False):
        key = f"{row.PeptideSequence}/{int(row.PrecursorCharge)}"
        group = int(row.WindowGroup)
        collision_count, competitor_precursors, collision_intensity, nearest_ppm = collision_metrics_for_transition(
            float(row.ProductMz),
            key,
            float(row.ReferenceRT),
            float(row.ReferenceIM),
            group,
            background_by_group.get(group, pd.DataFrame(columns=background.columns)),
            max_rt=max_rt,
            max_im=max_im,
            max_ppm=max_ppm,
        )
        rows.append(
            {
                "PeptideSequence": row.PeptideSequence,
                "PrecursorCharge": int(row.PrecursorCharge),
                "TransitionGroupId": row.TransitionGroupId,
                "TransitionOrdinal": int(row.TransitionOrdinal),
                "ProductMz": float(row.ProductMz),
                "LibraryIntensity": float(row.LibraryIntensity),
                "WindowGroup": group,
                "CollisionCount": collision_count,
                "CompetitorPrecursors": competitor_precursors,
                "CollisionPredictedIntensity": collision_intensity,
                "NearestCollisionPPM": nearest_ppm,
            }
        )
    result = pd.DataFrame(rows)
    # Rank independently within each target. Lower geometric collision burden is more specific;
    # predicted target intensity is used only as a deterministic tie-breaker.
    result = result.sort_values(
        [
            "PeptideSequence",
            "PrecursorCharge",
            "CompetitorPrecursors",
            "CollisionCount",
            "CollisionPredictedIntensity",
            "LibraryIntensity",
            "ProductMz",
        ],
        ascending=[True, True, True, True, True, False, True],
        kind="mergesort",
    )
    result["SpecificityRank"] = result.groupby(
        ["PeptideSequence", "PrecursorCharge"], sort=False
    ).cumcount() + 1
    return result.sort_values(
        ["PeptideSequence", "PrecursorCharge", "TransitionOrdinal"], kind="mergesort"
    ).reset_index(drop=True)


def integrate_transition_grid(
    sorted_mz: np.ndarray,
    sorted_intensity: np.ndarray,
    sorted_mobility: np.ndarray,
    *,
    product_mz: float,
    target_im: float,
    ppm_windows: tuple[float, ...],
    im_windows: tuple[float, ...],
) -> np.ndarray:
    result = np.zeros((len(ppm_windows), len(im_windows)), dtype=float)
    max_delta = product_mz * max(ppm_windows) * 1e-6
    left = int(np.searchsorted(sorted_mz, product_mz - max_delta, side="left"))
    right = int(np.searchsorted(sorted_mz, product_mz + max_delta, side="right"))
    if right <= left:
        return result
    mz = sorted_mz[left:right]
    intensity = sorted_intensity[left:right]
    mobility = sorted_mobility[left:right]
    ppm_error = np.abs(mz - product_mz) / product_mz * 1e6
    im_error = np.abs(mobility - target_im)
    for ppm_index, ppm in enumerate(ppm_windows):
        mz_mask = ppm_error <= ppm
        if not np.any(mz_mask):
            continue
        for im_index, im_half_window in enumerate(im_windows):
            mask = mz_mask & (im_error <= im_half_window)
            if np.any(mask):
                result[ppm_index, im_index] = float(np.sum(intensity[mask]))
    return result


def configuration_rows(
    rt_windows: tuple[float, ...],
    im_windows: tuple[float, ...],
    ppm_windows: tuple[float, ...],
    subsets: tuple[int, ...],
) -> list[dict[str, Any]]:
    rows = []
    for rt in rt_windows:
        for im in im_windows:
            for ppm in ppm_windows:
                for subset in subsets:
                    rows.append(
                        {
                            "RTWindow": rt,
                            "IMWindow": im,
                            "FragmentPPM": ppm,
                            "TransitionSubset": subset,
                            "ConfigId": f"rt{rt:g}_im{im:g}_ppm{ppm:g}_top{subset}",
                        }
                    )
    return rows


def condition_effects(
    measurements: pd.DataFrame,
    peptide_truth: pd.DataFrame,
    value_column: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (sequence, charge), group in measurements.groupby(["sequence", "charge"], sort=True):
        control = pd.to_numeric(
            group.loc[group["condition"].eq("control"), value_column], errors="coerce"
        )
        treatment = pd.to_numeric(
            group.loc[group["condition"].eq("treatment"), value_column], errors="coerce"
        )
        control = control[control > 0]
        treatment = treatment[treatment > 0]
        effect = (
            benchmark.safe_log2_ratio(float(treatment.mean()), float(control.mean()))
            if len(control) and len(treatment)
            else math.nan
        )
        rows.append(
            {
                "PeptideSequence": sequence,
                "PrecursorCharge": int(charge),
                "RawSweepLog2FC": effect,
            }
        )
    effects = pd.DataFrame(rows)
    return peptide_truth.merge(
        effects,
        on=["PeptideSequence", "PrecursorCharge"],
        how="left",
        validate="one_to_one",
    )


def summarize_configuration(
    table: pd.DataFrame,
    peptide_truth: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    run_rows = []
    for _, group in table.groupby("run_name", sort=True):
        proxy_metric = benchmark.metric(
            benchmark.positive_log2(group["RealizedEventProxy"]),
            benchmark.positive_log2(group["RawSweepIntensity"]),
        )
        raw_open_metric = benchmark.metric(
            benchmark.positive_log2(group["RawSweepIntensity"]),
            benchmark.positive_log2(group["intensity"]),
        )
        run_rows.append((proxy_metric, raw_open_metric))
    proxy_r = float(np.nanmean([item[0].pearson_r for item in run_rows]))
    proxy_slope = float(np.nanmean([item[0].slope for item in run_rows]))
    raw_open_r = float(np.nanmean([item[1].pearson_r for item in run_rows]))
    raw_open_slope = float(np.nanmean([item[1].slope for item in run_rows]))

    effects = condition_effects(table, peptide_truth, "RawSweepIntensity")
    realized_raw = benchmark.metric(effects["RealizedConditionLog2FC"], effects["RawSweepLog2FC"])
    raw_open_effect = benchmark.metric(effects["RawSweepLog2FC"], effects["ObservedLog2FC"])
    positive_fraction = float(pd.to_numeric(table["RawSweepIntensity"], errors="coerce").gt(0).mean())

    return {
        **config,
        "TargetRunRows": int(len(table)),
        "RawPositiveFraction": positive_fraction,
        "MeanProxyVsRawPearsonR": proxy_r,
        "MeanProxyVsRawSlope": proxy_slope,
        "MeanRawVsOpenDIAPearsonR": raw_open_r,
        "MeanRawVsOpenDIASlope": raw_open_slope,
        "RealizedVsRawEffectN": realized_raw.n,
        "RealizedVsRawEffectPearsonR": realized_raw.pearson_r,
        "RealizedVsRawEffectSlope": realized_raw.slope,
        "RealizedVsRawEffectMAE": realized_raw.mae,
        "RawVsOpenDIAEffectPearsonR": raw_open_effect.pearson_r,
        "RawVsOpenDIAEffectSlope": raw_open_effect.slope,
        "RawVsOpenDIAEffectMAE": raw_open_effect.mae,
    }


def factor_summary(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = [
        "MeanProxyVsRawPearsonR",
        "MeanProxyVsRawSlope",
        "RealizedVsRawEffectPearsonR",
        "RealizedVsRawEffectSlope",
        "MeanRawVsOpenDIASlope",
        "RawVsOpenDIAEffectSlope",
        "RawPositiveFraction",
    ]
    for factor in ("RTWindow", "IMWindow", "FragmentPPM", "TransitionSubset"):
        for value, group in table.groupby(factor, sort=True):
            row = {"Factor": factor, "Level": value, "Configurations": int(len(group))}
            for metric_name in metrics:
                row[metric_name] = float(pd.to_numeric(group[metric_name], errors="coerce").mean())
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--study-results-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rt-windows", default=",".join(str(value) for value in DEFAULT_RT_WINDOWS))
    parser.add_argument("--im-windows", default=",".join(str(value) for value in DEFAULT_IM_WINDOWS))
    parser.add_argument("--fragment-ppm", default=",".join(str(value) for value in DEFAULT_PPM_WINDOWS))
    parser.add_argument("--transition-subsets", default=",".join(str(value) for value in DEFAULT_SUBSETS))
    args = parser.parse_args()

    try:
        rt_windows = parse_float_grid(args.rt_windows, "--rt-windows")
        im_windows = parse_float_grid(args.im_windows, "--im-windows")
        ppm_windows = parse_float_grid(args.fragment_ppm, "--fragment-ppm")
        subsets = parse_int_grid(args.transition_subsets, "--transition-subsets")
    except ValueError as exc:
        parser.error(str(exc))
    if max(subsets) > 8:
        parser.error("transition subsets cannot exceed the eight frozen benchmark transitions")

    study = args.study_dir.resolve()
    study_results = args.study_results_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    manifest_path = study / "OpenSwathTimSim.study_manifest.tsv"
    truth_path = study / "OpenSwathTimSim.realized_truth.tsv"
    library_path = study / "OpenSwathTimSim.target_only.transitions.tsv"
    observed_path = study_results / "precursor_measurements.tsv"
    peptide_path = study_results / "peptide_effects.tsv"
    for path in (manifest_path, truth_path, library_path, observed_path, peptide_path):
        if not path.is_file():
            raise SystemExit(f"Missing sweep input: {path}")

    manifest = pd.read_csv(manifest_path, sep="\t")
    truth = pd.read_csv(truth_path, sep="\t", low_memory=False)
    library = load_library_table(library_path)
    observed = pd.read_csv(observed_path, sep="\t", low_memory=False)
    peptide_truth = pd.read_csv(peptide_path, sep="\t", low_memory=False)
    peptide_truth["PeptideSequence"] = peptide_truth["PeptideSequence"].astype(str).str.upper()
    peptide_truth["PrecursorCharge"] = pd.to_numeric(peptide_truth["PrecursorCharge"], errors="raise").astype(int)
    observed["sequence"] = observed["sequence"].astype(str).str.upper()
    observed["charge"] = pd.to_numeric(observed["charge"], errors="raise").astype(int)
    observed["run_name"] = observed["run_name"].astype(str)

    print("Ranking transition specificity from the complete TimSim blueprint fragment universe...")
    specificity = rank_transition_specificity(
        study,
        library,
        max_rt=max(rt_windows),
        max_im=max(im_windows),
        max_ppm=max(ppm_windows),
    )
    specificity.to_csv(out / "transition_specificity.tsv", sep="\t", index=False)
    rank_map = {
        (str(row.PeptideSequence), int(row.PrecursorCharge), int(row.TransitionOrdinal)): int(row.SpecificityRank)
        for row in specificity.itertuples(index=False)
    }

    transition_map = oracle.load_transition_map(library_path)
    run_names = [str(value) for value in manifest.sort_values("RunOrdinal")["RunName"]]
    config_list = configuration_rows(rt_windows, im_windows, ppm_windows, subsets)
    base_rows: list[dict[str, Any]] = []
    raw_blocks: list[np.ndarray] = []
    config_index = {row["ConfigId"]: index for index, row in enumerate(config_list)}

    for run_index, run_name in enumerate(run_names, start=1):
        d_dir = oracle.d_path_for_run(study, run_name)
        frame_records = oracle.load_dia_frames(d_dir)
        if any(frame.window_group is None for frame in frame_records):
            raise RuntimeError(
                f"DIA frame/window-group mapping unavailable for {run_name}; sweep refuses unmapped frames"
            )
        targets = oracle.load_targets_for_run(study, run_name, truth, transition_map)
        by_group: dict[int, list[int]] = defaultdict(list)
        for target_index, target in enumerate(targets):
            by_group[int(target.window_group)].append(target_index)

        # target x transition x RT x ppm x IM
        accumulator = np.zeros(
            (len(targets), 8, len(rt_windows), len(ppm_windows), len(im_windows)),
            dtype=np.float64,
        )
        dataset, backend = oracle.open_raw_dataset(d_dir)
        print(
            f"[{run_index}/{len(run_names)}] {run_name}: {len(targets)} targets, "
            f"{len(frame_records)} DIA frames, backend={backend}"
        )
        try:
            for frame_number, frame_record in enumerate(frame_records, start=1):
                candidate_indices = by_group.get(int(frame_record.window_group), [])
                active: list[tuple[int, float]] = []
                for target_index in candidate_indices:
                    rt_error = abs(targets[target_index].realized_rt - frame_record.rt_seconds)
                    if rt_error <= max(rt_windows):
                        active.append((target_index, rt_error))
                if not active:
                    continue
                frame = oracle.get_raw_frame(dataset, frame_record.frame_id)
                mz, intensity, mobility = oracle.frame_arrays(dataset, frame, frame_record.frame_id)
                order = np.argsort(mz, kind="mergesort")
                sorted_mz = mz[order]
                sorted_intensity = intensity[order]
                sorted_mobility = mobility[order]
                for target_index, rt_error in active:
                    target = targets[target_index]
                    if len(target.product_mz) != 8:
                        raise RuntimeError(
                            f"Sweep requires eight frozen transitions per target; "
                            f"{target.sequence}/{target.charge} has {len(target.product_mz)}"
                        )
                    rt_mask = np.asarray([rt_error <= window for window in rt_windows], dtype=float)
                    for transition_index, product_mz in enumerate(target.product_mz):
                        grid = integrate_transition_grid(
                            sorted_mz,
                            sorted_intensity,
                            sorted_mobility,
                            product_mz=product_mz,
                            target_im=target.realized_im,
                            ppm_windows=ppm_windows,
                            im_windows=im_windows,
                        )
                        accumulator[target_index, transition_index] += (
                            rt_mask[:, None, None] * grid[None, :, :]
                        )
                if frame_number % 50 == 0:
                    print(f"  processed {frame_number}/{len(frame_records)} DIA frames")
        finally:
            close = getattr(dataset, "close", None)
            if callable(close):
                close()

        observed_run = observed[observed["run_name"].eq(run_name)].copy()
        observed_by_key = {
            (str(row.sequence), int(row.charge)): row
            for row in observed_run.itertuples(index=False)
        }
        run_raw = np.zeros((len(targets), len(config_list)), dtype=np.float64)
        for target_index, target in enumerate(targets):
            key = (target.sequence, target.charge)
            observed_row = observed_by_key.get(key)
            if observed_row is None:
                raise RuntimeError(f"Missing OpenDIA study measurement for {run_name} {key}")
            ranks = np.asarray(
                [rank_map[(target.sequence, target.charge, index + 1)] for index in range(8)],
                dtype=int,
            )
            base_rows.append(
                {
                    "RunName": run_name,
                    "run_name": run_name,
                    "condition": str(observed_row.condition),
                    "PeptideSequence": target.sequence,
                    "sequence": target.sequence,
                    "PrecursorCharge": target.charge,
                    "charge": target.charge,
                    "ProteinId": target.protein_id,
                    "RealizedInputEvents": target.realized_input_events,
                    "RealizedEventProxy": target.realized_event_proxy,
                    "intensity": float(observed_row.intensity),
                }
            )
            for rt_index, rt in enumerate(rt_windows):
                for ppm_index, ppm in enumerate(ppm_windows):
                    for im_index, im in enumerate(im_windows):
                        transition_signal = accumulator[
                            target_index, :, rt_index, ppm_index, im_index
                        ]
                        for subset in subsets:
                            selected = ranks <= subset
                            config_id = f"rt{rt:g}_im{im:g}_ppm{ppm:g}_top{subset}"
                            run_raw[target_index, config_index[config_id]] = float(
                                np.sum(transition_signal[selected])
                            )
        raw_blocks.append(run_raw)

    base = pd.DataFrame(base_rows)
    raw_matrix = np.vstack(raw_blocks)
    if raw_matrix.shape != (len(base), len(config_list)):
        raise RuntimeError(
            f"Sweep matrix shape mismatch: raw={raw_matrix.shape}, "
            f"expected={(len(base), len(config_list))}"
        )

    summaries: list[dict[str, Any]] = []
    measurement_path = out / "sweep_measurements.tsv.gz"
    if measurement_path.exists():
        measurement_path.unlink()
    for config_number, config in enumerate(config_list):
        group = base.copy()
        group["RTWindow"] = config["RTWindow"]
        group["IMWindow"] = config["IMWindow"]
        group["FragmentPPM"] = config["FragmentPPM"]
        group["TransitionSubset"] = config["TransitionSubset"]
        group["RawSweepIntensity"] = raw_matrix[:, config_number]
        summaries.append(summarize_configuration(group, peptide_truth, config))
        group.to_csv(
            measurement_path,
            sep="\t",
            index=False,
            compression="gzip",
            mode="wt" if config_number == 0 else "at",
            header=config_number == 0,
        )
    summary = pd.DataFrame(summaries).sort_values(
        ["TransitionSubset", "RTWindow", "IMWindow", "FragmentPPM"],
        ascending=[False, True, True, True],
        kind="mergesort",
    )
    summary.to_csv(out / "sweep_summary.tsv", sep="\t", index=False)
    factors = factor_summary(summary)
    factors.to_csv(out / "sweep_factor_summary.tsv", sep="\t", index=False)

    baseline_mask = (
        np.isclose(summary["RTWindow"], max(rt_windows))
        & np.isclose(summary["IMWindow"], max(im_windows))
        & np.isclose(summary["FragmentPPM"], max(ppm_windows))
        & summary["TransitionSubset"].eq(max(subsets))
    )
    tight_subset = 4 if 4 in subsets else min(subsets)
    tight_mask = (
        np.isclose(summary["RTWindow"], min(rt_windows))
        & np.isclose(summary["IMWindow"], min(im_windows))
        & np.isclose(summary["FragmentPPM"], min(ppm_windows))
        & summary["TransitionSubset"].eq(tight_subset)
    )
    if not baseline_mask.any() or not tight_mask.any():
        raise RuntimeError("Sweep grid does not contain required baseline and tight high-specificity configurations")
    baseline = summary.loc[baseline_mask].iloc[0]
    tight = summary.loc[tight_mask].iloc[0]

    proxy_improvement = abs(1.0 - float(baseline.MeanProxyVsRawSlope)) - abs(
        1.0 - float(tight.MeanProxyVsRawSlope)
    )
    effect_improvement = abs(1.0 - float(baseline.RealizedVsRawEffectSlope)) - abs(
        1.0 - float(tight.RealizedVsRawEffectSlope)
    )
    if (
        float(tight.MeanProxyVsRawSlope) >= 0.85
        and float(tight.RealizedVsRawEffectSlope) >= 0.85
        and proxy_improvement >= 0.10
        and effect_improvement >= 0.10
    ):
        interpretation = "interference_supported"
    elif (
        float(tight.MeanProxyVsRawSlope) <= 0.75
        and float(tight.RealizedVsRawEffectSlope) <= 0.75
    ):
        interpretation = "persistent_compression"
    else:
        interpretation = "mixed"

    machine = {
        "sweep_version": "truth_centered_interference_sensitivity_v1",
        "runs": len(run_names),
        "targets": int(peptide_truth.shape[0]),
        "configurations": int(len(summary)),
        "grid": {
            "rt_half_windows_seconds": list(rt_windows),
            "im_half_windows_1_over_k0": list(im_windows),
            "fragment_ppm": list(ppm_windows),
            "transition_subsets": list(subsets),
        },
        "transition_specificity": {
            "source": "complete_TimSim_blueprint_fragment_geometry",
            "uses_opendia": False,
            "max_rt_seconds": max(rt_windows),
            "max_im_1_over_k0": max(im_windows),
            "max_fragment_ppm": max(ppm_windows),
        },
        "baseline": baseline.to_dict(),
        "tight_high_specificity": tight.to_dict(),
        "distance_to_unit_slope_improvement": {
            "proxy_vs_raw": proxy_improvement,
            "realized_effect_vs_raw": effect_improvement,
        },
        "interpretation": interpretation,
    }
    (out / "sweep_summary.json").write_text(
        json.dumps(machine, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# Raw-signal oracle interference sensitivity sweep",
        "",
        "This diagnostic scans the existing generated DIA-PASEF raw files once and evaluates a fixed grid of truth-centered extraction windows. Transition subsets are ranked only by collision geometry in the complete TimSim blueprint fragment universe; OpenDIA is not used for transition ranking or parameter selection.",
        "",
        "## Primary comparison",
        "",
        "| Configuration | Proxy→raw r | Proxy→raw slope | Realized effect→raw r | Realized effect→raw slope | Raw→OpenDIA slope | Raw effect→OpenDIA slope | Positive fraction |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        f"| Baseline: ±{baseline.RTWindow:g}s, ±{baseline.IMWindow:g} 1/K0, ±{baseline.FragmentPPM:g} ppm, all {int(baseline.TransitionSubset)} transitions | {baseline.MeanProxyVsRawPearsonR:.4f} | {baseline.MeanProxyVsRawSlope:.4f} | {baseline.RealizedVsRawEffectPearsonR:.4f} | {baseline.RealizedVsRawEffectSlope:.4f} | {baseline.MeanRawVsOpenDIASlope:.4f} | {baseline.RawVsOpenDIAEffectSlope:.4f} | {baseline.RawPositiveFraction:.4f} |",
        f"| Tight/high-specificity: ±{tight.RTWindow:g}s, ±{tight.IMWindow:g} 1/K0, ±{tight.FragmentPPM:g} ppm, top {int(tight.TransitionSubset)} transitions | {tight.MeanProxyVsRawPearsonR:.4f} | {tight.MeanProxyVsRawSlope:.4f} | {tight.RealizedVsRawEffectPearsonR:.4f} | {tight.RealizedVsRawEffectSlope:.4f} | {tight.MeanRawVsOpenDIASlope:.4f} | {tight.RawVsOpenDIAEffectSlope:.4f} | {tight.RawPositiveFraction:.4f} |",
        "",
        f"Distance-to-unit-slope improvement: proxy→raw **{proxy_improvement:.4f}**; realized-effect→raw **{effect_improvement:.4f}**.",
        "",
        f"Diagnostic interpretation: **{interpretation}**.",
        "",
        "`interference_supported` requires both tight slopes >=0.85 and at least 0.10 improvement toward unit slope. `persistent_compression` requires both tight slopes <=0.75. Other outcomes are reported as `mixed`; these labels are diagnostic gates, not benchmark tuning rules.",
        "",
        "## Factor averages",
        "",
        "| Factor | Level | n configs | Proxy→raw slope | Realized effect→raw slope | Raw→OpenDIA slope | Positive fraction |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in factors.itertuples(index=False):
        lines.append(
            f"| {row.Factor} | {row.Level} | {row.Configurations} | {row.MeanProxyVsRawSlope:.4f} | "
            f"{row.RealizedVsRawEffectSlope:.4f} | {row.MeanRawVsOpenDIASlope:.4f} | {row.RawPositiveFraction:.4f} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation guardrails",
            "",
            "- This sweep is diagnostic only. Do not choose canonical extraction windows by maximizing agreement with OpenDIA.",
            "- Transition specificity is computed from TimSim blueprint geometry before OpenDIA and without raw intensity, so it cannot be circular with OpenDIA scoring.",
            "- A systematic movement of proxy→raw and realized-effect→raw slopes toward one as windows tighten / collisions decrease supports DIA cofragmentation as the source of raw-oracle attenuation.",
            "- Persistent compression across tight, low-collision configurations points upstream toward TimSim fragment/raw response and should be resolved before the 25+25 canonical production study.",
            "",
            "## Outputs",
            "",
            "- `transition_specificity.tsv`: simulator-only per-transition collision burden and within-target specificity rank.",
            "- `sweep_summary.tsv`: one row per extraction/subset configuration.",
            "- `sweep_factor_summary.tsv`: factor-level averages across the fixed sweep grid.",
            "- `sweep_measurements.tsv.gz`: target/run raw intensity for every configuration.",
            "- `sweep_summary.json`: machine-readable baseline/tight comparison and diagnostic gate.",
            "",
        ]
    )
    (out / "raw_oracle_sweep_report.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"Raw-oracle sweep complete: {out}")
    print(
        "Baseline: "
        f"proxy->raw slope={baseline.MeanProxyVsRawSlope:.4f}, "
        f"effect slope={baseline.RealizedVsRawEffectSlope:.4f}"
    )
    print(
        "Tight/top-specific: "
        f"proxy->raw slope={tight.MeanProxyVsRawSlope:.4f}, "
        f"effect slope={tight.RealizedVsRawEffectSlope:.4f}"
    )
    print(f"Interpretation: {interpretation}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
