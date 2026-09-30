#!/usr/bin/env python3
"""Extract truth-centered raw DIA-PASEF fragment signal for frozen TimSim targets.

This diagnostic deliberately does not use OpenDIA-selected RT/IM coordinates, scores,
or peak boundaries. It integrates the frozen target transition m/z values around the
simulator-realized RT/IM coordinates directly from each generated Bruker .d run.
"""
from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

TOOLS_DIR = Path(__file__).resolve().parent
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))
import select_high_signal_precursors as shared


@dataclass(frozen=True)
class FrameRecord:
    frame_id: int
    rt_seconds: float
    window_group: int | None


@dataclass(frozen=True)
class TargetRecord:
    run_name: str
    sequence: str
    charge: int
    protein_id: str
    precursor_mz: float
    realized_rt: float
    realized_im: float
    realized_input_events: float
    realized_event_proxy: float
    transition_group_id: str
    product_mz: tuple[float, ...]
    window_group: int | None


@dataclass
class Accumulator:
    intensity: float
    matched_peaks: int
    frames_considered: int
    transition_intensity: np.ndarray
    transition_peaks: np.ndarray


def normalized_name(value: str) -> str:
    return "".join(ch for ch in str(value).lower() if ch.isalnum())


def first_name(options: Iterable[str], available: Iterable[str], label: str) -> str:
    by_normalized = {normalized_name(value): value for value in available}
    for option in options:
        key = normalized_name(option)
        if key in by_normalized:
            return by_normalized[key]
    raise RuntimeError(f"Could not identify {label}; available: {sorted(available)}")


def optional_name(options: Iterable[str], available: Iterable[str]) -> str | None:
    by_normalized = {normalized_name(value): value for value in available}
    for option in options:
        key = normalized_name(option)
        if key in by_normalized:
            return by_normalized[key]
    return None


def sqlite_tables(connection: sqlite3.Connection) -> dict[str, str]:
    return {
        normalized_name(name): str(name)
        for (name,) in connection.execute("SELECT name FROM sqlite_master WHERE type='table'")
    }


def sqlite_columns(connection: sqlite3.Connection, table: str) -> list[str]:
    quoted = '"' + table.replace('"', '""') + '"'
    return [str(row[1]) for row in connection.execute(f"PRAGMA table_info({quoted})")]


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def load_dia_frames(d_dir: Path) -> list[FrameRecord]:
    analysis = d_dir / "analysis.tdf"
    if not analysis.is_file():
        raise FileNotFoundError(f"Missing Bruker analysis.tdf: {analysis}")
    with sqlite3.connect(analysis) as connection:
        tables = sqlite_tables(connection)
        frames_table = tables.get(normalized_name("Frames"))
        if not frames_table:
            raise RuntimeError(f"Frames table not found in {analysis}")
        frame_columns = sqlite_columns(connection, frames_table)
        frame_id_col = first_name(("Id", "Frame", "FrameId", "frame_id"), frame_columns, "frame ID")
        time_col = first_name(("Time", "Rt", "RetentionTime", "retention_time"), frame_columns, "frame time")

        info_table = None
        for candidate in ("DiaFrameMsMsInfo", "dia_frame_ms_ms_info", "DiaFrameMsMsInfos"):
            info_table = tables.get(normalized_name(candidate))
            if info_table:
                break

        if info_table:
            info_columns = sqlite_columns(connection, info_table)
            info_frame_col = first_name(("Frame", "FrameId", "frame_id"), info_columns, "DIA frame ID")
            group_col = first_name(("WindowGroup", "window_group", "Group"), info_columns, "DIA window group")
            query = (
                f"SELECT f.{quote_identifier(frame_id_col)}, f.{quote_identifier(time_col)}, "
                f"d.{quote_identifier(group_col)} "
                f"FROM {quote_identifier(frames_table)} f "
                f"JOIN {quote_identifier(info_table)} d "
                f"ON f.{quote_identifier(frame_id_col)} = d.{quote_identifier(info_frame_col)} "
                f"ORDER BY f.{quote_identifier(frame_id_col)}"
            )
            return [
                FrameRecord(int(frame_id), float(rt), int(group))
                for frame_id, rt, group in connection.execute(query)
            ]

        # Fail-closed fallback: if there is no explicit frame->window-group mapping,
        # retain non-MS1 frames only. The extractor will require --allow-unmapped-dia-frames
        # before using these frames for quantitative integration.
        msms_col = optional_name(("MsMsType", "msms_type", "MsLevel", "ms_level"), frame_columns)
        if not msms_col:
            raise RuntimeError(
                f"No DiaFrameMsMsInfo table and no MS/MS frame-type column found in {analysis}"
            )
        query = (
            f"SELECT {quote_identifier(frame_id_col)}, {quote_identifier(time_col)}, "
            f"{quote_identifier(msms_col)} FROM {quote_identifier(frames_table)} "
            f"ORDER BY {quote_identifier(frame_id_col)}"
        )
        result = []
        for frame_id, rt, msms_type in connection.execute(query):
            if int(msms_type or 0) != 0:
                result.append(FrameRecord(int(frame_id), float(rt), None))
        return result


def load_transition_map(path: Path) -> dict[tuple[str, int], tuple[str, tuple[float, ...]]]:
    table = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"PeptideSequence", "PrecursorCharge", "ProductMz", "TransitionGroupId"}
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(f"Transition library missing columns: {sorted(missing)}")
    if "Decoy" in table.columns:
        table = table[pd.to_numeric(table["Decoy"], errors="coerce").fillna(0).eq(0)]
    if "QuantifyingTransition" in table.columns:
        quant = pd.to_numeric(table["QuantifyingTransition"], errors="coerce").fillna(0).ne(0)
        table = table[quant]
    result: dict[tuple[str, int], tuple[str, tuple[float, ...]]] = {}
    for key, group in table.groupby(["PeptideSequence", "PrecursorCharge"], sort=False):
        sequence = str(key[0]).strip().upper()
        charge = int(key[1])
        mz = tuple(float(value) for value in pd.to_numeric(group["ProductMz"], errors="raise"))
        if not mz:
            continue
        group_ids = set(group["TransitionGroupId"].astype(str))
        if len(group_ids) != 1:
            raise RuntimeError(f"Expected one transition group for {(sequence, charge)}; found {group_ids}")
        result[(sequence, charge)] = (next(iter(group_ids)), mz)
    if not result:
        raise RuntimeError(f"No target quantifying transitions found in {path}")
    return result


def assign_window_group(precursor_mz: float, mobility: float, windows: list[Any]) -> int | None:
    matches = []
    for window in windows:
        if (
            window.mz_lower <= precursor_mz <= window.mz_upper
            and window.im_lower <= mobility <= window.im_upper
        ):
            mz_center = 0.5 * (window.mz_lower + window.mz_upper)
            im_center = 0.5 * (window.im_lower + window.im_upper)
            mz_scale = max(window.mz_upper - window.mz_lower, 1e-12)
            im_scale = max(window.im_upper - window.im_lower, 1e-12)
            distance = abs(precursor_mz - mz_center) / mz_scale + abs(mobility - im_center) / im_scale
            matches.append((distance, int(window.window_group)))
    return min(matches)[1] if matches else None


def load_targets_for_run(
    study_dir: Path,
    run_name: str,
    truth: pd.DataFrame,
    transitions: dict[tuple[str, int], tuple[str, tuple[float, ...]]],
) -> list[TargetRecord]:
    run_truth = truth[truth["RunName"].astype(str).eq(run_name)].copy()
    if run_truth.empty:
        raise RuntimeError(f"No realized truth rows for {run_name}")
    run_db = study_dir / run_name / "synthetic_data.db"
    if not run_db.is_file():
        raise FileNotFoundError(f"Missing TimSim run database: {run_db}")
    with sqlite3.connect(run_db) as connection:
        windows = shared.load_windows_with_identity(connection)
    targets = []
    for row in run_truth.itertuples(index=False):
        sequence = str(row.PeptideSequence).strip().upper()
        charge = int(row.PrecursorCharge)
        key = (sequence, charge)
        transition = transitions.get(key)
        if transition is None:
            raise RuntimeError(f"Frozen target {key} is missing from target-only transition library")
        precursor_mz = float(row.PrecursorMz)
        mobility = float(row.RealizedIMApexSQLite)
        window_group = assign_window_group(precursor_mz, mobility, windows)
        if window_group is None:
            raise RuntimeError(
                f"Could not assign DIA window group for {run_name} {sequence}/{charge} "
                f"at m/z={precursor_mz:.6f}, IM={mobility:.6f}"
            )
        group_id, product_mz = transition
        targets.append(
            TargetRecord(
                run_name=run_name,
                sequence=sequence,
                charge=charge,
                protein_id=str(row.ProteinId),
                precursor_mz=precursor_mz,
                realized_rt=float(row.RealizedRTApex),
                realized_im=mobility,
                realized_input_events=float(row.RealizedInputEvents),
                realized_event_proxy=float(row.RealizedEventProxy),
                transition_group_id=group_id,
                product_mz=product_mz,
                window_group=window_group,
            )
        )
    return targets


def open_raw_dataset(d_dir: Path) -> tuple[Any, str]:
    errors: list[str] = []
    dataset_cls = None
    backend = ""
    try:
        from imspy_core.timstof import TimsDatasetDIA  # type: ignore

        dataset_cls = TimsDatasetDIA
        backend = "imspy_core.timstof.TimsDatasetDIA"
    except Exception as exc:  # pragma: no cover - depends on installed TimSim stack
        errors.append(f"imspy_core: {exc}")
    if dataset_cls is None:
        try:
            from imspy.timstof import TimsDatasetDIA  # type: ignore

            dataset_cls = TimsDatasetDIA
            backend = "imspy.timstof.TimsDatasetDIA"
        except Exception as exc:  # pragma: no cover - compatibility fallback
            errors.append(f"imspy: {exc}")
    if dataset_cls is None:
        raise RuntimeError(
            "Could not import a TimSim-compatible DIA raw reader. "
            "Run ./scripts/setup.sh and verify imspy-core/imspy-connector. "
            + " | ".join(errors)
        )
    for kwargs in ({"in_memory": False}, {}):
        try:
            return dataset_cls(str(d_dir), **kwargs), backend
        except TypeError:
            continue
    return dataset_cls(str(d_dir)), backend


def get_raw_frame(dataset: Any, frame_id: int) -> Any:
    for name in ("get_frame", "get_tims_frame"):
        method = getattr(dataset, name, None)
        if callable(method):
            return method(frame_id)
    try:
        return dataset[frame_id]
    except Exception as exc:
        raise RuntimeError(
            f"Raw reader exposes neither get_frame/get_tims_frame nor usable indexing; "
            f"dataset type={type(dataset)!r}"
        ) from exc


def dataframe_from_frame(frame: Any) -> pd.DataFrame | None:
    for name in ("df", "to_dataframe", "as_dataframe"):
        method = getattr(frame, name, None)
        if callable(method):
            try:
                value = method()
            except TypeError:
                continue
            if isinstance(value, pd.DataFrame):
                return value
            if isinstance(value, dict):
                return pd.DataFrame(value)
    return None


def attribute_array(obj: Any, aliases: Iterable[str]) -> np.ndarray | None:
    attrs = {normalized_name(name): name for name in dir(obj) if not name.startswith("__")}
    for alias in aliases:
        name = attrs.get(normalized_name(alias))
        if not name:
            continue
        value = getattr(obj, name)
        if callable(value):
            continue
        try:
            array = np.asarray(value)
        except Exception:
            continue
        if array.ndim == 1:
            return array
    return None


def resolve_scan_mobility(dataset: Any, frame: Any, frame_id: int, scans: np.ndarray) -> np.ndarray | None:
    for owner in (frame, dataset):
        for name in (
            "scan_to_inverse_mobility",
            "scan_to_inv_mobility",
            "scan_to_one_over_k0",
            "scan_num_to_one_over_k0",
        ):
            method = getattr(owner, name, None)
            if not callable(method):
                continue
            for args in ((scans,), (frame_id, scans)):
                try:
                    value = np.asarray(method(*args), dtype=float)
                except Exception:
                    continue
                if value.shape == scans.shape:
                    return value
    return None


def frame_arrays(dataset: Any, frame: Any, frame_id: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    table = dataframe_from_frame(frame)
    if table is not None:
        cols = list(table.columns)
        mz_col = first_name(("mz", "m/z", "mz_values", "mzs"), cols, "frame m/z")
        intensity_col = first_name(("intensity", "intensities", "intensity_values"), cols, "frame intensity")
        mobility_col = optional_name(
            ("mobility", "inv_mobility", "inverse_mobility", "inv_ion_mobility", "one_over_k0", "mobility_values"),
            cols,
        )
        mz = pd.to_numeric(table[mz_col], errors="coerce").to_numpy(dtype=float)
        intensity = pd.to_numeric(table[intensity_col], errors="coerce").to_numpy(dtype=float)
        mobility = (
            pd.to_numeric(table[mobility_col], errors="coerce").to_numpy(dtype=float)
            if mobility_col
            else None
        )
        if mobility is None:
            scan_col = optional_name(("scan", "scan_id", "scan_index"), cols)
            if scan_col:
                scans = pd.to_numeric(table[scan_col], errors="coerce").to_numpy(dtype=float)
                mobility = resolve_scan_mobility(dataset, frame, frame_id, scans)
    else:
        mz = attribute_array(frame, ("mz", "mz_values", "mzs"))
        intensity = attribute_array(frame, ("intensity", "intensities", "intensity_values"))
        mobility = attribute_array(
            frame,
            ("mobility", "inv_mobility", "inverse_mobility", "inv_ion_mobility", "one_over_k0", "mobility_values"),
        )
        if mobility is None:
            scans = attribute_array(frame, ("scan", "scans", "scan_ids", "scan_indices"))
            if scans is not None:
                mobility = resolve_scan_mobility(dataset, frame, frame_id, scans)
    if mz is None or intensity is None or mobility is None:
        available = sorted(name for name in dir(frame) if not name.startswith("__"))
        raise RuntimeError(
            "Could not normalize raw frame to m/z, intensity, and 1/K0 arrays. "
            f"frame type={type(frame)!r}; available attributes={available[:80]}"
        )
    if not (len(mz) == len(intensity) == len(mobility)):
        raise RuntimeError(
            f"Raw frame array length mismatch for frame {frame_id}: "
            f"mz={len(mz)}, intensity={len(intensity)}, mobility={len(mobility)}"
        )
    mask = np.isfinite(mz) & np.isfinite(intensity) & np.isfinite(mobility) & (intensity >= 0)
    return mz[mask], intensity[mask], mobility[mask]


def integrate_transition(
    sorted_mz: np.ndarray,
    sorted_intensity: np.ndarray,
    sorted_mobility: np.ndarray,
    product_mz: float,
    target_im: float,
    ppm: float,
    im_half_window: float,
) -> tuple[float, int]:
    delta = product_mz * ppm * 1e-6
    left = int(np.searchsorted(sorted_mz, product_mz - delta, side="left"))
    right = int(np.searchsorted(sorted_mz, product_mz + delta, side="right"))
    if right <= left:
        return 0.0, 0
    mobility = sorted_mobility[left:right]
    mask = np.abs(mobility - target_im) <= im_half_window
    if not np.any(mask):
        return 0.0, 0
    intensity = sorted_intensity[left:right][mask]
    return float(np.sum(intensity)), int(np.count_nonzero(mask))


def d_path_for_run(study_dir: Path, run_name: str) -> Path:
    candidates = [study_dir / "tdfs" / f"{run_name}.d", study_dir / f"{run_name}.d"]
    for path in candidates:
        if path.is_dir():
            return path
    raise FileNotFoundError(f"Could not locate generated .d for {run_name}; checked {candidates}")


def probe_reader(study_dir: Path, run_name: str) -> None:
    d_dir = d_path_for_run(study_dir, run_name)
    frames = load_dia_frames(d_dir)
    if not frames:
        raise RuntimeError(f"No DIA frames found in {d_dir}")
    dataset, backend = open_raw_dataset(d_dir)
    try:
        frame = get_raw_frame(dataset, frames[0].frame_id)
        mz, intensity, mobility = frame_arrays(dataset, frame, frames[0].frame_id)
    finally:
        close = getattr(dataset, "close", None)
        if callable(close):
            close()
    print(f"Raw reader backend: {backend}")
    print(f"Run: {run_name}")
    print(
        f"Frame {frames[0].frame_id}: peaks={len(mz)}, RT={frames[0].rt_seconds:.4f}s, "
        f"window_group={frames[0].window_group}"
    )
    if len(mz):
        print(
            f"m/z range={float(np.min(mz)):.4f}..{float(np.max(mz)):.4f}; "
            f"IM range={float(np.min(mobility)):.6f}..{float(np.max(mobility)):.6f}; "
            f"intensity sum={float(np.sum(intensity)):.1f}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rt-half-window", type=float, default=6.0)
    parser.add_argument("--im-half-window", type=float, default=0.03)
    parser.add_argument("--fragment-ppm", type=float, default=25.0)
    parser.add_argument("--probe-reader", action="store_true")
    parser.add_argument(
        "--allow-unmapped-dia-frames",
        action="store_true",
        help="Allow DIA frames without frame-to-window-group metadata (not recommended for canonical oracle runs).",
    )
    args = parser.parse_args()
    if args.rt_half_window <= 0 or args.im_half_window <= 0 or args.fragment_ppm <= 0:
        parser.error("raw-oracle extraction windows must be positive")

    study = args.study_dir.resolve()
    out = args.out_dir.resolve()
    manifest_path = study / "OpenSwathTimSim.study_manifest.tsv"
    truth_path = study / "OpenSwathTimSim.realized_truth.tsv"
    library_path = study / "OpenSwathTimSim.target_only.transitions.tsv"
    for path in (manifest_path, truth_path, library_path):
        if not path.is_file():
            raise SystemExit(f"Missing study input: {path}")

    manifest = pd.read_csv(manifest_path, sep="\t")
    truth = pd.read_csv(truth_path, sep="\t", low_memory=False)
    transitions = load_transition_map(library_path)
    run_names = [str(value) for value in manifest.sort_values("RunOrdinal")["RunName"]]
    if not run_names:
        raise SystemExit("Study manifest contains no runs")
    if args.probe_reader:
        probe_reader(study, run_names[0])
        return 0

    out.mkdir(parents=True, exist_ok=True)
    measurement_rows: list[dict[str, Any]] = []
    transition_rows: list[dict[str, Any]] = []
    backend_by_run: dict[str, str] = {}

    for run_index, run_name in enumerate(run_names, start=1):
        d_dir = d_path_for_run(study, run_name)
        frame_records = load_dia_frames(d_dir)
        if not frame_records:
            raise RuntimeError(f"No DIA frames found in {d_dir}")
        if any(frame.window_group is None for frame in frame_records) and not args.allow_unmapped_dia_frames:
            raise RuntimeError(
                f"DIA frame/window-group mapping is unavailable for {run_name}. "
                "Refusing canonical extraction; inspect analysis.tdf or pass "
                "--allow-unmapped-dia-frames only for a diagnostic fallback."
            )
        targets = load_targets_for_run(study, run_name, truth, transitions)
        by_group: dict[int | None, list[int]] = defaultdict(list)
        for index, target in enumerate(targets):
            by_group[target.window_group].append(index)
        accumulators = [
            Accumulator(
                intensity=0.0,
                matched_peaks=0,
                frames_considered=0,
                transition_intensity=np.zeros(len(target.product_mz), dtype=float),
                transition_peaks=np.zeros(len(target.product_mz), dtype=np.int64),
            )
            for target in targets
        ]
        dataset, backend = open_raw_dataset(d_dir)
        backend_by_run[run_name] = backend
        print(
            f"[{run_index}/{len(run_names)}] {run_name}: {len(targets)} targets, "
            f"{len(frame_records)} DIA frames, backend={backend}"
        )
        try:
            for frame_number, frame_record in enumerate(frame_records, start=1):
                if frame_record.window_group is None:
                    candidate_indices = range(len(targets))
                else:
                    candidate_indices = by_group.get(frame_record.window_group, [])
                active = [
                    index
                    for index in candidate_indices
                    if abs(targets[index].realized_rt - frame_record.rt_seconds) <= args.rt_half_window
                ]
                if not active:
                    continue
                frame = get_raw_frame(dataset, frame_record.frame_id)
                mz, intensity, mobility = frame_arrays(dataset, frame, frame_record.frame_id)
                order = np.argsort(mz, kind="mergesort")
                sorted_mz = mz[order]
                sorted_intensity = intensity[order]
                sorted_mobility = mobility[order]
                for index in active:
                    target = targets[index]
                    accumulator = accumulators[index]
                    accumulator.frames_considered += 1
                    for transition_index, product_mz in enumerate(target.product_mz):
                        signal, peaks = integrate_transition(
                            sorted_mz,
                            sorted_intensity,
                            sorted_mobility,
                            product_mz,
                            target.realized_im,
                            args.fragment_ppm,
                            args.im_half_window,
                        )
                        accumulator.intensity += signal
                        accumulator.matched_peaks += peaks
                        accumulator.transition_intensity[transition_index] += signal
                        accumulator.transition_peaks[transition_index] += peaks
                if frame_number % 50 == 0:
                    print(f"  processed {frame_number}/{len(frame_records)} DIA frames")
        finally:
            close = getattr(dataset, "close", None)
            if callable(close):
                close()

        for target, accumulator in zip(targets, accumulators):
            matched_transitions = int(np.count_nonzero(accumulator.transition_intensity > 0))
            measurement_rows.append(
                {
                    "RunName": run_name,
                    "PeptideSequence": target.sequence,
                    "PrecursorCharge": target.charge,
                    "ProteinId": target.protein_id,
                    "TransitionGroupId": target.transition_group_id,
                    "PrecursorMz": target.precursor_mz,
                    "RealizedRTApex": target.realized_rt,
                    "RealizedIMApex": target.realized_im,
                    "RealizedInputEvents": target.realized_input_events,
                    "RealizedEventProxy": target.realized_event_proxy,
                    "RawOracleIntensity": accumulator.intensity,
                    "MatchedTransitions": matched_transitions,
                    "TotalTransitions": len(target.product_mz),
                    "MatchedRawPeaks": accumulator.matched_peaks,
                    "FramesConsidered": accumulator.frames_considered,
                    "WindowGroup": target.window_group,
                }
            )
            for transition_index, product_mz in enumerate(target.product_mz):
                transition_rows.append(
                    {
                        "RunName": run_name,
                        "PeptideSequence": target.sequence,
                        "PrecursorCharge": target.charge,
                        "TransitionGroupId": target.transition_group_id,
                        "TransitionIndex": transition_index + 1,
                        "ProductMz": product_mz,
                        "RawOracleTransitionIntensity": float(accumulator.transition_intensity[transition_index]),
                        "MatchedRawPeaks": int(accumulator.transition_peaks[transition_index]),
                    }
                )

    measurements = pd.DataFrame(measurement_rows)
    transitions_out = pd.DataFrame(transition_rows)
    measurements.to_csv(out / "raw_signal_oracle_measurements.tsv", sep="\t", index=False)
    transitions_out.to_csv(out / "raw_signal_oracle_transitions.tsv", sep="\t", index=False)
    metadata = {
        "oracle_version": "truth_centered_raw_ms2_v1",
        "study_dir": str(study),
        "runs": len(run_names),
        "targets_per_run": int(len(measurements) / max(len(run_names), 1)),
        "rt_half_window_seconds": args.rt_half_window,
        "im_half_window_1_over_k0": args.im_half_window,
        "fragment_mz_tolerance_ppm": args.fragment_ppm,
        "reader_backends": backend_by_run,
        "uses_opendia_coordinates": False,
        "uses_opendia_scores": False,
    }
    (out / "raw_signal_oracle_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    positive = pd.to_numeric(measurements["RawOracleIntensity"], errors="coerce").gt(0)
    complete = pd.to_numeric(measurements["MatchedTransitions"], errors="coerce").eq(
        pd.to_numeric(measurements["TotalTransitions"], errors="coerce")
    )
    print(f"Raw oracle measurements: {out / 'raw_signal_oracle_measurements.tsv'}")
    print(f"Positive raw signal: {int(positive.sum())}/{len(measurements)} ({float(positive.mean()):.2%})")
    print(f"All transitions matched: {int(complete.sum())}/{len(measurements)} ({float(complete.mean()):.2%})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
