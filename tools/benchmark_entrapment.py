#!/usr/bin/env python3
"""Benchmark OpenDIA identity-error calibration with external entrapment truth.

The benchmark separates three distinct identification contexts instead of
conflating them:

1. ``run_level``: best precursor/run peak group evaluated with
   ``score_ms2_qvalue`` (or a user-selected run-level q column).
2. ``global_peptide``: one best precursor identity across the experiment,
   evaluated with ``score_peptide_global_qvalue``.
3. ``actual_export``: every target-labelled row that OpenDIA actually writes to
   ``OpenDIA.results.tsv``. This is the most literal empirical FDP of the user-
   visible result set; it is reported both pooled and per run.

True targets are high-signal TimSim precursors. Entrapments are target-labelled
assays whose peptide sequences are known absent from the TimSim runs.
OpenDIA-generated decoys remain a third, independent class used only as the
internal target-decoy null model.

For a single synthetic realization the directly observed quantity is the false
discovery proportion (FDP), not the population expectation FDR.
"""
from __future__ import annotations

import argparse
import io
import json
import math
import re
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DEFAULT_Q_GRID = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10)


class BenchmarkError(RuntimeError):
    pass


@dataclass(frozen=True)
class CurveMetrics:
    roc_auc: float
    average_precision: float
    n_positive: int
    n_negative: int


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise BenchmarkError(f"Required TSV does not exist: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def _normalise_sequence(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.upper()


def _first_existing(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    available = set(columns)
    return next((candidate for candidate in candidates if candidate in available), None)


def _normalise_decoy(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    text = series.astype(str).str.strip().str.lower()
    return numeric.fillna(0).ne(0) | text.isin({"true", "yes", "decoy"})


def _find_workflow(opendia_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.is_file():
            raise BenchmarkError(f"workflow.oswpq does not exist: {path}")
        return path
    candidates = [
        opendia_dir / "intermediates" / "workflow.oswpq",
        opendia_dir / "OpenDIA_intermediates" / "workflow.oswpq",
        opendia_dir / "workflow.oswpq",
    ]
    for path in candidates:
        if path.is_file():
            return path

    discovered = sorted(opendia_dir.rglob("workflow.oswpq"))
    if len(discovered) == 1:
        return discovered[0]
    if len(discovered) > 1:
        raise BenchmarkError(
            "Multiple workflow.oswpq archives were found under the OpenDIA output. "
            "Pass --workflow-oswpq explicitly: " + ", ".join(str(path) for path in discovered)
        )

    raise BenchmarkError(
        "Could not locate workflow.oswpq. Re-run OpenDIA with "
        "workflow:keep_intermediate_files=true (optionally together with "
        "workflow:intermediate_dir) or pass --workflow-oswpq explicitly. Tried: "
        + ", ".join(str(path) for path in candidates)
    )


def _find_results(opendia_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.is_file():
            raise BenchmarkError(f"OpenDIA results TSV does not exist: {path}")
        return path
    path = opendia_dir / "OpenDIA.results.tsv"
    if not path.is_file():
        raise BenchmarkError(f"OpenDIA.results.tsv does not exist: {path}")
    return path


def _read_parquet_member(archive: zipfile.ZipFile, member: str) -> pd.DataFrame:
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise BenchmarkError(
            "pyarrow is required to read workflow.oswpq. Re-run ./scripts/setup.sh."
        ) from exc
    if member not in archive.namelist():
        raise BenchmarkError(f"workflow.oswpq is missing {member}")
    return pq.read_table(io.BytesIO(archive.read(member))).to_pandas()


def _load_truth(
    build_dir: Path,
    explicit_entrapment_truth: Path | None,
) -> tuple[set[tuple[str, int]], set[tuple[str, int]], str, Path]:
    selected = _read_tsv(build_dir / "OpenSwathTimSim.high_signal_selection.tsv")
    entrapment_path = (
        explicit_entrapment_truth.expanduser().resolve()
        if explicit_entrapment_truth is not None
        else build_dir / "OpenSwathTimSim.entrapment_truth.tsv"
    )
    entrapment = _read_tsv(entrapment_path)

    if not {"sequence", "charge"}.issubset(selected.columns):
        raise BenchmarkError("High-signal selection TSV must contain sequence and charge")
    if not {"EntrapmentSequence", "PrecursorCharge"}.issubset(entrapment.columns):
        raise BenchmarkError("Entrapment truth TSV is missing EntrapmentSequence/PrecursorCharge")

    true_keys = set(
        zip(
            _normalise_sequence(selected["sequence"]),
            pd.to_numeric(selected["charge"], errors="raise").astype(int),
        )
    )
    entrapment_keys = set(
        zip(
            _normalise_sequence(entrapment["EntrapmentSequence"]),
            pd.to_numeric(entrapment["PrecursorCharge"], errors="raise").astype(int),
        )
    )
    overlap = true_keys & entrapment_keys
    if overlap:
        raise BenchmarkError(f"True-target and entrapment keys overlap: {list(overlap)[:5]}")

    if "EntrapmentMode" not in entrapment.columns:
        raise BenchmarkError("Entrapment truth TSV is missing EntrapmentMode")
    modes = sorted(entrapment["EntrapmentMode"].dropna().astype(str).unique())
    if not modes:
        raise BenchmarkError("Entrapment truth TSV has no EntrapmentMode value")
    entrapment_mode = ",".join(modes)
    return true_keys, entrapment_keys, entrapment_mode, entrapment_path


def _classify_precursors(
    precursors: pd.DataFrame,
    true_keys: set[tuple[str, int]],
    entrapment_keys: set[tuple[str, int]],
) -> pd.DataFrame:
    required = {"precursor_id", "charge", "decoy", "unmodified_sequence"}
    missing = required - set(precursors.columns)
    if missing:
        raise BenchmarkError(f"Prepared precursor table is missing columns: {sorted(missing)}")

    out = precursors.copy()
    out["sequence"] = _normalise_sequence(out["unmodified_sequence"])
    out["precursor_charge"] = pd.to_numeric(out["charge"], errors="raise").astype(int)
    out["is_decoy"] = _normalise_decoy(out["decoy"])
    keys = list(zip(out["sequence"], out["precursor_charge"]))
    classes = []
    for is_decoy, key in zip(out["is_decoy"], keys):
        if is_decoy:
            classes.append("decoy")
        elif key in true_keys:
            classes.append("true_target")
        elif key in entrapment_keys:
            classes.append("entrapment")
        else:
            classes.append("unknown_target")
    out["truth_class"] = classes

    unknown = out[(~out["is_decoy"]) & (out["truth_class"] == "unknown_target")]
    if not unknown.empty:
        examples = unknown[["sequence", "precursor_charge"]].head(10).to_dict("records")
        raise BenchmarkError(
            f"Prepared library contains {len(unknown)} target precursors absent from true/entrapment truth; examples={examples}"
        )
    return out


def _load_workflow(
    workflow: Path,
    true_keys: set[tuple[str, int]],
    entrapment_keys: set[tuple[str, int]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    with zipfile.ZipFile(workflow) as archive:
        precursors = _classify_precursors(
            _read_parquet_member(archive, "library/precursors.parquet"),
            true_keys,
            entrapment_keys,
        )
        runs = _read_parquet_member(archive, "runs/runs.parquet")
        feature_members = [
            name
            for name in archive.namelist()
            if re.fullmatch(r"runs/run_id=[^/]+/features\.parquet", name)
        ]
        if not feature_members:
            raise BenchmarkError("workflow.oswpq contains no run feature tables")
        frames: list[pd.DataFrame] = []
        for member in feature_members:
            frame = _read_parquet_member(archive, member)
            if "run_id" not in frame.columns:
                archive_run_id = member.split("run_id=", 1)[1].split("/", 1)[0]
                frame["run_id"] = archive_run_id
            frames.append(frame)
        features = pd.concat(frames, ignore_index=True)
    return precursors, runs, features


def _run_role(filename: Any) -> str:
    name = Path(str(filename)).stem
    if "01_baseline" in name:
        return "baseline"
    if "02_biological_replicate" in name:
        return "control_replicate"
    if "03_treatment" in name:
        return "treatment"
    return name


def _run_table(runs: pd.DataFrame) -> pd.DataFrame:
    if "run_id" not in runs.columns or "filename" not in runs.columns:
        raise BenchmarkError(f"runs/runs.parquet must contain run_id/filename; columns={list(runs.columns)}")
    out = runs[["run_id", "filename"]].copy()
    out["run_id_key"] = out["run_id"].astype(str)
    out["run_name"] = out["filename"].map(lambda value: Path(str(value)).stem)
    out["run_role"] = out["filename"].map(_run_role)
    return out


def _numeric_or_nan(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _prepare_scoring_tables(
    precursors: pd.DataFrame,
    runs: pd.DataFrame,
    features: pd.DataFrame,
    *,
    run_score_column: str,
    run_q_column: str,
    run_pep_column: str,
    global_score_column: str,
    global_q_column: str,
    global_pep_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    for column in ("precursor_id", "run_id"):
        if column not in features.columns:
            raise BenchmarkError(f"Feature table lacks {column}; columns={list(features.columns)}")
    if run_score_column not in features.columns:
        raise BenchmarkError(f"Feature run-level score column is absent: {run_score_column}")
    if run_q_column not in features.columns:
        raise BenchmarkError(f"Feature run-level q-value column is absent: {run_q_column}")
    if global_q_column not in features.columns:
        raise BenchmarkError(
            f"Feature global peptide q-value column is absent: {global_q_column}. "
            "Use --global-q-column if the workflow uses a different name."
        )

    run_info = _run_table(runs)
    feature = features.copy()
    feature["run_id_key"] = feature["run_id"].astype(str)
    feature["run_score"] = _numeric_or_nan(feature, run_score_column)
    feature["run_q"] = _numeric_or_nan(feature, run_q_column)
    feature["run_pep"] = _numeric_or_nan(feature, run_pep_column)
    feature["global_score"] = _numeric_or_nan(feature, global_score_column)
    feature["global_q"] = _numeric_or_nan(feature, global_q_column)
    feature["global_pep"] = _numeric_or_nan(feature, global_pep_column)
    feature["candidate_rows"] = 1
    feature["scored_candidate_rows"] = feature["run_score"].notna().astype(int)

    lib = precursors[
        ["precursor_id", "sequence", "precursor_charge", "truth_class", "is_decoy"]
    ].copy()

    # --- run-level universe: one best peak group per precursor/run ---
    run_counts = feature.groupby(["run_id_key", "precursor_id"], as_index=False).agg(
        candidate_rows=("candidate_rows", "sum"),
        scored_candidate_rows=("scored_candidate_rows", "sum"),
    )
    run_sorted = feature.sort_values(
        ["run_id_key", "precursor_id", "run_score"],
        ascending=[True, True, False],
        na_position="last",
    )
    run_best = run_sorted.drop_duplicates(["run_id_key", "precursor_id"], keep="first")
    run_keep = [
        "run_id_key", "precursor_id", "run_score", "run_q", "run_pep",
        "global_score", "global_q", "global_pep",
    ]
    for optional in ("feature_id", "exp_rt", "exp_im", "ms2_area_intensity"):
        if optional in run_best.columns:
            run_keep.append(optional)
    run_best = run_best[run_keep].merge(
        run_counts, on=["run_id_key", "precursor_id"], how="left"
    )
    run_universe = run_info[["run_id_key", "run_name", "run_role"]].assign(_join=1).merge(
        lib.assign(_join=1), on="_join", how="inner"
    ).drop(columns="_join")
    run_universe = run_universe.merge(
        run_best,
        on=["run_id_key", "precursor_id"],
        how="left",
        validate="one_to_one",
    )
    run_universe["candidate_rows"] = run_universe["candidate_rows"].fillna(0).astype(int)
    run_universe["scored_candidate_rows"] = run_universe["scored_candidate_rows"].fillna(0).astype(int)
    # Common aliases used by null/ROC plotting.
    run_universe["score"] = run_universe["run_score"]
    run_universe["reported_q"] = run_universe["run_q"]
    run_universe["reported_pep"] = run_universe["run_pep"]

    # --- global peptide universe: one identity per prepared precursor ---
    global_counts = feature.groupby("precursor_id", as_index=False).agg(
        candidate_rows=("candidate_rows", "sum"),
        scored_candidate_rows=("scored_candidate_rows", "sum"),
    )
    # Global score is the natural ranking when available. If a workflow only
    # synchronized q/PEP but not a global discriminant score, prefer lowest q
    # and then the best run-level score to choose a representative row.
    global_sort = feature.copy()
    if global_sort["global_score"].notna().any():
        global_sort = global_sort.sort_values(
            ["precursor_id", "global_score", "global_q", "run_score"],
            ascending=[True, False, True, False],
            na_position="last",
        )
    else:
        global_sort = global_sort.sort_values(
            ["precursor_id", "global_q", "run_score"],
            ascending=[True, True, False],
            na_position="last",
        )
    global_best = global_sort.drop_duplicates("precursor_id", keep="first")
    global_keep = [
        "precursor_id", "run_id_key", "global_score", "global_q", "global_pep",
        "run_score", "run_q", "run_pep",
    ]
    global_best = global_best[global_keep].merge(global_counts, on="precursor_id", how="left")
    global_best = global_best.merge(
        run_info[["run_id_key", "run_name", "run_role"]],
        on="run_id_key", how="left"
    )
    global_universe = lib.merge(global_best, on="precursor_id", how="left", validate="one_to_one")
    global_universe["candidate_rows"] = global_universe["candidate_rows"].fillna(0).astype(int)
    global_universe["scored_candidate_rows"] = global_universe["scored_candidate_rows"].fillna(0).astype(int)
    global_universe["score"] = global_universe["global_score"].where(
        global_universe["global_score"].notna(), global_universe["run_score"]
    )
    global_universe["reported_q"] = global_universe["global_q"]
    global_universe["reported_pep"] = global_universe["global_pep"]

    return run_universe, global_universe


def _score_sentinel(values: pd.Series) -> float:
    finite = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if finite.size == 0:
        return -1.0
    span = float(np.max(finite) - np.min(finite))
    return float(np.min(finite) - max(1.0, 0.05 * span))


def _roc_pr_points(
    frame: pd.DataFrame,
    *,
    score_column: str = "score",
) -> tuple[pd.DataFrame, pd.DataFrame, CurveMetrics]:
    work = frame[frame["truth_class"].isin(["true_target", "entrapment"])].copy()
    y = (work["truth_class"] == "true_target").astype(int).to_numpy()
    sentinel = _score_sentinel(work[score_column])
    scores = pd.to_numeric(work[score_column], errors="coerce").fillna(sentinel).to_numpy(dtype=float)
    positives = int(y.sum())
    negatives = int(len(y) - positives)
    if positives == 0 or negatives == 0:
        raise BenchmarkError("ROC/PR requires both true targets and entrapment negatives")

    order = np.argsort(-scores, kind="stable")
    y_sorted = y[order]
    score_sorted = scores[order]
    tp = np.cumsum(y_sorted == 1)
    fp = np.cumsum(y_sorted == 0)
    distinct = np.r_[np.where(np.diff(score_sorted) != 0)[0], len(score_sorted) - 1]
    tp_d = tp[distinct].astype(float)
    fp_d = fp[distinct].astype(float)
    thresholds = score_sorted[distinct]

    tpr = tp_d / positives
    fpr = fp_d / negatives
    roc = pd.DataFrame(
        {
            "threshold": np.r_[np.inf, thresholds],
            "fpr": np.r_[0.0, fpr],
            "tpr": np.r_[0.0, tpr],
        }
    )
    # Use NumPy 2.x-compatible trapezoidal integration.
    roc_auc = float(np.trapezoid(roc["tpr"].to_numpy(), roc["fpr"].to_numpy()))

    precision = tp_d / np.maximum(tp_d + fp_d, 1.0)
    recall = tp_d / positives
    pr = pd.DataFrame(
        {
            "threshold": np.r_[np.inf, thresholds],
            "recall": np.r_[0.0, recall],
            "precision": np.r_[1.0, precision],
        }
    )
    tp_increment = np.diff(np.r_[0.0, tp_d])
    average_precision = float(np.sum(precision * tp_increment) / positives)
    return roc, pr, CurveMetrics(roc_auc, average_precision, positives, negatives)


def _confusion_at_q(
    frame: pd.DataFrame,
    q_threshold: float,
    *,
    q_column: str = "reported_q",
) -> dict[str, Any]:
    work = frame[frame["truth_class"].isin(["true_target", "entrapment"])].copy()
    accepted = pd.to_numeric(work[q_column], errors="coerce").le(q_threshold)
    true = work["truth_class"] == "true_target"
    ent = work["truth_class"] == "entrapment"
    tp = int((accepted & true).sum())
    fp = int((accepted & ent).sum())
    fn = int((~accepted & true).sum())
    tn = int((~accepted & ent).sum())
    precision = tp / (tp + fp) if tp + fp else math.nan
    recall = tp / (tp + fn) if tp + fn else math.nan
    specificity = tn / (tn + fp) if tn + fp else math.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else math.nan
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = (tp * tn - fp * fn) / denominator if denominator else math.nan
    # FDP is V / max(R, 1): when there are no discoveries (R=0), FDP is 0.
    empirical_fdp = fp / (tp + fp) if tp + fp else 0.0
    return {
        "q_threshold": q_threshold,
        "tp": tp,
        "fp_entrapment": fp,
        "fn": fn,
        "tn": tn,
        "accepted": tp + fp,
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
        "mcc": mcc,
        "empirical_fdp": empirical_fdp,
    }


def _fdr_calibration(
    frame: pd.DataFrame,
    q_grid: Iterable[float],
    *,
    level: str,
    include_run_scopes: bool,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    groups: list[tuple[str, pd.DataFrame]] = [("pooled", frame)]
    if include_run_scopes:
        run_roles = sorted(str(value) for value in frame["run_role"].dropna().unique())
        groups.extend((role, frame[frame["run_role"] == role]) for role in run_roles)
    for scope, subset in groups:
        for threshold in q_grid:
            row = _confusion_at_q(subset, threshold)
            row["level"] = level
            row["scope"] = scope
            row["fdp_to_nominal_ratio"] = (
                row["empirical_fdp"] / threshold
                if threshold > 0 and np.isfinite(row["empirical_fdp"])
                else math.nan
            )
            rows.append(row)
    return pd.DataFrame(rows)


def _pep_calibration(
    frame: pd.DataFrame,
    *,
    pep_column: str = "reported_pep",
) -> tuple[pd.DataFrame, dict[str, float]]:
    work = frame[frame["truth_class"].isin(["true_target", "entrapment"])].copy()
    work["pep"] = pd.to_numeric(work[pep_column], errors="coerce")
    work = work[work["pep"].notna() & (work["pep"] >= 0) & (work["pep"] <= 1)].copy()
    if work.empty:
        return pd.DataFrame(), {
            "brier_score": math.nan,
            "expected_calibration_error": math.nan,
            "n": 0,
        }
    work["observed_error"] = (work["truth_class"] == "entrapment").astype(float)
    edges = np.array([0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 1.0000001])
    work["pep_bin"] = pd.cut(work["pep"], bins=edges, include_lowest=True, right=False)
    grouped = work.groupby("pep_bin", observed=True).agg(
        n=("observed_error", "size"),
        mean_reported_pep=("pep", "mean"),
        observed_entrapment_error_rate=("observed_error", "mean"),
    ).reset_index()
    grouped["pep_bin"] = grouped["pep_bin"].astype(str)
    brier = float(np.mean((work["pep"].to_numpy() - work["observed_error"].to_numpy()) ** 2))
    ece = float(
        np.sum(
            grouped["n"].to_numpy()
            * np.abs(
                grouped["mean_reported_pep"].to_numpy()
                - grouped["observed_entrapment_error_rate"].to_numpy()
            )
        )
        / len(work)
    )
    return grouped, {"brier_score": brier, "expected_calibration_error": ece, "n": int(len(work))}


def _null_distribution_metrics(frame: pd.DataFrame) -> dict[str, float]:
    ent = pd.to_numeric(
        frame.loc[frame["truth_class"] == "entrapment", "score"], errors="coerce"
    ).dropna().to_numpy(dtype=float)
    decoy = pd.to_numeric(
        frame.loc[frame["truth_class"] == "decoy", "score"], errors="coerce"
    ).dropna().to_numpy(dtype=float)
    if len(ent) == 0 or len(decoy) == 0:
        return {
            "n_entrapment_scored": int(len(ent)),
            "n_decoy_scored": int(len(decoy)),
            "ks_statistic": math.nan,
            "ks_pvalue": math.nan,
            "wasserstein_distance": math.nan,
        }
    ks = stats.ks_2samp(ent, decoy, alternative="two-sided", method="auto")
    return {
        "n_entrapment_scored": int(len(ent)),
        "n_decoy_scored": int(len(decoy)),
        "ks_statistic": float(ks.statistic),
        "ks_pvalue": float(ks.pvalue),
        "wasserstein_distance": float(stats.wasserstein_distance(ent, decoy)),
    }


def _class_summary(frame: pd.DataFrame, level: str) -> pd.DataFrame:
    grouping = ["truth_class"] if level == "global_peptide" else ["run_role", "truth_class"]
    out = frame.groupby(grouping, as_index=False).agg(
        assays=("precursor_id", "size"),
        assays_with_candidate=("candidate_rows", lambda x: int((x > 0).sum())),
        assays_with_scored_candidate=("scored_candidate_rows", lambda x: int((x > 0).sum())),
        median_candidate_rows=("candidate_rows", "median"),
        median_best_score=("score", "median"),
        median_reported_q=("reported_q", "median"),
    )
    out.insert(0, "level", level)
    return out


def _classify_export_results(
    results_path: Path,
    true_keys: set[tuple[str, int]],
    entrapment_keys: set[tuple[str, int]],
    requested_q_column: str | None,
) -> tuple[pd.DataFrame, str | None]:
    frame = _read_tsv(results_path)
    sequence_col = _first_existing(
        frame.columns,
        ("Sequence", "sequence", "PeptideSequence", "unmodified_sequence"),
    )
    charge_col = _first_existing(frame.columns, ("Charge", "charge", "PrecursorCharge"))
    run_col = _first_existing(frame.columns, ("run_name", "Run", "filename"))
    if sequence_col is None or charge_col is None:
        raise BenchmarkError(
            f"Could not classify OpenDIA.results.tsv; sequence/charge columns are missing. Columns={list(frame.columns)}"
        )

    work = frame.copy()
    work["sequence"] = _normalise_sequence(work[sequence_col])
    work["precursor_charge"] = pd.to_numeric(work[charge_col], errors="raise").astype(int)
    if run_col is not None:
        work["run_name"] = work[run_col].astype(str).map(lambda value: Path(value).stem)
        work["run_role"] = work["run_name"].map(_run_role)
    else:
        work["run_name"] = ""
        work["run_role"] = ""

    is_decoy = _normalise_decoy(work["decoy"]) if "decoy" in work.columns else pd.Series(False, index=work.index)
    keys = list(zip(work["sequence"], work["precursor_charge"]))
    classes = []
    for decoy, key in zip(is_decoy, keys):
        if decoy:
            classes.append("decoy")
        elif key in true_keys:
            classes.append("true_target")
        elif key in entrapment_keys:
            classes.append("entrapment")
        else:
            classes.append("unknown_target")
    work["truth_class"] = classes
    unknown = work[(work["truth_class"] == "unknown_target") & (~is_decoy)]
    if not unknown.empty:
        examples = unknown[["sequence", "precursor_charge"]].head(10).to_dict("records")
        raise BenchmarkError(
            f"OpenDIA.results.tsv contains {len(unknown)} target rows absent from truth; examples={examples}"
        )

    q_column = requested_q_column
    if q_column:
        if q_column not in work.columns:
            raise BenchmarkError(f"Requested export q-value column is absent: {q_column}")
    else:
        q_column = _first_existing(
            work.columns,
            (
                "m_score_peptide_global",
                "m_score_peptide_experiment_wide",
                "m_score",
                "ms2_m_score",
                "m_score_peptide_run_specific",
                "precursor_qvalue",
            ),
        )
    work["export_reported_q"] = (
        pd.to_numeric(work[q_column], errors="coerce") if q_column else math.nan
    )
    return work, q_column


def _accepted_metrics(frame: pd.DataFrame) -> dict[str, Any]:
    work = frame[frame["truth_class"].isin(["true_target", "entrapment"])].copy()
    tp = int((work["truth_class"] == "true_target").sum())
    fp = int((work["truth_class"] == "entrapment").sum())
    accepted = tp + fp
    return {
        "tp": tp,
        "fp_entrapment": fp,
        "accepted": accepted,
        "precision": tp / accepted if accepted else math.nan,
        # FDP is zero when the exported result set contains no discoveries.
        "empirical_fdp": fp / accepted if accepted else 0.0,
    }


def _export_context_metrics(
    export: pd.DataFrame,
    q_threshold: float,
    export_q_column: str | None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    groups = [("pooled", export)]
    run_roles = sorted(str(value) for value in export["run_role"].dropna().unique())
    groups.extend((role, export[export["run_role"] == role]) for role in run_roles)
    for scope, subset in groups:
        all_metrics = _accepted_metrics(subset)
        rows.append({"scope": scope, "selection": "all_exported_rows", **all_metrics})
        if export_q_column and subset["export_reported_q"].notna().any():
            filtered = subset[pd.to_numeric(subset["export_reported_q"], errors="coerce").le(q_threshold)]
            rows.append({
                "scope": scope,
                "selection": f"exported_rows_{export_q_column}_le_{q_threshold:g}",
                **_accepted_metrics(filtered),
            })

    target_export = export[export["truth_class"].isin(["true_target", "entrapment"])].copy()
    unique = target_export.drop_duplicates(["sequence", "precursor_charge"])
    unique_metrics = _accepted_metrics(unique)
    summary = {
        "all_rows_pooled": _accepted_metrics(target_export),
        "unique_precursors_pooled": unique_metrics,
        "detected_export_q_column": export_q_column,
        "export_rows": int(len(export)),
        "target_labelled_export_rows": int(len(target_export)),
    }
    return pd.DataFrame(rows), summary


def _save(path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


def _make_plots(
    run_best: pd.DataFrame,
    run_fdr: pd.DataFrame,
    global_fdr: pd.DataFrame,
    run_roc: pd.DataFrame,
    run_pr: pd.DataFrame,
    run_pep: pd.DataFrame,
    global_pep: pd.DataFrame,
    export_context: pd.DataFrame,
    q_threshold: float,
    out_dir: Path,
    dpi: int,
) -> list[Path]:
    fig_dir = out_dir / "figures"
    outputs: list[Path] = []

    path = fig_dir / "01_target_entrapment_decoy_score_distribution.png"
    plt.figure(figsize=(8.2, 5.5))
    finite_scores = pd.to_numeric(run_best["score"], errors="coerce").dropna()
    if not finite_scores.empty:
        bins = np.linspace(float(finite_scores.min()), float(finite_scores.max()), 45)
        for label in ("true_target", "entrapment", "decoy"):
            values = pd.to_numeric(
                run_best.loc[run_best["truth_class"] == label, "score"], errors="coerce"
            ).dropna()
            if len(values):
                plt.hist(
                    values,
                    bins=bins,
                    density=True,
                    histtype="step",
                    linewidth=1.8,
                    label=label.replace("_", " "),
                )
        plt.xlabel("Best run-level OpenDIA/Percolator discriminant score")
        plt.ylabel("Density")
        plt.title("True-target, entrapment and generated-decoy score distributions")
        plt.legend()
        _save(path, dpi)
        outputs.append(path)

    path = fig_dir / "02_empirical_fdr_calibration.png"
    run_pooled = run_fdr[run_fdr["scope"] == "pooled"].sort_values("q_threshold")
    global_pooled = global_fdr[global_fdr["scope"] == "pooled"].sort_values("q_threshold")
    values = pd.concat([run_pooled["empirical_fdp"], global_pooled["empirical_fdp"]], ignore_index=True)
    upper = max(0.10, float(values.dropna().max()) if values.notna().any() else 0.10)
    plt.figure(figsize=(6.6, 5.6))
    plt.plot([0, upper], [0, upper], "--", linewidth=1.1, label="ideal")
    plt.plot(run_pooled["q_threshold"], run_pooled["empirical_fdp"], marker="o", label="run-level pooled")
    plt.plot(global_pooled["q_threshold"], global_pooled["empirical_fdp"], marker="s", label="global peptide")
    run_roles = sorted(scope for scope in run_fdr["scope"].dropna().astype(str).unique() if scope != "pooled")
    # Individual-run curves are useful for small fixtures; for large studies only
    # show them when the legend remains readable. Per-run tables are always written.
    plot_run_roles = run_roles if len(run_roles) <= 12 else []
    for role in plot_run_roles:
        group = run_fdr[run_fdr["scope"] == role].sort_values("q_threshold")
        if not group.empty:
            plt.plot(group["q_threshold"], group["empirical_fdp"], marker=".", alpha=0.4, linewidth=0.8, label=f"run: {role.replace('_', ' ')}")
    export_pooled = export_context[
        (export_context["scope"] == "pooled") & (export_context["selection"] == "all_exported_rows")
    ]
    if not export_pooled.empty and np.isfinite(float(export_pooled.iloc[0]["empirical_fdp"])):
        plt.scatter([q_threshold], [float(export_pooled.iloc[0]["empirical_fdp"])], marker="x", s=70, label="actual final export")
    plt.xlabel("Nominal reported q-value threshold")
    plt.ylabel("Observed entrapment FDP")
    plt.title("FDR calibration by identification context")
    plt.xlim(0, upper)
    plt.ylim(0, upper)
    plt.legend(fontsize=8)
    _save(path, dpi)
    outputs.append(path)

    path = fig_dir / "03_entrapment_roc_curve.png"
    plt.figure(figsize=(6.0, 5.4))
    plt.plot(run_roc["fpr"], run_roc["tpr"])
    plt.plot([0, 1], [0, 1], "--", linewidth=1.0)
    plt.xlabel("False-positive rate (entrapment)")
    plt.ylabel("True-positive rate")
    plt.title("Run-level true-target versus entrapment ROC")
    _save(path, dpi)
    outputs.append(path)

    path = fig_dir / "04_entrapment_precision_recall_curve.png"
    plt.figure(figsize=(6.0, 5.4))
    plt.plot(run_pr["recall"], run_pr["precision"])
    prevalence = float(
        (run_best["truth_class"] == "true_target").sum()
        / run_best["truth_class"].isin(["true_target", "entrapment"]).sum()
    )
    plt.axhline(prevalence, linestyle="--", linewidth=1.0, label=f"prevalence={prevalence:.2f}")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Run-level true-target versus entrapment precision-recall")
    plt.legend()
    _save(path, dpi)
    outputs.append(path)

    if not run_pep.empty or not global_pep.empty:
        path = fig_dir / "05_pep_entrapment_calibration.png"
        plt.figure(figsize=(6.2, 5.4))
        upper_values = []
        for table in (run_pep, global_pep):
            if not table.empty:
                upper_values.extend(
                    table[["mean_reported_pep", "observed_entrapment_error_rate"]].to_numpy(dtype=float).ravel()
                )
        upper = max(0.10, float(np.nanmax(upper_values)))
        plt.plot([0, upper], [0, upper], "--", linewidth=1.0, label="ideal")
        if not run_pep.empty:
            plt.scatter(
                run_pep["mean_reported_pep"], run_pep["observed_entrapment_error_rate"],
                s=np.maximum(20, np.sqrt(run_pep["n"]) * 5), alpha=0.75, label="run-level PEP",
            )
        if not global_pep.empty:
            plt.scatter(
                global_pep["mean_reported_pep"], global_pep["observed_entrapment_error_rate"],
                s=np.maximum(20, np.sqrt(global_pep["n"]) * 5), marker="s", alpha=0.75, label="global peptide PEP",
            )
        plt.xlabel("Mean reported PEP")
        plt.ylabel("Observed entrapment identity-error rate")
        plt.title("PEP reliability against entrapment truth")
        plt.xlim(0, upper)
        plt.ylim(0, upper)
        plt.legend()
        _save(path, dpi)
        outputs.append(path)

    path = fig_dir / "06_fdr_context_comparison.png"
    context_rows = []
    run_at = run_pooled[np.isclose(run_pooled["q_threshold"], q_threshold)]
    global_at = global_pooled[np.isclose(global_pooled["q_threshold"], q_threshold)]
    if not run_at.empty:
        context_rows.append(("run-level q", float(run_at.iloc[0]["empirical_fdp"])))
    if not global_at.empty:
        context_rows.append(("global peptide q", float(global_at.iloc[0]["empirical_fdp"])))
    if not export_pooled.empty:
        context_rows.append(("actual export", float(export_pooled.iloc[0]["empirical_fdp"])))
    if context_rows:
        plt.figure(figsize=(6.5, 4.8))
        labels = [row[0] for row in context_rows]
        vals = [row[1] for row in context_rows]
        bars = plt.bar(labels, vals)
        plt.axhline(q_threshold, linestyle="--", linewidth=1.1, label=f"nominal {q_threshold:g}")
        plt.ylabel("Empirical entrapment FDP")
        plt.title(f"Identity-error context comparison at nominal q={q_threshold:g}")
        plt.xticks(rotation=15, ha="right")
        for bar, value in zip(bars, vals):
            plt.text(bar.get_x() + bar.get_width()/2, value, f"{100*value:.2f}%", ha="center", va="bottom")
        plt.legend()
        _save(path, dpi)
        outputs.append(path)

    return outputs


def _fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return "NA" if not np.isfinite(number) else f"{number:.{digits}f}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--opendia-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--workflow-oswpq", type=Path)
    parser.add_argument("--results-tsv", type=Path)
    parser.add_argument("--entrapment-truth", type=Path)

    parser.add_argument("--run-score-column", default="score_ms2_score")
    parser.add_argument("--run-q-column", default="score_ms2_qvalue")
    parser.add_argument("--run-pep-column", default="score_ms2_pep")
    parser.add_argument("--global-score-column", default="score_peptide_global_score")
    parser.add_argument("--global-q-column", default="score_peptide_global_qvalue")
    parser.add_argument("--global-pep-column", default="score_peptide_global_pep")
    parser.add_argument("--export-q-column")
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()

    build_dir = args.build_dir.expanduser().resolve()
    opendia_dir = args.opendia_dir.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    workflow = _find_workflow(opendia_dir, args.workflow_oswpq)
    results_path = _find_results(opendia_dir, args.results_tsv)

    true_keys, entrapment_keys, entrapment_mode, entrapment_truth_path = _load_truth(
        build_dir, args.entrapment_truth
    )
    precursors, runs, features = _load_workflow(workflow, true_keys, entrapment_keys)
    run_best, global_best = _prepare_scoring_tables(
        precursors,
        runs,
        features,
        run_score_column=args.run_score_column,
        run_q_column=args.run_q_column,
        run_pep_column=args.run_pep_column,
        global_score_column=args.global_score_column,
        global_q_column=args.global_q_column,
        global_pep_column=args.global_pep_column,
    )

    class_counts = precursors["truth_class"].value_counts().to_dict()
    if int(class_counts.get("true_target", 0)) != len(true_keys):
        raise BenchmarkError(
            f"Prepared library contains {class_counts.get('true_target', 0)} true targets; expected {len(true_keys)}"
        )
    if int(class_counts.get("entrapment", 0)) != len(entrapment_keys):
        raise BenchmarkError(
            f"Prepared library contains {class_counts.get('entrapment', 0)} entrapments; expected {len(entrapment_keys)}"
        )

    run_roc, run_pr, run_curve_metrics = _roc_pr_points(run_best)
    global_roc, global_pr, global_curve_metrics = _roc_pr_points(global_best)
    run_fdr = _fdr_calibration(run_best, DEFAULT_Q_GRID, level="run_level", include_run_scopes=True)
    global_fdr = _fdr_calibration(global_best, DEFAULT_Q_GRID, level="global_peptide", include_run_scopes=False)
    combined_fdr = pd.concat([run_fdr, global_fdr], ignore_index=True)
    run_pep, run_pep_metrics = _pep_calibration(run_best)
    global_pep, global_pep_metrics = _pep_calibration(global_best)
    null_metrics = _null_distribution_metrics(run_best)
    class_summary = pd.concat(
        [_class_summary(run_best, "run_level"), _class_summary(global_best, "global_peptide")],
        ignore_index=True,
    )

    run_pooled_q = _confusion_at_q(run_best, args.q_threshold)
    run_per_run = []
    run_roles = sorted(str(value) for value in run_best["run_role"].dropna().unique())
    for role in run_roles:
        subset = run_best[run_best["run_role"] == role]
        metrics = _confusion_at_q(subset, args.q_threshold)
        _, _, curve_metrics = _roc_pr_points(subset)
        metrics["run_role"] = role
        metrics["roc_auc"] = curve_metrics.roc_auc
        metrics["average_precision"] = curve_metrics.average_precision
        run_per_run.append(metrics)
    run_metrics = pd.DataFrame(run_per_run)
    global_q = _confusion_at_q(global_best, args.q_threshold)

    export, export_q_column = _classify_export_results(
        results_path, true_keys, entrapment_keys, args.export_q_column
    )
    export_context, export_summary = _export_context_metrics(
        export, args.q_threshold, export_q_column
    )

    context_rows = [
        {"context": "run_level", "scope": "pooled", **run_pooled_q},
        {"context": "global_peptide", "scope": "pooled", **global_q},
        {
            "context": "actual_export",
            "scope": "pooled",
            "q_threshold": math.nan,
            **export_summary["all_rows_pooled"],
        },
        {
            "context": "actual_export_unique_precursor",
            "scope": "pooled",
            "q_threshold": math.nan,
            **export_summary["unique_precursors_pooled"],
        },
    ]
    context_summary = pd.DataFrame(context_rows)

    # Output detailed tables.
    run_best.to_csv(out_dir / "entrapment_run_best_scores.tsv", sep="\t", index=False)
    global_best.to_csv(out_dir / "entrapment_global_peptide_scores.tsv", sep="\t", index=False)
    export.to_csv(out_dir / "actual_export_truth.tsv", sep="\t", index=False)
    export_context.to_csv(out_dir / "actual_export_fdp.tsv", sep="\t", index=False)
    context_summary.to_csv(out_dir / "fdr_context_summary.tsv", sep="\t", index=False)

    run_fdr.to_csv(out_dir / "run_level_empirical_fdr_calibration.tsv", sep="\t", index=False)
    global_fdr.to_csv(out_dir / "global_peptide_empirical_fdr_calibration.tsv", sep="\t", index=False)
    combined_fdr.to_csv(out_dir / "empirical_fdr_calibration.tsv", sep="\t", index=False)

    run_roc.to_csv(out_dir / "roc_curve_run_level.tsv", sep="\t", index=False)
    run_pr.to_csv(out_dir / "precision_recall_curve_run_level.tsv", sep="\t", index=False)
    global_roc.to_csv(out_dir / "roc_curve_global_peptide.tsv", sep="\t", index=False)
    global_pr.to_csv(out_dir / "precision_recall_curve_global_peptide.tsv", sep="\t", index=False)

    run_pep.to_csv(out_dir / "pep_calibration_run_level.tsv", sep="\t", index=False)
    global_pep.to_csv(out_dir / "pep_calibration_global_peptide.tsv", sep="\t", index=False)
    run_metrics.to_csv(out_dir / "entrapment_run_metrics.tsv", sep="\t", index=False)
    class_summary.to_csv(out_dir / "score_class_summary.tsv", sep="\t", index=False)

    plots = _make_plots(
        run_best,
        run_fdr,
        global_fdr,
        run_roc,
        run_pr,
        run_pep,
        global_pep,
        export_context,
        args.q_threshold,
        out_dir,
        args.dpi,
    )

    summary = {
        "workflow_oswpq": str(workflow),
        "results_tsv": str(results_path),
        "entrapment_truth_tsv": str(entrapment_truth_path),
        "entrapment_mode": entrapment_mode,
        "columns": {
            "run_score": args.run_score_column,
            "run_q": args.run_q_column,
            "run_pep": args.run_pep_column,
            "global_score": args.global_score_column,
            "global_q": args.global_q_column,
            "global_pep": args.global_pep_column,
            "export_q": export_q_column,
        },
        "prepared_library_class_counts": {key: int(value) for key, value in class_counts.items()},
        "operating_point_q_threshold": args.q_threshold,
        "run_level": {
            "roc_auc": run_curve_metrics.roc_auc,
            "average_precision": run_curve_metrics.average_precision,
            "pooled_operating_point": run_pooled_q,
            "per_run_operating_points": run_per_run,
        },
        "global_peptide": {
            "roc_auc": global_curve_metrics.roc_auc,
            "average_precision": global_curve_metrics.average_precision,
            "operating_point": global_q,
        },
        "actual_export": export_summary,
        "null_entrapment_vs_decoy_run_level": null_metrics,
        "pep_identity_error_calibration": {
            "run_level": run_pep_metrics,
            "global_peptide": global_pep_metrics,
        },
        "interpretation": (
            "Run-level and global-peptide q-values are evaluated in their own statistical contexts. "
            "Actual-export FDP classifies the rows OpenDIA really emitted and is not silently substituted for either q-value context. "
            "Observed entrapment FP/(true target TP + entrapment FP) is an empirical FDP for this realization."
        ),
        "figures": [str(path) for path in plots],
    }
    (out_dir / "entrapment_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    export_pooled = export_summary["all_rows_pooled"]
    report = [
        "# OpenDIA entrapment/FDR benchmark",
        "",
        f"Entrapment construction: **{entrapment_mode}**.",
        "",
        f"Prepared library truth: **{class_counts.get('true_target', 0)} true targets**, "
        f"**{class_counts.get('entrapment', 0)} known-absent entrapments**, and "
        f"**{class_counts.get('decoy', 0)} OpenDIA-generated decoys**.",
        "",
        "Entrapments are external target-labelled negatives. Generated decoys remain the internal null used by OpenDIA/Percolator.",
        "",
        "## Identity-error contexts",
        "",
        f"At nominal q <= **{args.q_threshold:g}**, the contexts are intentionally reported separately:",
        "",
        "| Context | TP | Entrapment FP | Accepted | Empirical FDP | Recall/precision note |",
        "|---|---:|---:|---:|---:|---|",
        f"| Run-level peak groups (pooled target/run observations) | {run_pooled_q['tp']} | {run_pooled_q['fp_entrapment']} | {run_pooled_q['accepted']} | {_fmt(run_pooled_q['empirical_fdp'])} | recall={_fmt(run_pooled_q['recall'])} |",
        f"| Global peptide identities | {global_q['tp']} | {global_q['fp_entrapment']} | {global_q['accepted']} | {_fmt(global_q['empirical_fdp'])} | recall={_fmt(global_q['recall'])} |",
        f"| Actual OpenDIA export rows | {export_pooled['tp']} | {export_pooled['fp_entrapment']} | {export_pooled['accepted']} | {_fmt(export_pooled['empirical_fdp'])} | final emitted rows |",
        "",
        f"The export q-value column detected in `OpenDIA.results.tsv` is **{export_q_column or 'none'}**. "
        "`actual_export_fdp.tsv` additionally reports the export after applying the configured q threshold when that column is available.",
        "",
        "## Discrimination",
        "",
        f"Run-level true-target versus entrapment ROC AUC: **{_fmt(run_curve_metrics.roc_auc)}**; average precision: **{_fmt(run_curve_metrics.average_precision)}**.",
        f"Global-peptide true-target versus entrapment ROC AUC: **{_fmt(global_curve_metrics.roc_auc)}**; average precision: **{_fmt(global_curve_metrics.average_precision)}**.",
        "",
        "## Null-model comparison",
        "",
        f"Run-level entrapment-vs-generated-decoy KS statistic: **{_fmt(null_metrics['ks_statistic'])}**; Wasserstein distance: **{_fmt(null_metrics['wasserstein_distance'])}**.",
        "",
        "## PEP calibration",
        "",
        f"Run-level Brier score: **{_fmt(run_pep_metrics['brier_score'])}**; ECE: **{_fmt(run_pep_metrics['expected_calibration_error'])}**.",
        f"Global-peptide Brier score: **{_fmt(global_pep_metrics['brier_score'])}**; ECE: **{_fmt(global_pep_metrics['expected_calibration_error'])}**.",
        "",
        "## Interpretation",
        "",
        (
            "The observed entrapment fraction is an empirical **false discovery proportion (FDP)** for this one simulation. "
            "For `independent` entrapments this is the canonical external-null calibration measurement; FDR is the expectation of FDP over repeated independent realizations. "
            "For `paired_hard`, the same arithmetic is reported only as an adversarial false-assignment stress metric because those negatives are intentionally colocated with real analyte coordinates and are not the canonical FDR null."
        ),
        "For calibration summaries, FDP is defined as `V / max(R, 1)`, so a threshold with zero discoveries contributes FDP = 0 rather than a missing value.",
        "",
        "The run-level q-value, global peptide q-value, and final export answer different questions and must not be interchanged. This benchmark therefore preserves all three contexts explicitly.",
        "",
        "## Outputs",
        "",
        "- `entrapment_run_best_scores.tsv`: best peak group per precursor/run.",
        "- `entrapment_global_peptide_scores.tsv`: one global identity per prepared precursor.",
        "- `actual_export_truth.tsv`: every OpenDIA result row classified by external truth.",
        "- `actual_export_fdp.tsv`: final-export FDP pooled and by run.",
        "- `fdr_context_summary.tsv`: side-by-side run/global/export comparison.",
        "- `run_level_empirical_fdr_calibration.tsv` and `global_peptide_empirical_fdr_calibration.tsv`.",
        "- `empirical_fdr_calibration.tsv`: combined calibration table with a `level` column.",
        "- run/global ROC, PR and PEP calibration TSVs.",
        "- `figures/01_target_entrapment_decoy_score_distribution.png`.",
        "- `figures/02_empirical_fdr_calibration.png`.",
        "- `figures/03_entrapment_roc_curve.png`.",
        "- `figures/04_entrapment_precision_recall_curve.png`.",
        "- `figures/05_pep_entrapment_calibration.png`.",
        "- `figures/06_fdr_context_comparison.png`.",
    ]
    (out_dir / "entrapment_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    print("Entrapment/FDR benchmark complete")
    print(f"Entrapment mode: {entrapment_mode}")
    print(
        f"Run-level q<={args.q_threshold:g}: TP={run_pooled_q['tp']}, "
        f"entrapment FP={run_pooled_q['fp_entrapment']}, FDP={run_pooled_q['empirical_fdp']:.4g}"
    )
    print(
        f"Global peptide q<={args.q_threshold:g}: TP={global_q['tp']}, "
        f"entrapment FP={global_q['fp_entrapment']}, FDP={global_q['empirical_fdp']:.4g}"
    )
    print(
        f"Actual export: TP={export_pooled['tp']}, entrapment FP={export_pooled['fp_entrapment']}, "
        f"FDP={export_pooled['empirical_fdp']:.4g}"
    )
    print(f"Run-level ROC AUC: {run_curve_metrics.roc_auc:.4f}")
    print(f"Global-peptide ROC AUC: {global_curve_metrics.roc_auc:.4f}")
    print(f"Report: {out_dir / 'entrapment_report.md'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BenchmarkError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
