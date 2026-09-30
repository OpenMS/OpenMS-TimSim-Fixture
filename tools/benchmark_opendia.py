#!/usr/bin/env python3
"""Benchmark OpenDIA output against TimSim synthetic ground truth.

The benchmark evaluates:
  * target precursor recovery per run;
  * retention-time and ion-mobility accuracy against realized TimSim signal;
  * recovery as a function of realized precursor signal;
  * control-replicate quantitative agreement;
  * peptide- and protein-level recovery of known treatment effects;
  * missingness and treatment-class detection coverage.

The transition library contains predictor RT/IM coordinates, while TimSim stores
the realized chromatographic and mobility distributions separately. By default
this benchmark reads each run's ``synthetic_data.db`` and evaluates coordinate
accuracy against the realized signal apex. Peptide-level event counts are retained
for fold-change truth, but they are not treated as precursor-level detectability.

The three-run fixture contains two controls and one treatment. It is suitable for
coordinate, recovery, and fold-change validation, but not for calibrated
statistical differential-expression testing because the treatment group has one
replicate.
"""
from __future__ import annotations

import argparse
import json
import math
import sqlite3
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

RUN_ROLES = ("baseline", "control_replicate", "treatment")
DEFAULT_RUN_NAMES = {
    "baseline": "OpenSwathTimSim_01_baseline",
    "control_replicate": "OpenSwathTimSim_02_biological_replicate",
    "treatment": "OpenSwathTimSim_03_treatment",
}


@dataclass(frozen=True)
class BenchmarkConfig:
    build_dir: Path
    opendia_dir: Path
    output_dir: Path
    results_tsv: Path | None = None
    peptide_matrix_tsv: Path | None = None
    protein_matrix_tsv: Path | None = None
    precursor_truth_tsv: Path | None = None
    condition_truth_tsv: Path | None = None
    synthetic_db_root: Path | None = None
    realized_coordinate: str = "apex"
    use_realized_truth: bool = True
    q_threshold: float = 0.01
    q_column: str | None = None
    intensity_column: str = "Intensity"
    pseudocount: float = 1.0
    fold_change_threshold: float = 0.5
    baseline_run: str = DEFAULT_RUN_NAMES["baseline"]
    control_replicate_run: str = DEFAULT_RUN_NAMES["control_replicate"]
    treatment_run: str = DEFAULT_RUN_NAMES["treatment"]
    dpi: int = 170


@dataclass(frozen=True)
class MetricResult:
    n: int
    bias: float
    median_error: float
    mae: float
    rmse: float
    pearson_r: float
    spearman_rho: float
    r_squared: float
    p95_absolute_error: float


def _finite_pairs(x: Iterable[Any], y: Iterable[Any]) -> tuple[np.ndarray, np.ndarray]:
    x_array = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    y_array = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(x_array) & np.isfinite(y_array)
    return x_array[mask], y_array[mask]


def _safe_correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 2 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return math.nan
    if method == "pearson":
        return float(stats.pearsonr(x, y).statistic)
    return float(stats.spearmanr(x, y).statistic)


def calculate_metrics(expected: Iterable[Any], observed: Iterable[Any]) -> MetricResult:
    expected_array, observed_array = _finite_pairs(expected, observed)
    if len(expected_array) == 0:
        return MetricResult(0, *(math.nan for _ in range(8)))
    errors = observed_array - expected_array
    denominator = float(np.sum((expected_array - np.mean(expected_array)) ** 2))
    r_squared = (
        1.0 - float(np.sum(errors ** 2)) / denominator
        if denominator > 0
        else math.nan
    )
    return MetricResult(
        n=int(len(expected_array)),
        bias=float(np.mean(errors)),
        median_error=float(np.median(errors)),
        mae=float(np.mean(np.abs(errors))),
        rmse=float(np.sqrt(np.mean(errors ** 2))),
        pearson_r=_safe_correlation(expected_array, observed_array, "pearson"),
        spearman_rho=_safe_correlation(expected_array, observed_array, "spearman"),
        r_squared=r_squared,
        p95_absolute_error=float(np.quantile(np.abs(errors), 0.95)),
    )


def _read_tsv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required TSV does not exist: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _first_existing(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    available = set(columns)
    return next((candidate for candidate in candidates if candidate in available), None)


def _normalise_decoy(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.lower()
    return text.isin({"1", "true", "yes", "decoy"}) | (_numeric(series).fillna(0) != 0)


def _resolve_paths(config: BenchmarkConfig) -> dict[str, Path]:
    return {
        "results": config.results_tsv or config.opendia_dir / "OpenDIA.results.tsv",
        "peptide_matrix": config.peptide_matrix_tsv or config.opendia_dir / "OpenDIA.peptide.matrix.tsv",
        "protein_matrix": config.protein_matrix_tsv or config.opendia_dir / "OpenDIA.protein.matrix.tsv",
        "precursor_truth": config.precursor_truth_tsv or config.build_dir / "OpenSwathTimSim.ground_truth.tsv",
        "condition_truth": config.condition_truth_tsv or config.build_dir / "OpenSwathTimSim.condition_truth.tsv",
    }


def _run_mapping(config: BenchmarkConfig) -> dict[str, str]:
    return {
        "baseline": config.baseline_run,
        "control_replicate": config.control_replicate_run,
        "treatment": config.treatment_run,
    }


def _detect_q_column(results: pd.DataFrame, requested: str | None) -> str | None:
    if requested:
        if requested.lower() in {"none", "off", "disable"}:
            return None
        if requested not in results.columns:
            raise ValueError(f"Requested q-value column is absent: {requested}")
        return requested
    return _first_existing(
        results.columns,
        (
            "m_score",
            "ms2_m_score",
            "m_score_peptide_run_specific",
            "m_score_peptide_experiment_wide",
            "m_score_peptide_global",
            "precursor_pep",
            "pep",
        ),
    )


def prepare_results(
    results: pd.DataFrame,
    config: BenchmarkConfig,
) -> tuple[pd.DataFrame, str | None, dict[str, int]]:
    required = {"run_name", "Sequence", "Charge", "RT", "EXP_IM", config.intensity_column}
    missing = required - set(results.columns)
    if missing:
        raise ValueError(f"OpenDIA results are missing required columns: {sorted(missing)}")

    work = results.copy()
    work["run_name"] = work["run_name"].astype(str)
    work["sequence"] = work["Sequence"].astype(str).str.strip()
    work["charge"] = _numeric(work["Charge"]).astype("Int64")
    work["intensity"] = _numeric(work[config.intensity_column])
    work["observed_rt"] = _numeric(work["RT"])
    work["observed_im"] = _numeric(work["EXP_IM"])
    work["d_score_numeric"] = _numeric(work["d_score"]) if "d_score" in work else math.nan
    work["is_decoy"] = _normalise_decoy(work["decoy"]) if "decoy" in work else False

    q_column = _detect_q_column(work, config.q_column)
    if q_column:
        work["q_value"] = _numeric(work[q_column])
    else:
        work["q_value"] = math.nan

    run_names = set(_run_mapping(config).values())
    unknown_runs = sorted(set(work["run_name"]) - run_names)
    work = work[work["run_name"].isin(run_names)].copy()

    target_rows = work[~work["is_decoy"]].copy()
    decoy_rows = work[work["is_decoy"]].copy()
    counts = {
        "input_rows": int(len(results)),
        "target_rows_before_q_filter": int(len(target_rows)),
        "decoy_rows_before_q_filter": int(len(decoy_rows)),
        "unknown_run_rows_ignored": int(len(results) - len(work)),
    }

    if q_column and target_rows["q_value"].notna().any():
        target_rows = target_rows[target_rows["q_value"].notna() & (target_rows["q_value"] <= config.q_threshold)]
        decoy_rows = decoy_rows[decoy_rows["q_value"].notna() & (decoy_rows["q_value"] <= config.q_threshold)]
    counts["target_rows_after_q_filter"] = int(len(target_rows))
    counts["decoy_rows_after_q_filter"] = int(len(decoy_rows))

    sort_columns = ["run_name", "sequence", "charge", "q_value", "d_score_numeric", "intensity"]
    target_rows = target_rows.sort_values(
        sort_columns,
        ascending=[True, True, True, True, False, False],
        na_position="last",
    )
    best = target_rows.drop_duplicates(["run_name", "sequence", "charge"], keep="first")
    selected_columns = [
        "run_name", "sequence", "charge", "intensity", "observed_rt", "observed_im",
        "q_value", "d_score_numeric",
    ]
    for optional in ("feature_id", "ProteinName", "FullPeptideName", "mz", "delta_rt", "delta_iRT"):
        if optional in best.columns:
            selected_columns.append(optional)
    if unknown_runs:
        print(f"WARNING: Ignored result rows from unconfigured runs: {unknown_runs}", file=sys.stderr)
    return best[selected_columns].copy(), q_column, counts



def _parse_numeric_array(value: Any) -> np.ndarray:
    """Parse TimSim's JSON-like numeric array columns."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return np.empty(0, dtype=float)
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    text = str(value).strip()
    if not text:
        return np.empty(0, dtype=float)
    try:
        return np.asarray(json.loads(text), dtype=float).reshape(-1)
    except json.JSONDecodeError:
        stripped = text.strip("[]")
        if not stripped:
            return np.empty(0, dtype=float)
        return np.asarray(
            [float(token) for token in stripped.replace(",", " ").split()],
            dtype=float,
        )


def _realized_coordinate(
    ids: np.ndarray,
    weights: np.ndarray,
    coordinate_lookup: dict[int, float],
) -> tuple[float, float, int]:
    """Return weighted apex, centroid, and retained point count."""
    if ids.size == 0 or weights.size == 0 or ids.size != weights.size:
        return math.nan, math.nan, 0

    coordinates = np.asarray(
        [coordinate_lookup.get(int(identifier), math.nan) for identifier in ids],
        dtype=float,
    )
    valid = np.isfinite(coordinates) & np.isfinite(weights)
    if not valid.any():
        return math.nan, math.nan, 0

    coordinates = coordinates[valid]
    retained_weights = weights[valid]
    if retained_weights.size == 0:
        return math.nan, math.nan, 0

    apex = float(coordinates[int(np.argmax(retained_weights))])
    total = float(retained_weights.sum())
    centroid = (
        float(np.average(coordinates, weights=retained_weights))
        if total > 0.0
        else math.nan
    )
    return apex, centroid, int(len(coordinates))


def _load_realized_truth_for_run(
    db_path: Path,
    truth: pd.DataFrame,
    *,
    run_role: str,
    run_name: str,
) -> pd.DataFrame:
    """Read precursor-level realized signal information from one TimSim database."""
    if not db_path.is_file():
        raise FileNotFoundError(f"TimSim synthetic database does not exist: {db_path}")

    connection = sqlite3.connect(db_path)
    try:
        peptides = pd.read_sql_query(
            """
            SELECT
                peptide_id,
                sequence,
                events,
                retention_time_gru_predictor,
                rt_mu,
                rt_sigma,
                rt_lambda,
                frame_occurrence,
                frame_abundance
            FROM peptides
            """,
            connection,
        )
        ions = pd.read_sql_query(
            """
            SELECT
                ion_id,
                peptide_id,
                charge,
                mz,
                relative_abundance,
                inv_mobility_gru_predictor,
                inv_mobility_gru_predictor_std,
                scan_occurrence,
                scan_abundance
            FROM ions
            """,
            connection,
        )
        frames = pd.read_sql_query("SELECT frame_id, time FROM frames", connection)
        scans = pd.read_sql_query("SELECT scan, mobility FROM scans", connection)
    finally:
        connection.close()

    peptides["sequence"] = peptides["sequence"].astype(str).str.strip()
    ions["charge"] = _numeric(ions["charge"]).astype("Int64")
    frame_time = {
        int(row.frame_id): float(row.time)
        for row in frames.itertuples(index=False)
    }
    scan_im = {
        int(row.scan): float(row.mobility)
        for row in scans.itertuples(index=False)
    }

    selected = truth[["sequence", "charge", "precursor_key"]].drop_duplicates().copy()
    selected = selected.merge(peptides, on="sequence", how="left", validate="many_to_one")
    selected = selected.merge(
        ions,
        on=["peptide_id", "charge"],
        how="left",
        validate="one_to_one",
        suffixes=("_peptide", "_ion"),
    )

    missing = selected["ion_id"].isna()
    if missing.any():
        examples = selected.loc[missing, ["sequence", "charge"]].head(10).to_dict("records")
        raise ValueError(
            f"{int(missing.sum())} selected precursor groups are absent from {db_path}; "
            f"examples: {examples}"
        )

    rows: list[dict[str, Any]] = []
    for row in selected.itertuples(index=False):
        frame_ids = _parse_numeric_array(row.frame_occurrence).astype(int)
        frame_abundance = _parse_numeric_array(row.frame_abundance)
        scan_ids = _parse_numeric_array(row.scan_occurrence).astype(int)
        scan_abundance = _parse_numeric_array(row.scan_abundance)

        frame_sum = float(frame_abundance.sum()) if frame_abundance.size else 0.0
        scan_sum = float(scan_abundance.sum()) if scan_abundance.size else 0.0
        rt_apex, rt_centroid, frame_points = _realized_coordinate(
            frame_ids, frame_abundance, frame_time
        )
        im_apex, im_centroid, scan_points = _realized_coordinate(
            scan_ids, scan_abundance, scan_im
        )

        peptide_events = float(row.events)
        ion_fraction = float(row.relative_abundance)
        realized_event_proxy = (
            peptide_events * frame_sum * ion_fraction * scan_sum
        )

        rows.append(
            {
                "run_role": run_role,
                "run_name": run_name,
                "precursor_key": row.precursor_key,
                "sequence": row.sequence,
                "charge": int(row.charge),
                "peptide_id": int(row.peptide_id),
                "ion_id": int(row.ion_id),
                "simulated_peptide_events": peptide_events,
                "ion_relative_abundance": ion_fraction,
                "frame_abundance_sum": frame_sum,
                "scan_abundance_sum": scan_sum,
                "realized_event_proxy": realized_event_proxy,
                "db_predictor_rt": float(row.retention_time_gru_predictor),
                "rt_mu": float(row.rt_mu),
                "rt_sigma": float(row.rt_sigma),
                "rt_lambda": float(row.rt_lambda),
                "frame_points": frame_points,
                "realized_rt_apex": rt_apex,
                "realized_rt_centroid": rt_centroid,
                "db_predictor_im": float(row.inv_mobility_gru_predictor),
                "db_predictor_im_std": float(row.inv_mobility_gru_predictor_std),
                "scan_points": scan_points,
                "realized_im_apex": im_apex,
                "realized_im_centroid": im_centroid,
            }
        )

    return pd.DataFrame(rows)


def build_realized_truth(
    truth: pd.DataFrame,
    config: BenchmarkConfig,
    run_mapping: dict[str, str],
) -> pd.DataFrame:
    """Build run-specific realized precursor truth from TimSim synthetic_data.db files."""
    if not config.use_realized_truth:
        return pd.DataFrame()

    root = config.synthetic_db_root or config.build_dir
    frames: list[pd.DataFrame] = []
    for role, run_name in run_mapping.items():
        db_path = root / run_name / "synthetic_data.db"
        frames.append(
            _load_realized_truth_for_run(
                db_path,
                truth,
                run_role=role,
                run_name=run_name,
            )
        )
    return pd.concat(frames, ignore_index=True)


def prepare_truth(precursor_truth: pd.DataFrame, condition_truth: pd.DataFrame) -> pd.DataFrame:
    required_precursor = {
        "TransitionGroupId", "PeptideSequence", "PrecursorCharge",
        "PrecursorMz", "PrecursorIonMobility", "RetentionTimeSeconds",
    }
    missing = required_precursor - set(precursor_truth.columns)
    if missing:
        raise ValueError(f"Precursor truth is missing columns: {sorted(missing)}")

    required_condition = {
        "PeptideSequence", "ProteinId", "TreatmentClass", "BaselineEvents",
        "ControlReplicateEvents", "TreatmentEvents",
    }
    missing = required_condition - set(condition_truth.columns)
    if missing:
        raise ValueError(f"Condition truth is missing columns: {sorted(missing)}")

    precursor = precursor_truth.copy()
    precursor["sequence"] = precursor["PeptideSequence"].astype(str).str.strip()
    precursor["charge"] = _numeric(precursor["PrecursorCharge"]).astype("Int64")
    # These are assay/predictor coordinates supplied to OpenDIA, not necessarily
    # the realized signal apex in the simulated frames/scans.
    precursor["assay_rt"] = _numeric(precursor["RetentionTimeSeconds"])
    precursor["assay_im"] = _numeric(precursor["PrecursorIonMobility"])
    # Expose stable assay-coordinate aliases used by benchmark outputs.
    precursor["truth_rt"] = precursor["assay_rt"]
    precursor["truth_im"] = precursor["assay_im"]
    precursor["truth_precursor_mz"] = _numeric(precursor["PrecursorMz"])

    condition = condition_truth.copy()
    condition["sequence"] = condition["PeptideSequence"].astype(str).str.strip()
    for column in (
        "BaselineEvents", "ControlReplicateEvents", "TreatmentEvents",
        "ControlReplicateFactor", "TreatmentBiologicalFactor", "TreatmentNoiseFactor",
        "TreatmentTotalFactor",
    ):
        if column in condition:
            condition[column] = _numeric(condition[column])

    duplicated = condition[condition.duplicated("sequence", keep=False)]
    if not duplicated.empty:
        conflicting = duplicated.groupby("sequence")["ProteinId"].nunique()
        if (conflicting > 1).any():
            raise ValueError("Condition truth contains peptide sequences assigned to multiple proteins")
        condition = condition.drop_duplicates("sequence", keep="first")

    columns = [
        "sequence", "ProteinId", "TreatmentClass", "BaselineEvents",
        "ControlReplicateEvents", "TreatmentEvents",
    ]
    for optional in (
        "ControlReplicateFactor", "TreatmentBiologicalFactor", "TreatmentNoiseFactor",
        "TreatmentTotalFactor",
    ):
        if optional in condition.columns:
            columns.append(optional)

    truth = precursor.merge(condition[columns], on="sequence", how="left", validate="many_to_one")
    missing_condition = truth["ProteinId"].isna().sum()
    if missing_condition:
        raise ValueError(f"{missing_condition} selected precursor groups lack condition truth")
    truth["precursor_key"] = truth["sequence"] + "/" + truth["charge"].astype(str)
    return truth


def build_precursor_long(
    truth: pd.DataFrame,
    best_results: pd.DataFrame,
    run_mapping: dict[str, str],
    realized_truth: pd.DataFrame,
    config: BenchmarkConfig,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []

    coordinate_suffix = config.realized_coordinate
    if coordinate_suffix not in {"apex", "centroid"}:
        raise ValueError("realized_coordinate must be 'apex' or 'centroid'")

    for role, run_name in run_mapping.items():
        run = best_results[best_results["run_name"] == run_name].copy()
        merged = truth.merge(run, on=["sequence", "charge"], how="left", validate="one_to_one")

        if not realized_truth.empty:
            realized_run = realized_truth[realized_truth["run_role"] == role].copy()
            realized_columns = [
                column
                for column in realized_run.columns
                if column not in {"run_role", "run_name", "sequence", "charge"}
            ]
            merged = merged.merge(
                realized_run[["precursor_key", *[c for c in realized_columns if c != "precursor_key"]]],
                on="precursor_key",
                how="left",
                validate="one_to_one",
            )

        merged["run_role"] = role
        merged["configured_run_name"] = run_name
        merged["detected"] = merged["intensity"].notna() & (merged["intensity"] > 0)

        rt_column = f"realized_rt_{coordinate_suffix}"
        im_column = f"realized_im_{coordinate_suffix}"
        if rt_column in merged:
            merged["evaluation_rt"] = _numeric(merged[rt_column]).where(
                _numeric(merged[rt_column]).notna(), merged["assay_rt"]
            )
            merged["evaluation_im"] = _numeric(merged[im_column]).where(
                _numeric(merged[im_column]).notna(), merged["assay_im"]
            )
            merged["coordinate_truth_source"] = f"realized_{coordinate_suffix}"
        else:
            merged["evaluation_rt"] = merged["assay_rt"]
            merged["evaluation_im"] = merged["assay_im"]
            merged["coordinate_truth_source"] = "assay_predictor"

        merged["rt_error_seconds"] = merged["observed_rt"] - merged["evaluation_rt"]
        merged["im_error"] = merged["observed_im"] - merged["evaluation_im"]
        merged["assay_rt_error_seconds"] = merged["observed_rt"] - merged["assay_rt"]
        merged["assay_im_error"] = merged["observed_im"] - merged["assay_im"]

        if "realized_rt_apex" in merged:
            merged["assay_to_realized_rt_apex"] = merged["realized_rt_apex"] - merged["assay_rt"]
            merged["assay_to_realized_im_apex"] = merged["realized_im_apex"] - merged["assay_im"]

        rows.append(merged)
    return pd.concat(rows, ignore_index=True)


def _run_metrics(precursor_long: pd.DataFrame, q_column: str | None) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for role in RUN_ROLES:
        subset = precursor_long[precursor_long["run_role"] == role]
        detected = subset[subset["detected"]]

        rt = calculate_metrics(detected["evaluation_rt"], detected["observed_rt"])
        im = calculate_metrics(detected["evaluation_im"], detected["observed_im"])
        assay_rt = calculate_metrics(detected["assay_rt"], detected["observed_rt"])
        assay_im = calculate_metrics(detected["assay_im"], detected["observed_im"])

        row: dict[str, Any] = {
            "run_role": role,
            "run_name": subset["configured_run_name"].iloc[0],
            "coordinate_truth_source": subset["coordinate_truth_source"].iloc[0],
            "target_precursors": int(len(subset)),
            "detected_precursors": int(subset["detected"].sum()),
            "recovery_rate": float(subset["detected"].mean()),
            "median_intensity": float(detected["intensity"].median()) if len(detected) else math.nan,
            "q_column": q_column or "",
            "median_q_value": float(detected["q_value"].median()) if detected["q_value"].notna().any() else math.nan,
            "median_realized_event_proxy": (
                float(subset["realized_event_proxy"].median())
                if "realized_event_proxy" in subset
                else math.nan
            ),
        }
        row.update({f"rt_{key}": value for key, value in asdict(rt).items()})
        row.update({f"im_{key}": value for key, value in asdict(im).items()})
        row.update({f"assay_rt_{key}": value for key, value in asdict(assay_rt).items()})
        row.update({f"assay_im_{key}": value for key, value in asdict(assay_im).items()})
        rows.append(row)
    return pd.DataFrame(rows)


def build_signal_recovery_summary(precursor_long: pd.DataFrame) -> pd.DataFrame:
    """Summarize detection recovery across realized-signal quartiles for each run."""
    if "realized_event_proxy" not in precursor_long:
        return pd.DataFrame()

    rows: list[dict[str, Any]] = []
    for role in RUN_ROLES:
        subset = precursor_long[precursor_long["run_role"] == role].copy()
        valid = subset["realized_event_proxy"].notna()
        if valid.sum() < 4:
            continue

        # Rank-based qcut avoids duplicate-bin-edge failures for tied TimSim values.
        subset.loc[valid, "realized_signal_quartile"] = pd.qcut(
            subset.loc[valid, "realized_event_proxy"].rank(method="first"),
            4,
            labels=["Q1", "Q2", "Q3", "Q4"],
        ).astype(str)

        for quartile in ("Q1", "Q2", "Q3", "Q4"):
            group = subset[subset["realized_signal_quartile"] == quartile]
            if group.empty:
                continue
            rows.append(
                {
                    "run_role": role,
                    "run_name": group["configured_run_name"].iloc[0],
                    "realized_signal_quartile": quartile,
                    "targets": int(len(group)),
                    "detected": int(group["detected"].sum()),
                    "recovery_rate": float(group["detected"].mean()),
                    "median_realized_event_proxy": float(group["realized_event_proxy"].median()),
                    "median_peptide_events": float(group["simulated_peptide_events"].median()),
                    "median_frame_abundance_sum": float(group["frame_abundance_sum"].median()),
                    "median_ion_relative_abundance": float(group["ion_relative_abundance"].median()),
                    "median_intensity": float(group["intensity"].median()),
                    "median_d_score": float(group["d_score_numeric"].median()),
                }
            )
    return pd.DataFrame(rows)



def _log2_ratio(numerator: pd.Series, denominator: pd.Series, pseudocount: float) -> pd.Series:
    return np.log2(_numeric(numerator) + pseudocount) - np.log2(_numeric(denominator) + pseudocount)


def build_peptide_quant_benchmark(
    precursor_long: pd.DataFrame,
    config: BenchmarkConfig,
) -> pd.DataFrame:
    metadata_columns = [
        "precursor_key", "sequence", "charge", "ProteinId", "TreatmentClass",
        "BaselineEvents", "ControlReplicateEvents", "TreatmentEvents",
    ]
    metadata = precursor_long[metadata_columns].drop_duplicates("precursor_key", keep="first")
    intensity = precursor_long.pivot(
        index="precursor_key", columns="run_role", values="intensity"
    ).reset_index()
    pivot = metadata.merge(intensity, on="precursor_key", how="left", validate="one_to_one")
    for role in RUN_ROLES:
        if role not in pivot.columns:
            pivot[role] = math.nan

    if "realized_event_proxy" in precursor_long.columns:
        realized = precursor_long.pivot(
            index="precursor_key",
            columns="run_role",
            values="realized_event_proxy",
        ).rename(columns={role: f"{role}_realized_event_proxy" for role in RUN_ROLES})
        pivot = pivot.merge(realized.reset_index(), on="precursor_key", how="left", validate="one_to_one")

    pivot["complete_controls"] = pivot[["baseline", "control_replicate"]].notna().all(axis=1)
    pivot["complete_all_runs"] = pivot[list(RUN_ROLES)].notna().all(axis=1)
    pivot["observed_control_log2fc"] = _log2_ratio(
        pivot["control_replicate"], pivot["baseline"], config.pseudocount
    )
    pivot["expected_control_log2fc"] = np.log2(
        pivot["ControlReplicateEvents"] / pivot["BaselineEvents"]
    )
    pivot["control_log2fc_residual"] = (
        pivot["observed_control_log2fc"] - pivot["expected_control_log2fc"]
    )

    pivot["observed_treatment_log2fc"] = (
        np.log2(pivot["treatment"] + config.pseudocount)
        - 0.5 * (
            np.log2(pivot["baseline"] + config.pseudocount)
            + np.log2(pivot["control_replicate"] + config.pseudocount)
        )
    )
    pivot["expected_treatment_log2fc"] = (
        np.log2(pivot["TreatmentEvents"])
        - 0.5 * (
            np.log2(pivot["BaselineEvents"])
            + np.log2(pivot["ControlReplicateEvents"])
        )
    )
    pivot["treatment_log2fc_residual"] = (
        pivot["observed_treatment_log2fc"] - pivot["expected_treatment_log2fc"]
    )
    pivot["predicted_class"] = classify_fold_change(
        pivot["observed_treatment_log2fc"], config.fold_change_threshold
    )
    return pivot


def _matrix_run_columns(matrix: pd.DataFrame, run_mapping: dict[str, str]) -> dict[str, str]:
    columns = set(matrix.columns)
    missing = [name for name in run_mapping.values() if name not in columns]
    if missing:
        raise ValueError(f"Matrix is missing configured run columns: {missing}")
    return run_mapping


def build_protein_quant_benchmark(
    protein_matrix: pd.DataFrame,
    condition_truth: pd.DataFrame,
    config: BenchmarkConfig,
    run_mapping: dict[str, str],
) -> pd.DataFrame:
    protein_column = _first_existing(protein_matrix.columns, ("ProteinName", "ProteinId", "protein"))
    if not protein_column:
        raise ValueError("Could not identify the protein identifier column in OpenDIA.protein.matrix.tsv")
    _matrix_run_columns(protein_matrix, run_mapping)

    matrix = protein_matrix.copy()
    matrix["ProteinId"] = matrix[protein_column].astype(str)
    for role, run_name in run_mapping.items():
        matrix[role] = _numeric(matrix[run_name])

    truth = condition_truth.copy()
    for column in ("BaselineEvents", "ControlReplicateEvents", "TreatmentEvents"):
        truth[column] = _numeric(truth[column])
    class_counts = truth.groupby("ProteinId")["TreatmentClass"].nunique()
    if (class_counts > 1).any():
        raise ValueError("Treatment class is inconsistent within at least one protein")
    protein_truth = truth.groupby("ProteinId", as_index=False).agg(
        TreatmentClass=("TreatmentClass", "first"),
        BaselineEvents=("BaselineEvents", "sum"),
        ControlReplicateEvents=("ControlReplicateEvents", "sum"),
        TreatmentEvents=("TreatmentEvents", "sum"),
        PeptidesInTruth=("PeptideSequence", "nunique"),
    )

    merged = protein_truth.merge(matrix[["ProteinId", *RUN_ROLES]], on="ProteinId", how="left")
    merged["complete_controls"] = merged[["baseline", "control_replicate"]].notna().all(axis=1)
    merged["complete_all_runs"] = merged[list(RUN_ROLES)].notna().all(axis=1)
    merged["observed_control_log2fc"] = _log2_ratio(
        merged["control_replicate"], merged["baseline"], config.pseudocount
    )
    merged["expected_control_log2fc"] = np.log2(
        merged["ControlReplicateEvents"] / merged["BaselineEvents"]
    )
    merged["observed_treatment_log2fc"] = (
        np.log2(merged["treatment"] + config.pseudocount)
        - 0.5 * (
            np.log2(merged["baseline"] + config.pseudocount)
            + np.log2(merged["control_replicate"] + config.pseudocount)
        )
    )
    merged["expected_treatment_log2fc"] = (
        np.log2(merged["TreatmentEvents"])
        - 0.5 * (
            np.log2(merged["BaselineEvents"])
            + np.log2(merged["ControlReplicateEvents"])
        )
    )
    merged["treatment_log2fc_residual"] = (
        merged["observed_treatment_log2fc"] - merged["expected_treatment_log2fc"]
    )
    merged["predicted_class"] = classify_fold_change(
        merged["observed_treatment_log2fc"], config.fold_change_threshold
    )
    return merged


def classify_fold_change(values: pd.Series, threshold: float) -> pd.Series:
    numeric = _numeric(values)
    result = pd.Series(pd.NA, index=values.index, dtype="string")
    result.loc[numeric >= threshold] = "up"
    result.loc[numeric <= -threshold] = "down"
    result.loc[numeric.notna() & (numeric.abs() < threshold)] = "unchanged"
    return result


def _binary_auc(labels: pd.Series, scores: pd.Series) -> float:
    label_array = labels.astype(bool).to_numpy()
    score_array = _numeric(scores).to_numpy(dtype=float)
    mask = np.isfinite(score_array)
    label_array = label_array[mask]
    score_array = score_array[mask]
    positives = int(label_array.sum())
    negatives = int((~label_array).sum())
    if positives == 0 or negatives == 0:
        return math.nan
    ranks = stats.rankdata(score_array, method="average")
    rank_sum = float(ranks[label_array].sum())
    return (rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def classification_metrics(table: pd.DataFrame) -> tuple[dict[str, float], pd.DataFrame]:
    complete = table.dropna(subset=["TreatmentClass", "predicted_class"]).copy()
    labels = ["down", "unchanged", "up"]
    confusion = pd.crosstab(
        complete["TreatmentClass"],
        complete["predicted_class"],
        rownames=["truth"],
        colnames=["predicted"],
        dropna=False,
    ).reindex(index=labels, columns=labels, fill_value=0)
    total = int(confusion.to_numpy().sum())
    accuracy = float(np.trace(confusion.to_numpy()) / total) if total else math.nan
    f1_scores: list[float] = []
    for label in labels:
        tp = float(confusion.loc[label, label])
        fp = float(confusion[label].sum() - tp)
        fn = float(confusion.loc[label].sum() - tp)
        precision = tp / (tp + fp) if tp + fp else math.nan
        recall = tp / (tp + fn) if tp + fn else math.nan
        if np.isfinite(precision) and np.isfinite(recall) and precision + recall:
            f1_scores.append(2.0 * precision * recall / (precision + recall))
    changed = complete[complete["TreatmentClass"].isin(["up", "down"])]
    direction_accuracy = (
        float((changed["TreatmentClass"] == changed["predicted_class"]).mean())
        if len(changed)
        else math.nan
    )
    altered_auc = _binary_auc(
        complete["TreatmentClass"].ne("unchanged"),
        complete["observed_treatment_log2fc"].abs(),
    )
    return {
        "n_classified": int(len(complete)),
        "three_class_accuracy": accuracy,
        "macro_f1": float(np.mean(f1_scores)) if f1_scores else math.nan,
        "changed_direction_accuracy": direction_accuracy,
        "altered_vs_unchanged_auc": altered_auc,
    }, confusion.reset_index()


def _quant_summary(peptide: pd.DataFrame, protein: pd.DataFrame) -> dict[str, Any]:
    control_peptide = peptide[peptide["complete_controls"]]
    treatment_peptide = peptide[peptide["complete_all_runs"]]
    control_protein = protein[protein["complete_controls"]]
    treatment_protein = protein[protein["complete_all_runs"]]

    summaries: dict[str, Any] = {
        "peptide_control": asdict(calculate_metrics(
            control_peptide["expected_control_log2fc"],
            control_peptide["observed_control_log2fc"],
        )),
        "peptide_treatment": asdict(calculate_metrics(
            treatment_peptide["expected_treatment_log2fc"],
            treatment_peptide["observed_treatment_log2fc"],
        )),
        "protein_control": asdict(calculate_metrics(
            control_protein["expected_control_log2fc"],
            control_protein["observed_control_log2fc"],
        )),
        "protein_treatment": asdict(calculate_metrics(
            treatment_protein["expected_treatment_log2fc"],
            treatment_protein["observed_treatment_log2fc"],
        )),
    }
    peptide_classification, _ = classification_metrics(treatment_peptide)
    protein_classification, _ = classification_metrics(treatment_protein)
    summaries["peptide_classification"] = peptide_classification
    summaries["protein_classification"] = protein_classification
    return summaries


def _save_figure(path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


def _identity_scatter(
    table: pd.DataFrame,
    expected_column: str,
    observed_column: str,
    xlabel: str,
    ylabel: str,
    title: str,
    output: Path,
    dpi: int,
) -> None:
    expected, observed = _finite_pairs(table[expected_column], table[observed_column])
    if len(expected) == 0:
        return
    plt.figure(figsize=(6.3, 5.5))
    plt.scatter(expected, observed, s=14, alpha=0.65)
    lower = float(min(np.min(expected), np.min(observed)))
    upper = float(max(np.max(expected), np.max(observed)))
    plt.plot([lower, upper], [lower, upper], "--", linewidth=1.2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    _save_figure(output, dpi)


def create_plots(
    precursor_long: pd.DataFrame,
    run_metrics: pd.DataFrame,
    signal_recovery: pd.DataFrame,
    peptide: pd.DataFrame,
    protein: pd.DataFrame,
    config: BenchmarkConfig,
) -> list[Path]:
    figure_dir = config.output_dir / "figures"
    outputs: list[Path] = []

    path = figure_dir / "01_target_recovery_per_run.png"
    plt.figure(figsize=(7.2, 4.8))
    labels = [role.replace("_", " ") for role in run_metrics["run_role"]]
    values = run_metrics["recovery_rate"] * 100.0
    bars = plt.bar(labels, values)
    plt.ylabel("Recovered target precursors (%)")
    plt.ylim(0, 105)
    plt.title("OpenDIA recovery of simulated target precursors")
    for bar, value, count in zip(bars, values, run_metrics["detected_precursors"]):
        plt.text(bar.get_x() + bar.get_width() / 2, value + 2, f"{int(count)}\n({value:.1f}%)", ha="center")
    _save_figure(path, config.dpi)
    outputs.append(path)

    path = figure_dir / "02_detection_heatmap.png"
    detection = precursor_long.pivot(index="precursor_key", columns="run_role", values="detected")
    order = (
        precursor_long[precursor_long["run_role"] == "baseline"]
        .sort_values("evaluation_rt")["precursor_key"]
        .tolist()
    )
    detection = detection.reindex(index=order, columns=list(RUN_ROLES)).fillna(False)
    plt.figure(figsize=(5.8, 9.0))
    plt.imshow(detection.to_numpy(dtype=float), aspect="auto", interpolation="nearest", vmin=0, vmax=1)
    plt.xticks(range(len(RUN_ROLES)), [role.replace("_", "\n") for role in RUN_ROLES])
    plt.yticks([])
    plt.xlabel("Run")
    plt.ylabel("Target precursors ordered by realized RT")
    plt.title("Target detection matrix")
    plt.colorbar(label="Detected")
    _save_figure(path, config.dpi)
    outputs.append(path)

    coordinate_label = f"Realized TimSim {config.realized_coordinate}"
    for role in RUN_ROLES:
        subset = precursor_long[(precursor_long["run_role"] == role) & precursor_long["detected"]]
        path = figure_dir / f"03_rt_realized_vs_observed_{role}.png"
        _identity_scatter(
            subset, "evaluation_rt", "observed_rt",
            f"{coordinate_label} RT (s)", "Observed OpenDIA RT (s)",
            f"RT accuracy: {role.replace('_', ' ')}", path, config.dpi,
        )
        if path.exists():
            outputs.append(path)

    for role in RUN_ROLES:
        subset = precursor_long[(precursor_long["run_role"] == role) & precursor_long["detected"]]
        path = figure_dir / f"04_im_realized_vs_observed_{role}.png"
        _identity_scatter(
            subset, "evaluation_im", "observed_im",
            f"{coordinate_label} inverse mobility (1/K₀)", "Observed EXP_IM (1/K₀)",
            f"Ion-mobility accuracy: {role.replace('_', ' ')}", path, config.dpi,
        )
        if path.exists():
            outputs.append(path)

    complete_control = peptide[peptide["complete_controls"]]
    path = figure_dir / "05_control_expected_vs_observed_log2fc.png"
    _identity_scatter(
        complete_control, "expected_control_log2fc", "observed_control_log2fc",
        "Expected control-replicate log₂ ratio", "Observed control-replicate log₂ ratio",
        "Control replicate quantitative agreement", path, config.dpi,
    )
    if path.exists():
        outputs.append(path)

    complete_peptide = peptide[peptide["complete_all_runs"]]
    path = figure_dir / "06_peptide_treatment_expected_vs_observed_log2fc.png"
    _identity_scatter(
        complete_peptide, "expected_treatment_log2fc", "observed_treatment_log2fc",
        "Expected peptide treatment log₂ fold change",
        "Observed peptide treatment log₂ fold change",
        "Peptide-level treatment-effect recovery", path, config.dpi,
    )
    if path.exists():
        outputs.append(path)

    complete_protein = protein[protein["complete_all_runs"]]
    path = figure_dir / "07_protein_treatment_expected_vs_observed_log2fc.png"
    _identity_scatter(
        complete_protein, "expected_treatment_log2fc", "observed_treatment_log2fc",
        "Expected protein treatment log₂ fold change",
        "Observed protein treatment log₂ fold change",
        "Protein-level treatment-effect recovery", path, config.dpi,
    )
    if path.exists():
        outputs.append(path)

    path = figure_dir / "08_protein_observed_log2fc_by_truth_class.png"
    classes = ["down", "unchanged", "up"]
    values = [
        complete_protein.loc[complete_protein["TreatmentClass"] == treatment_class, "observed_treatment_log2fc"].dropna().to_numpy()
        for treatment_class in classes
    ]
    if any(len(value) for value in values):
        plt.figure(figsize=(6.4, 5.0))
        plt.boxplot(values, showfliers=True)
        plt.xticks(np.arange(1, len(classes) + 1), classes)
        plt.axhline(0, linestyle="--", linewidth=1.0)
        plt.xlabel("Ground-truth protein class")
        plt.ylabel("Observed treatment log₂ fold change")
        plt.title("Observed protein response by ground-truth class")
        _save_figure(path, config.dpi)
        outputs.append(path)

    path = figure_dir / "09_log2_intensity_distributions.png"
    intensity_values = [
        np.log2(precursor_long.loc[
            (precursor_long["run_role"] == role) & precursor_long["detected"], "intensity"
        ].to_numpy(dtype=float) + config.pseudocount)
        for role in RUN_ROLES
    ]
    plt.figure(figsize=(7.2, 5.0))
    plt.boxplot(intensity_values, showfliers=False)
    plt.xticks(
        np.arange(1, len(RUN_ROLES) + 1),
        [role.replace("_", "\n") for role in RUN_ROLES],
    )
    plt.ylabel(f"log₂({config.intensity_column} + {config.pseudocount:g})")
    plt.title("Detected precursor intensity distributions")
    _save_figure(path, config.dpi)
    outputs.append(path)

    path = figure_dir / "10_detection_rate_by_treatment_class.png"
    coverage = (
        precursor_long.groupby(["run_role", "TreatmentClass"], observed=True)["detected"]
        .mean().mul(100).unstack("run_role").reindex(index=["down", "unchanged", "up"])
    )
    if not coverage.empty:
        plt.figure(figsize=(7.5, 5.0))
        x = np.arange(len(coverage.index))
        width = 0.24
        for index, role in enumerate(RUN_ROLES):
            plt.bar(x + (index - 1) * width, coverage.get(role, pd.Series(index=coverage.index, dtype=float)), width, label=role.replace("_", " "))
        plt.xticks(x, coverage.index)
        plt.ylabel("Detected precursors (%)")
        plt.ylim(0, 105)
        plt.xlabel("Ground-truth treatment class")
        plt.title("Detection coverage by treatment class")
        plt.legend()
        _save_figure(path, config.dpi)
        outputs.append(path)

    if not signal_recovery.empty:
        path = figure_dir / "11_recovery_by_realized_signal_quartile.png"
        pivot = signal_recovery.pivot(
            index="realized_signal_quartile",
            columns="run_role",
            values="recovery_rate",
        ).reindex(["Q1", "Q2", "Q3", "Q4"])
        plt.figure(figsize=(7.5, 5.0))
        x = np.arange(len(pivot.index))
        width = 0.24
        for index, role in enumerate(RUN_ROLES):
            plt.bar(
                x + (index - 1) * width,
                100.0 * pivot.get(role, pd.Series(index=pivot.index, dtype=float)),
                width,
                label=role.replace("_", " "),
            )
        plt.xticks(x, pivot.index)
        plt.ylim(0, 105)
        plt.xlabel("Realized precursor-signal quartile")
        plt.ylabel("Recovered target precursors (%)")
        plt.title("OpenDIA recovery versus realized TimSim precursor signal")
        plt.legend()
        _save_figure(path, config.dpi)
        outputs.append(path)

        baseline = precursor_long[
            (precursor_long["run_role"] == "baseline")
            & precursor_long["realized_event_proxy"].notna()
        ].copy()
        if not baseline.empty:
            path = figure_dir / "12_realized_signal_vs_opendia_intensity_baseline.png"
            plt.figure(figsize=(7.0, 5.3))
            for detected, group in baseline.groupby("detected"):
                plt.scatter(
                    np.log10(group["realized_event_proxy"].clip(lower=1.0)),
                    np.log10(group["intensity"].fillna(0).clip(lower=1.0)),
                    s=18,
                    alpha=0.65,
                    label=f"detected={bool(detected)}",
                )
            plt.xlabel("log₁₀ realized precursor event proxy")
            plt.ylabel(f"log₁₀ OpenDIA {config.intensity_column}")
            plt.title("Realized TimSim signal versus OpenDIA intensity: baseline")
            plt.legend()
            _save_figure(path, config.dpi)
            outputs.append(path)

    if "realized_rt_apex" in precursor_long.columns:
        baseline_truth = precursor_long[precursor_long["run_role"] == "baseline"].copy()
        path = figure_dir / "13_assay_rt_vs_realized_rt_apex_baseline.png"
        _identity_scatter(
            baseline_truth,
            "assay_rt",
            "realized_rt_apex",
            "Assay/predictor RT (s)",
            "Realized TimSim RT apex (s)",
            "Assay RT prior versus realized chromatographic apex",
            path,
            config.dpi,
        )
        if path.exists():
            outputs.append(path)

    return outputs



def _format_float(value: Any, digits: int = 4) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    return "NA" if not np.isfinite(numeric) else f"{numeric:.{digits}f}"


def _markdown_table(frame: pd.DataFrame, columns: Sequence[str]) -> str:
    subset = frame[list(columns)].copy()
    header = "| " + " | ".join(columns) + " |"
    divider = "|" + "|".join("---" for _ in columns) + "|"
    rows = []
    for _, row in subset.iterrows():
        rows.append("| " + " | ".join(str(row[column]) for column in columns) + " |")
    return "\n".join([header, divider, *rows])


def write_report(
    config: BenchmarkConfig,
    run_metrics: pd.DataFrame,
    signal_recovery: pd.DataFrame,
    quant_summary: dict[str, Any],
    peptide: pd.DataFrame,
    protein: pd.DataFrame,
    result_counts: dict[str, int],
    q_column: str | None,
    plot_paths: list[Path],
) -> Path:
    report_path = config.output_dir / "benchmark_report.md"
    run_table = run_metrics.copy()
    run_table["recovery"] = run_table.apply(
        lambda row: f"{int(row.detected_precursors)}/{int(row.target_precursors)} ({100 * row.recovery_rate:.1f}%)",
        axis=1,
    )
    run_table["RT MAE (s)"] = run_table["rt_mae"].map(_format_float)
    run_table["RT R²"] = run_table["rt_r_squared"].map(_format_float)
    run_table["IM MAE"] = run_table["im_mae"].map(lambda value: _format_float(value, 6))
    run_table["IM R²"] = run_table["im_r_squared"].map(_format_float)

    peptide_control = quant_summary["peptide_control"]
    peptide_treatment = quant_summary["peptide_treatment"]
    protein_treatment = quant_summary["protein_treatment"]
    peptide_classification = quant_summary["peptide_classification"]
    protein_classification = quant_summary["protein_classification"]

    realized_label = (
        f"realized TimSim signal {config.realized_coordinate}"
        if config.use_realized_truth
        else "assay/predictor coordinates"
    )

    lines = [
        "# TimSim/OpenDIA synthetic benchmark",
        "",
        "## Scope",
        "",
        f"The transition library contains **{int(run_metrics.target_precursors.iloc[0])} simulated target precursor groups**. "
        "All selected peptides are present in the simulation, but their realized precursor signal spans a broad dynamic range. "
        "Peptide-level TimSim event counts therefore must not be interpreted as charge-state-specific detectability.",
        "",
        "The assay library contains predicted RT/IM coordinates. Coordinate accuracy in this report is evaluated against "
        f"the **{realized_label}**. Peptide-level event counts remain the quantitative truth for between-run fold changes.",
        "",
        "This three-run design has two controls and one treatment. It can validate fold-change recovery, but it cannot validate "
        "differential-expression p-values, empirical power, or treatment-group variance. Those require multiple independently "
        "simulated treatment replicates.",
        "",
        "## Identification and coordinate accuracy",
        "",
        _markdown_table(run_table, ["run_role", "recovery", "RT MAE (s)", "RT R²", "IM MAE", "IM R²"]),
        "",
        f"OpenDIA result rows before benchmark filtering: **{result_counts['input_rows']}**. "
        f"The benchmark selected one best target feature per precursor/run using **{q_column or 'no q-value column'}**"
        + (f" at q ≤ {config.q_threshold:g}." if q_column else "."),
        "",
    ]

    if not signal_recovery.empty:
        recovery_table = signal_recovery.copy()
        recovery_table["recovery"] = recovery_table.apply(
            lambda row: f"{int(row.detected)}/{int(row.targets)} ({100 * row.recovery_rate:.1f}%)",
            axis=1,
        )
        recovery_table["median realized signal"] = recovery_table[
            "median_realized_event_proxy"
        ].map(lambda value: _format_float(value, 1))
        lines.extend(
            [
                "## Recovery versus realized precursor signal",
                "",
                "The realized precursor event proxy is "
                "`peptide.events × sum(frame_abundance) × ion.relative_abundance × sum(scan_abundance)`. "
                "It is a relative detectability proxy, not a claim about TimSim's internal absolute ion-count arithmetic.",
                "",
                _markdown_table(
                    recovery_table,
                    [
                        "run_role",
                        "realized_signal_quartile",
                        "recovery",
                        "median realized signal",
                    ],
                ),
                "",
            ]
        )

    lines.extend(
        [
            "## Control-replicate quantitative agreement",
            "",
            f"Complete precursor pairs: **{peptide_control['n']}**. Expected-versus-observed log₂-ratio "
            f"Pearson r = **{_format_float(peptide_control['pearson_r'])}**, "
            f"MAE = **{_format_float(peptide_control['mae'])}**, "
            f"bias = **{_format_float(peptide_control['bias'])}**.",
            "",
            "## Treatment-effect recovery",
            "",
            f"Peptides quantified in all three runs: **{int(peptide.complete_all_runs.sum())}**. "
            f"Expected-versus-observed peptide log₂ fold-change Pearson r = **{_format_float(peptide_treatment['pearson_r'])}**, "
            f"MAE = **{_format_float(peptide_treatment['mae'])}**.",
            "",
            f"Proteins quantified in all three runs: **{int(protein.complete_all_runs.sum())}** of "
            f"**{len(protein)}** ground-truth proteins. Expected-versus-observed protein log₂ fold-change "
            f"Pearson r = **{_format_float(protein_treatment['pearson_r'])}**, "
            f"MAE = **{_format_float(protein_treatment['mae'])}**.",
            "",
            f"Using an absolute log₂ fold-change threshold of **{config.fold_change_threshold:g}**, peptide three-class "
            f"accuracy = **{_format_float(peptide_classification['three_class_accuracy'])}** and changed-protein direction "
            f"accuracy = **{_format_float(protein_classification['changed_direction_accuracy'])}**. Protein altered-versus-unchanged "
            f"AUC = **{_format_float(protein_classification['altered_vs_unchanged_auc'])}**.",
            "",
            "## Interpretation cautions",
            "",
            "- `BaselineEvents`, `ControlReplicateEvents`, and `TreatmentEvents` are peptide-level simulation inputs. "
            "They are appropriate for expected fold-change calculations, but not as precursor-level detectability truth.",
            "- The realized precursor event proxy incorporates retained chromatographic frame mass, selected ion charge-state abundance, "
            "and retained mobility scan mass. It should be interpreted comparatively rather than as an exact absolute ion count.",
            "- Assay RT/IM values are predictor coordinates supplied to OpenDIA. Realized RT/IM apex/centroid values are computed from "
            "TimSim's stored frame/scan abundance distributions and are used for coordinate-error evaluation.",
            "- Missing values are treated as failed recovery, not as zero abundance, in coordinate and fold-change comparisons.",
            "- The final `OpenDIA.results.tsv` commonly excludes decoys; inspect the internal scored tables if decoy-score calibration is required.",
            "",
            "## Output files",
            "",
            "- `run_metrics.tsv`: per-run recovery plus realized-coordinate and assay-prior RT/IM metrics.",
            "- `realized_signal_recovery.tsv`: recovery by realized precursor-signal quartile.",
            "- `precursor_benchmark.tsv`: every ground-truth precursor/run with assay priors, realized truth, OpenDIA results, and errors.",
            "- `peptide_quant_benchmark.tsv`: control and treatment fold-change recovery at precursor/peptide level.",
            "- `protein_quant_benchmark.tsv`: protein-level effect recovery.",
            "- `peptide_confusion_matrix.tsv` and `protein_confusion_matrix.tsv`: up/down/unchanged classification.",
            "- `benchmark_summary.json`: machine-readable summary.",
            "",
            "## Figures",
            "",
        ]
    )
    for plot in plot_paths:
        relative = plot.relative_to(config.output_dir)
        lines.append(f"- `{relative}`")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path



def run_benchmark(config: BenchmarkConfig) -> dict[str, Any]:
    if config.pseudocount <= 0:
        raise ValueError("pseudocount must be positive")
    if config.fold_change_threshold < 0:
        raise ValueError("fold-change threshold must be non-negative")
    config.output_dir.mkdir(parents=True, exist_ok=True)
    paths = _resolve_paths(config)

    results = _read_tsv(paths["results"])
    precursor_truth = _read_tsv(paths["precursor_truth"])
    condition_truth = _read_tsv(paths["condition_truth"])
    protein_matrix = _read_tsv(paths["protein_matrix"])
    # Read the peptide matrix for input validation and provenance. Quantitative
    # benchmarking uses the richer results TSV so charge and coordinates remain available.
    peptide_matrix = _read_tsv(paths["peptide_matrix"])

    run_mapping = _run_mapping(config)
    for matrix_name, matrix in (("peptide", peptide_matrix), ("protein", protein_matrix)):
        missing_runs = [run for run in run_mapping.values() if run not in matrix.columns]
        if missing_runs:
            raise ValueError(f"OpenDIA {matrix_name} matrix is missing run columns: {missing_runs}")

    best_results, q_column, result_counts = prepare_results(results, config)
    truth = prepare_truth(precursor_truth, condition_truth)
    realized_truth = build_realized_truth(truth, config, run_mapping)
    precursor_long = build_precursor_long(
        truth, best_results, run_mapping, realized_truth, config
    )
    run_metrics = _run_metrics(precursor_long, q_column)
    signal_recovery = build_signal_recovery_summary(precursor_long)
    peptide_quant = build_peptide_quant_benchmark(precursor_long, config)
    protein_quant = build_protein_quant_benchmark(
        protein_matrix, condition_truth, config, run_mapping
    )
    quant_summary = _quant_summary(peptide_quant, protein_quant)

    peptide_classification, peptide_confusion = classification_metrics(
        peptide_quant[peptide_quant["complete_all_runs"]]
    )
    protein_classification, protein_confusion = classification_metrics(
        protein_quant[protein_quant["complete_all_runs"]]
    )
    quant_summary["peptide_classification"] = peptide_classification
    quant_summary["protein_classification"] = protein_classification

    run_metrics_path = config.output_dir / "run_metrics.tsv"
    precursor_path = config.output_dir / "precursor_benchmark.tsv"
    peptide_path = config.output_dir / "peptide_quant_benchmark.tsv"
    protein_path = config.output_dir / "protein_quant_benchmark.tsv"
    peptide_confusion_path = config.output_dir / "peptide_confusion_matrix.tsv"
    protein_confusion_path = config.output_dir / "protein_confusion_matrix.tsv"
    signal_recovery_path = config.output_dir / "realized_signal_recovery.tsv"
    run_metrics.to_csv(run_metrics_path, sep="\t", index=False)
    precursor_long.to_csv(precursor_path, sep="\t", index=False)
    peptide_quant.to_csv(peptide_path, sep="\t", index=False)
    protein_quant.to_csv(protein_path, sep="\t", index=False)
    peptide_confusion.to_csv(peptide_confusion_path, sep="\t", index=False)
    protein_confusion.to_csv(protein_confusion_path, sep="\t", index=False)
    signal_recovery.to_csv(signal_recovery_path, sep="\t", index=False)

    plot_paths = create_plots(
        precursor_long, run_metrics, signal_recovery, peptide_quant, protein_quant, config
    )

    summary: dict[str, Any] = {
        "configuration": {
            **asdict(config),
            **{key: str(value) for key, value in asdict(config).items() if isinstance(value, Path)},
        },
        "input_paths": {key: str(value) for key, value in paths.items()},
        "q_column": q_column,
        "result_counts": result_counts,
        "run_metrics": run_metrics.to_dict(orient="records"),
        "realized_signal_recovery": signal_recovery.to_dict(orient="records"),
        "realized_truth": {
            "enabled": config.use_realized_truth,
            "coordinate": config.realized_coordinate,
            "synthetic_db_root": str(config.synthetic_db_root or config.build_dir),
        },
        "quantitative_metrics": quant_summary,
        "ground_truth": {
            "selected_precursors": int(len(truth)),
            "selected_peptides": int(truth["sequence"].nunique()),
            "selected_proteins": int(truth["ProteinId"].nunique()),
            "all_condition_truth_proteins": int(condition_truth["ProteinId"].nunique()),
            "treatment_class_counts_all_proteins": condition_truth.drop_duplicates("ProteinId")["TreatmentClass"].value_counts().to_dict(),
        },
        "outputs": {
            "run_metrics": str(run_metrics_path),
            "precursor_benchmark": str(precursor_path),
            "peptide_quant_benchmark": str(peptide_path),
            "protein_quant_benchmark": str(protein_path),
            "peptide_confusion_matrix": str(peptide_confusion_path),
            "protein_confusion_matrix": str(protein_confusion_path),
            "realized_signal_recovery": str(signal_recovery_path),
            "figures": [str(path) for path in plot_paths],
        },
    }
    summary_path = config.output_dir / "benchmark_summary.json"
    report_path = write_report(
        config, run_metrics, signal_recovery, quant_summary, peptide_quant, protein_quant,
        result_counts, q_column, plot_paths,
    )
    summary["outputs"]["summary"] = str(summary_path)
    summary["outputs"]["report"] = str(report_path)
    summary_path.write_text(
        json.dumps(summary, indent=2, allow_nan=True, default=str) + "\n",
        encoding="utf-8",
    )

    print(f"Benchmark complete: {config.output_dir}")
    print(run_metrics[["run_role", "detected_precursors", "target_precursors", "recovery_rate", "rt_mae", "im_mae"]].to_string(index=False))
    print(
        "Protein treatment effect: "
        f"n={quant_summary['protein_treatment']['n']}, "
        f"Pearson r={_format_float(quant_summary['protein_treatment']['pearson_r'])}, "
        f"MAE={_format_float(quant_summary['protein_treatment']['mae'])}"
    )
    print(f"Report: {report_path}")
    return summary


def _path_or_none(value: str | None) -> Path | None:
    return Path(value).expanduser().resolve() if value else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--opendia-dir", type=Path, default=Path("opendia_test"))
    parser.add_argument("--out-dir", type=Path, default=Path("benchmark_results"))
    parser.add_argument("--results-tsv")
    parser.add_argument("--peptide-matrix-tsv")
    parser.add_argument("--protein-matrix-tsv")
    parser.add_argument("--precursor-truth-tsv")
    parser.add_argument("--condition-truth-tsv")
    parser.add_argument(
        "--synthetic-db-root",
        type=Path,
        help=(
            "Directory containing <run_name>/synthetic_data.db. "
            "Defaults to --build-dir."
        ),
    )
    parser.add_argument(
        "--realized-coordinate",
        choices=("apex", "centroid"),
        default="apex",
        help="Realized TimSim coordinate used for RT/IM accuracy metrics.",
    )
    parser.add_argument(
        "--skip-realized-truth",
        action="store_true",
        help="Use assay/predictor RT/IM instead of reading realized TimSim frame/scan truth.",
    )
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument("--q-column", help="Column used for filtering; use 'none' to disable")
    parser.add_argument("--intensity-column", default="Intensity")
    parser.add_argument("--pseudocount", type=float, default=1.0)
    parser.add_argument("--fold-change-threshold", type=float, default=0.5)
    parser.add_argument("--baseline-run", default=DEFAULT_RUN_NAMES["baseline"])
    parser.add_argument("--control-replicate-run", default=DEFAULT_RUN_NAMES["control_replicate"])
    parser.add_argument("--treatment-run", default=DEFAULT_RUN_NAMES["treatment"])
    parser.add_argument("--dpi", type=int, default=170)
    args = parser.parse_args()

    config = BenchmarkConfig(
        build_dir=args.build_dir.expanduser().resolve(),
        opendia_dir=args.opendia_dir.expanduser().resolve(),
        output_dir=args.out_dir.expanduser().resolve(),
        results_tsv=_path_or_none(args.results_tsv),
        peptide_matrix_tsv=_path_or_none(args.peptide_matrix_tsv),
        protein_matrix_tsv=_path_or_none(args.protein_matrix_tsv),
        precursor_truth_tsv=_path_or_none(args.precursor_truth_tsv),
        condition_truth_tsv=_path_or_none(args.condition_truth_tsv),
        synthetic_db_root=args.synthetic_db_root.expanduser().resolve()
        if args.synthetic_db_root
        else None,
        realized_coordinate=args.realized_coordinate,
        use_realized_truth=not args.skip_realized_truth,
        q_threshold=args.q_threshold,
        q_column=args.q_column,
        intensity_column=args.intensity_column,
        pseudocount=args.pseudocount,
        fold_change_threshold=args.fold_change_threshold,
        baseline_run=args.baseline_run,
        control_replicate_run=args.control_replicate_run,
        treatment_run=args.treatment_run,
        dpi=args.dpi,
    )
    run_benchmark(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
