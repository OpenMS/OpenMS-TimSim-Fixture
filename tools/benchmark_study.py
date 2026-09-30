#!/usr/bin/env python3
"""Benchmark a multi-run control/treatment OpenDIA study against TimSim truth."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class MetricResult:
    n: int
    bias: float
    mae: float
    rmse: float
    pearson_r: float
    spearman_rho: float


def first_existing(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    available = set(columns)
    return next((value for value in candidates if value in available), None)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def normalise_decoy(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.lower()
    return text.isin({"1", "true", "yes", "decoy"}) | numeric(series).fillna(0).ne(0)


def metric(expected: Iterable[float], observed: Iterable[float]) -> MetricResult:
    x = pd.to_numeric(pd.Series(expected), errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(observed), errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) == 0:
        return MetricResult(0, *(math.nan for _ in range(5)))
    error = y - x
    pearson = float(stats.pearsonr(x, y).statistic) if len(x) > 1 and not np.allclose(x, x[0]) and not np.allclose(y, y[0]) else math.nan
    spearman = float(stats.spearmanr(x, y).statistic) if len(x) > 1 and not np.allclose(x, x[0]) and not np.allclose(y, y[0]) else math.nan
    return MetricResult(
        n=int(len(x)),
        bias=float(np.mean(error)),
        mae=float(np.mean(np.abs(error))),
        rmse=float(np.sqrt(np.mean(error ** 2))),
        pearson_r=pearson,
        spearman_rho=spearman,
    )


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(values), np.nan, dtype=float)
    finite = np.isfinite(values)
    if not finite.any():
        return pd.Series(result, index=p_values.index)
    idx = np.flatnonzero(finite)
    order = idx[np.argsort(values[idx])]
    ranked = values[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    result[order] = adjusted
    return pd.Series(result, index=p_values.index)


def safe_log2_ratio(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or numerator <= 0 or denominator <= 0:
        return math.nan
    return float(math.log2(numerator / denominator))


def welch(control: pd.Series, treatment: pd.Series) -> float:
    a = pd.to_numeric(control, errors="coerce").dropna().to_numpy(dtype=float)
    b = pd.to_numeric(treatment, errors="coerce").dropna().to_numpy(dtype=float)
    if len(a) < 2 or len(b) < 2:
        return math.nan
    return float(stats.ttest_ind(b, a, equal_var=False).pvalue)


def prepare_results(results: pd.DataFrame, run_names: set[str], q_threshold: float) -> tuple[pd.DataFrame, str | None]:
    required = {"run_name", "Sequence", "Charge", "RT", "EXP_IM", "Intensity"}
    missing = required - set(results.columns)
    if missing:
        raise SystemExit(f"OpenDIA results are missing required columns: {sorted(missing)}")
    work = results.copy()
    work["run_name"] = work["run_name"].astype(str)
    work = work[work["run_name"].isin(run_names)].copy()
    work["sequence"] = work["Sequence"].astype(str).str.strip().str.upper()
    work["charge"] = numeric(work["Charge"]).astype("Int64")
    work["intensity"] = numeric(work["Intensity"])
    work["observed_rt"] = numeric(work["RT"])
    work["observed_im"] = numeric(work["EXP_IM"])
    work["is_decoy"] = normalise_decoy(work["decoy"]) if "decoy" in work else False
    work = work[~work["is_decoy"]].copy()
    q_column = first_existing(
        work.columns,
        ("m_score", "ms2_m_score", "m_score_peptide_run_specific", "m_score_peptide_experiment_wide", "m_score_peptide_global"),
    )
    if q_column:
        work["q_value"] = numeric(work[q_column])
        work = work[work["q_value"].isna() | work["q_value"].le(q_threshold)].copy()
    else:
        work["q_value"] = math.nan
    work["d_score_numeric"] = numeric(work["d_score"]) if "d_score" in work else math.nan
    work = work.sort_values(
        ["run_name", "sequence", "charge", "q_value", "d_score_numeric"],
        ascending=[True, True, True, True, False],
        na_position="last",
    ).drop_duplicates(["run_name", "sequence", "charge"], keep="first")
    return work, q_column


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--opendia-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--q-threshold", type=float, default=0.01)
    parser.add_argument("--differential-q-threshold", type=float, default=0.05)
    parser.add_argument("--effect-threshold", type=float, default=0.5)
    args = parser.parse_args()

    build = args.build_dir.resolve()
    opendia = args.opendia_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    study_manifest = pd.read_csv(build / "OpenSwathTimSim.study_manifest.tsv", sep="\t")
    selected = pd.read_csv(build / "OpenSwathTimSim.high_signal_selection.tsv", sep="\t")
    realized = pd.read_csv(build / "OpenSwathTimSim.realized_truth.tsv", sep="\t")
    design = pd.read_csv(build / "OpenSwathTimSim.study_design_truth.tsv", sep="\t")
    results = pd.read_csv(opendia / "OpenDIA.results.tsv", sep="\t", low_memory=False)

    expected_runs = set(study_manifest["RunName"].astype(str))
    target_keys = {
        (str(row.sequence).upper(), int(row.charge))
        for row in selected.itertuples(index=False)
    }
    prepared, q_column = prepare_results(results, expected_runs, args.q_threshold)
    prepared = prepared[
        prepared.apply(lambda row: (str(row["sequence"]), int(row["charge"])) in target_keys, axis=1)
    ].copy()

    truth = realized.copy()
    truth["sequence"] = truth["PeptideSequence"].astype(str).str.upper()
    truth["charge"] = numeric(truth["PrecursorCharge"]).astype(int)
    truth["run_name"] = truth["RunName"].astype(str)
    truth["realized_rt"] = numeric(truth["RealizedRTApex"])
    truth["realized_im"] = numeric(truth["RealizedIMApexSQLite"])
    truth["realized_input_events"] = numeric(truth["RealizedInputEvents"])
    truth["design_log2fc"] = numeric(truth["DesignLog2FC"])
    truth["condition"] = truth["Condition"].astype(str)
    truth["protein_id"] = truth["ProteinId"].astype(str)
    truth["treatment_class"] = truth["TreatmentClass"].astype(str)

    measurement = truth.merge(
        prepared[["run_name", "sequence", "charge", "intensity", "observed_rt", "observed_im", "q_value"]],
        on=["run_name", "sequence", "charge"],
        how="left",
        validate="one_to_one",
    )
    measurement["detected"] = measurement["intensity"].notna()
    measurement["rt_error_seconds"] = measurement["observed_rt"] - measurement["realized_rt"]
    measurement["im_error"] = measurement["observed_im"] - measurement["realized_im"]
    measurement.to_csv(out / "precursor_measurements.tsv", sep="\t", index=False)

    run_rows = []
    for row in study_manifest.sort_values("RunOrdinal").itertuples(index=False):
        subset = measurement[measurement["run_name"] == row.RunName]
        detected = subset[subset["detected"]]
        run_rows.append({
            "RunName": row.RunName,
            "Condition": row.Condition,
            "Replicate": int(row.Replicate),
            "targets": len(subset),
            "identified": int(subset["detected"].sum()),
            "recovery": float(subset["detected"].mean()),
            "rt_mae_seconds": float(detected["rt_error_seconds"].abs().mean()) if len(detected) else math.nan,
            "im_mae": float(detected["im_error"].abs().mean()) if len(detected) else math.nan,
        })
    run_metrics = pd.DataFrame(run_rows)
    run_metrics.to_csv(out / "run_metrics.tsv", sep="\t", index=False)

    peptide_rows = []
    for (sequence, charge), group in measurement.groupby(["sequence", "charge"], sort=True):
        control = group[group["condition"] == "control"]
        treatment = group[group["condition"] == "treatment"]
        control_obs = numeric(control["intensity"]).dropna()
        treatment_obs = numeric(treatment["intensity"]).dropna()
        control_truth = numeric(control["realized_input_events"]).dropna()
        treatment_truth = numeric(treatment["realized_input_events"]).dropna()
        design_log2fc = float(group["design_log2fc"].iloc[0])
        observed_log2fc = safe_log2_ratio(float(treatment_obs.mean()), float(control_obs.mean())) if len(control_obs) and len(treatment_obs) else math.nan
        realized_log2fc = safe_log2_ratio(float(treatment_truth.mean()), float(control_truth.mean()))
        p_value = welch(np.log2(control_obs[control_obs > 0]), np.log2(treatment_obs[treatment_obs > 0]))
        peptide_rows.append({
            "PeptideSequence": sequence,
            "PrecursorCharge": int(charge),
            "ProteinId": group["protein_id"].iloc[0],
            "TreatmentClass": group["treatment_class"].iloc[0],
            "DesignLog2FC": design_log2fc,
            "RealizedConditionLog2FC": realized_log2fc,
            "ObservedLog2FC": observed_log2fc,
            "ControlDetected": int(control["detected"].sum()),
            "TreatmentDetected": int(treatment["detected"].sum()),
            "ControlRuns": int(len(control)),
            "TreatmentRuns": int(len(treatment)),
            "WelchPValue": p_value,
        })
    peptide = pd.DataFrame(peptide_rows)
    peptide["BH_QValue"] = bh_adjust(peptide["WelchPValue"])
    peptide["CalledDifferential"] = peptide["BH_QValue"].le(args.differential_q_threshold) & peptide["ObservedLog2FC"].abs().ge(args.effect_threshold)
    peptide["TruthDifferential"] = peptide["TreatmentClass"].ne("unchanged")
    peptide.to_csv(out / "peptide_effects.tsv", sep="\t", index=False)

    # Aggregate selected precursor signal to protein per run, then perform the same two-group comparison.
    protein_run = (
        measurement.groupby(["run_name", "condition", "protein_id"], as_index=False)
        .agg(
            observed_intensity=("intensity", "sum"),
            realized_events=("realized_input_events", "sum"),
            detected_precursors=("detected", "sum"),
            selected_precursors=("detected", "size"),
        )
    )
    design_map = design.set_index("ProteinId")
    protein_rows = []
    for protein_id, group in protein_run.groupby("protein_id", sort=True):
        control = group[group["condition"] == "control"]
        treatment = group[group["condition"] == "treatment"]
        obs_control = numeric(control["observed_intensity"]).replace(0, np.nan).dropna()
        obs_treatment = numeric(treatment["observed_intensity"]).replace(0, np.nan).dropna()
        truth_control = numeric(control["realized_events"]).replace(0, np.nan).dropna()
        truth_treatment = numeric(treatment["realized_events"]).replace(0, np.nan).dropna()
        if protein_id in design_map.index:
            drow = design_map.loc[protein_id]
            design_log2fc = float(drow["DesignLog2FC"])
            treatment_class = str(drow["TreatmentClass"])
        else:
            design_log2fc = math.nan
            treatment_class = "unknown"
        protein_rows.append({
            "ProteinId": protein_id,
            "TreatmentClass": treatment_class,
            "DesignLog2FC": design_log2fc,
            "RealizedConditionLog2FC": safe_log2_ratio(float(truth_treatment.mean()), float(truth_control.mean())) if len(truth_control) and len(truth_treatment) else math.nan,
            "ObservedLog2FC": safe_log2_ratio(float(obs_treatment.mean()), float(obs_control.mean())) if len(obs_control) and len(obs_treatment) else math.nan,
            "WelchPValue": welch(np.log2(obs_control), np.log2(obs_treatment)),
            "ControlRuns": int(len(control)),
            "TreatmentRuns": int(len(treatment)),
            "SelectedPrecursors": int(group["selected_precursors"].max()),
        })
    protein = pd.DataFrame(protein_rows)
    protein["BH_QValue"] = bh_adjust(protein["WelchPValue"])
    protein["CalledDifferential"] = protein["BH_QValue"].le(args.differential_q_threshold) & protein["ObservedLog2FC"].abs().ge(args.effect_threshold)
    protein["TruthDifferential"] = protein["TreatmentClass"].ne("unchanged")
    protein.to_csv(out / "protein_effects.tsv", sep="\t", index=False)

    peptide_design_realized = metric(peptide["DesignLog2FC"], peptide["RealizedConditionLog2FC"])
    peptide_realized_observed = metric(peptide["RealizedConditionLog2FC"], peptide["ObservedLog2FC"])
    protein_design_realized = metric(protein["DesignLog2FC"], protein["RealizedConditionLog2FC"])
    protein_realized_observed = metric(protein["RealizedConditionLog2FC"], protein["ObservedLog2FC"])

    def differential_summary(table: pd.DataFrame) -> dict[str, float | int]:
        truth = table["TruthDifferential"].fillna(False)
        call = table["CalledDifferential"].fillna(False)
        tp = int((truth & call).sum())
        fp = int((~truth & call).sum())
        fn = int((truth & ~call).sum())
        tn = int((~truth & ~call).sum())
        return {
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "sensitivity": tp / (tp + fn) if tp + fn else math.nan,
            "false_positive_rate": fp / (fp + tn) if fp + tn else math.nan,
            "precision": tp / (tp + fp) if tp + fp else math.nan,
        }

    summary = {
        "q_column": q_column,
        "identification_q_threshold": args.q_threshold,
        "differential_q_threshold": args.differential_q_threshold,
        "effect_threshold_abs_log2fc": args.effect_threshold,
        "runs": {
            "total": int(len(study_manifest)),
            "control": int((study_manifest["Condition"] == "control").sum()),
            "treatment": int((study_manifest["Condition"] == "treatment").sum()),
        },
        "targets": int(len(selected)),
        "overall_recovery": float(measurement["detected"].mean()),
        "peptide_design_vs_realized": asdict(peptide_design_realized),
        "peptide_realized_vs_observed": asdict(peptide_realized_observed),
        "protein_design_vs_realized": asdict(protein_design_realized),
        "protein_realized_vs_observed": asdict(protein_realized_observed),
        "peptide_differential": differential_summary(peptide),
        "protein_differential": differential_summary(protein),
    }
    (out / "study_benchmark_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")

    report = [
        "# OpenDIA multi-run TimSim study benchmark",
        "",
        f"Study: **{summary['runs']['control']} control + {summary['runs']['treatment']} treatment runs**; "
        f"**{summary['targets']} frozen target precursors**.",
        "",
        f"Overall target/run recovery at q <= {args.q_threshold:g}: **{100.0 * summary['overall_recovery']:.2f}%**.",
        "",
        "## Quantitative truth layers",
        "",
        "`DesignLog2FC` is the predeclared treatment effect. `RealizedConditionLog2FC` is computed from the "
        "actual TimSim input abundance realized across replicate source databases. `ObservedLog2FC` is derived from OpenDIA intensity.",
        "",
        "| Level | Comparison | n | Pearson r | MAE | Bias |",
        "|---|---|---:|---:|---:|---:|",
        f"| Peptide | design vs realized | {peptide_design_realized.n} | {peptide_design_realized.pearson_r:.4f} | {peptide_design_realized.mae:.4f} | {peptide_design_realized.bias:.4f} |",
        f"| Peptide | realized vs OpenDIA | {peptide_realized_observed.n} | {peptide_realized_observed.pearson_r:.4f} | {peptide_realized_observed.mae:.4f} | {peptide_realized_observed.bias:.4f} |",
        f"| Protein | design vs realized | {protein_design_realized.n} | {protein_design_realized.pearson_r:.4f} | {protein_design_realized.mae:.4f} | {protein_design_realized.bias:.4f} |",
        f"| Protein | realized vs OpenDIA | {protein_realized_observed.n} | {protein_realized_observed.pearson_r:.4f} | {protein_realized_observed.mae:.4f} | {protein_realized_observed.bias:.4f} |",
        "",
        "## Differential abundance",
        "",
        f"Calls use BH q <= {args.differential_q_threshold:g} and |observed log2FC| >= {args.effect_threshold:g}.",
        "",
        f"Peptide: `{summary['peptide_differential']}`",
        "",
        f"Protein: `{summary['protein_differential']}`",
        "",
        "## Outputs",
        "",
        "- `run_metrics.tsv`: per-run recovery and RT/IM accuracy.",
        "- `precursor_measurements.tsv`: target/run truth joined to OpenDIA measurements.",
        "- `peptide_effects.tsv`: design, realized, observed effects and Welch/BH statistics.",
        "- `protein_effects.tsv`: selected-protein aggregate effects and statistics.",
        "- `study_benchmark_summary.json`: machine-readable summary.",
        "",
    ]
    (out / "study_benchmark_report.md").write_text("\n".join(report), encoding="utf-8")

    print(f"Study benchmark complete: {out}")
    print(f"Overall target/run recovery: {100.0 * summary['overall_recovery']:.2f}%")
    print(
        "Peptide realized-vs-observed effect: "
        f"r={peptide_realized_observed.pearson_r:.4f}, MAE={peptide_realized_observed.mae:.4f}"
    )
    print(
        "Protein realized-vs-observed effect: "
        f"r={protein_realized_observed.pearson_r:.4f}, MAE={protein_realized_observed.mae:.4f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
