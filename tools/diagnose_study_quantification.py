#!/usr/bin/env python3
"""Diagnose multi-run OpenDIA quantification errors against realized TimSim truth."""
from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


@dataclass(frozen=True)
class EffectMetric:
    n: int
    bias: float
    mae: float
    rmse: float
    pearson_r: float
    spearman_rho: float
    slope: float
    intercept: float


def finite_xy(x: pd.Series, y: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    xa = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    ya = pd.to_numeric(y, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(xa) & np.isfinite(ya)
    return xa[mask], ya[mask]


def effect_metric(expected: pd.Series, observed: pd.Series) -> EffectMetric:
    x, y = finite_xy(expected, observed)
    if len(x) == 0:
        return EffectMetric(0, *(math.nan for _ in range(7)))
    error = y - x
    if len(x) > 1 and not np.allclose(x, x[0]) and not np.allclose(y, y[0]):
        pearson = float(stats.pearsonr(x, y).statistic)
        spearman = float(stats.spearmanr(x, y).statistic)
        regression = stats.linregress(x, y)
        slope = float(regression.slope)
        intercept = float(regression.intercept)
    else:
        pearson = spearman = slope = intercept = math.nan
    return EffectMetric(
        n=int(len(x)),
        bias=float(np.mean(error)),
        mae=float(np.mean(np.abs(error))),
        rmse=float(np.sqrt(np.mean(error**2))),
        pearson_r=pearson,
        spearman_rho=spearman,
        slope=slope,
        intercept=intercept,
    )


def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    xa, ya = finite_xy(x, y)
    if len(xa) < 3 or np.allclose(xa, xa[0]) or np.allclose(ya, ya[0]):
        return math.nan
    return float(stats.spearmanr(xa, ya).statistic)


def trimmed_metric(table: pd.DataFrame, expected: str, observed: str, fraction: float = 0.01) -> EffectMetric:
    work = table[[expected, observed]].copy()
    work["residual"] = pd.to_numeric(work[observed], errors="coerce") - pd.to_numeric(work[expected], errors="coerce")
    work = work[np.isfinite(work["residual"])].copy()
    if len(work) < 2:
        return effect_metric(work[expected], work[observed])
    keep = max(1, int(math.floor(len(work) * (1.0 - fraction))))
    work = work.assign(abs_residual=work["residual"].abs()).nsmallest(keep, "abs_residual")
    return effect_metric(work[expected], work[observed])


def residual_summary(residual: pd.Series) -> dict[str, float | int]:
    values = pd.to_numeric(residual, errors="coerce").dropna().to_numpy(dtype=float)
    if len(values) == 0:
        return {"n": 0}
    abs_values = np.abs(values)
    return {
        "n": int(len(values)),
        "median_abs": float(np.median(abs_values)),
        "p90_abs": float(np.quantile(abs_values, 0.90)),
        "p95_abs": float(np.quantile(abs_values, 0.95)),
        "p99_abs": float(np.quantile(abs_values, 0.99)),
        "gt_0_25": int(np.sum(abs_values > 0.25)),
        "gt_0_5": int(np.sum(abs_values > 0.5)),
        "gt_1_0": int(np.sum(abs_values > 1.0)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-results-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--top", type=int, default=50)
    args = parser.parse_args()

    study = args.study_results_dir.resolve()
    out = (args.out_dir or (study / "quantification_diagnostics")).resolve()
    out.mkdir(parents=True, exist_ok=True)

    measurement_path = study / "precursor_measurements.tsv"
    peptide_path = study / "peptide_effects.tsv"
    protein_path = study / "protein_effects.tsv"
    for path in (measurement_path, peptide_path, protein_path):
        if not path.is_file():
            raise SystemExit(f"Missing study benchmark output: {path}")

    measurement = pd.read_csv(measurement_path, sep="\t", low_memory=False)
    peptide = pd.read_csv(peptide_path, sep="\t", low_memory=False)
    protein = pd.read_csv(protein_path, sep="\t", low_memory=False)

    required_measurement = {
        "run_name", "condition", "sequence", "charge", "realized_input_events", "intensity",
        "rt_error_seconds", "im_error",
    }
    missing = required_measurement - set(measurement.columns)
    if missing:
        raise SystemExit(f"precursor_measurements.tsv missing columns: {sorted(missing)}")

    # Per-run absolute-signal fidelity. The arbitrary intensity scale changes the intercept,
    # while correlation and slope quantify whether OpenDIA tracks TimSim abundance monotonically.
    run_rows: list[dict[str, float | int | str]] = []
    for run_name, group in measurement.groupby("run_name", sort=True):
        x = pd.to_numeric(group["realized_input_events"], errors="coerce")
        y = pd.to_numeric(group["intensity"], errors="coerce")
        mask = x.gt(0) & y.gt(0) & np.isfinite(x) & np.isfinite(y)
        lx = np.log2(x[mask])
        ly = np.log2(y[mask])
        metric = effect_metric(lx, ly)
        run_rows.append({
            "RunName": str(run_name),
            "Condition": str(group["condition"].iloc[0]),
            "N": metric.n,
            "Log2SignalPearsonR": metric.pearson_r,
            "Log2SignalSpearmanRho": metric.spearman_rho,
            "Log2SignalSlope": metric.slope,
            "MedianAbsRTErrorSeconds": float(pd.to_numeric(group["rt_error_seconds"], errors="coerce").abs().median()),
            "MedianAbsIMError": float(pd.to_numeric(group["im_error"], errors="coerce").abs().median()),
        })
    run_signal = pd.DataFrame(run_rows)
    run_signal.to_csv(out / "run_signal_fidelity.tsv", sep="\t", index=False)

    # Add per-precursor localization/signal diagnostics to effect residuals.
    measurement["realized_input_events"] = pd.to_numeric(measurement["realized_input_events"], errors="coerce")
    measurement["intensity"] = pd.to_numeric(measurement["intensity"], errors="coerce")
    per_precursor = (
        measurement.groupby(["sequence", "charge"], as_index=False)
        .agg(
            MedianRealizedEvents=("realized_input_events", "median"),
            MedianObservedIntensity=("intensity", "median"),
            MeanAbsRTErrorSeconds=("rt_error_seconds", lambda x: pd.to_numeric(x, errors="coerce").abs().mean()),
            MeanAbsIMError=("im_error", lambda x: pd.to_numeric(x, errors="coerce").abs().mean()),
        )
        .rename(columns={"sequence": "PeptideSequence", "charge": "PrecursorCharge"})
    )

    peptide["EffectResidual"] = pd.to_numeric(peptide["ObservedLog2FC"], errors="coerce") - pd.to_numeric(peptide["RealizedConditionLog2FC"], errors="coerce")
    peptide["AbsEffectResidual"] = peptide["EffectResidual"].abs()
    peptide = peptide.merge(per_precursor, on=["PeptideSequence", "PrecursorCharge"], how="left", validate="one_to_one")
    peptide["MedianLog2RealizedEvents"] = np.log2(pd.to_numeric(peptide["MedianRealizedEvents"], errors="coerce").where(lambda s: s > 0))
    peptide_outliers = peptide.sort_values("AbsEffectResidual", ascending=False).head(max(1, args.top))
    peptide_outliers.to_csv(out / "peptide_quant_outliers.tsv", sep="\t", index=False)

    protein["EffectResidual"] = pd.to_numeric(protein["ObservedLog2FC"], errors="coerce") - pd.to_numeric(protein["RealizedConditionLog2FC"], errors="coerce")
    protein["AbsEffectResidual"] = protein["EffectResidual"].abs()
    protein_outliers = protein.sort_values("AbsEffectResidual", ascending=False).head(max(1, args.top))
    protein_outliers.to_csv(out / "protein_quant_outliers.tsv", sep="\t", index=False)

    peptide_metric = effect_metric(peptide["RealizedConditionLog2FC"], peptide["ObservedLog2FC"])
    peptide_trimmed = trimmed_metric(peptide, "RealizedConditionLog2FC", "ObservedLog2FC")
    protein_metric = effect_metric(protein["RealizedConditionLog2FC"], protein["ObservedLog2FC"])
    protein_trimmed = trimmed_metric(protein, "RealizedConditionLog2FC", "ObservedLog2FC")

    protein_strata_rows = []
    selected = pd.to_numeric(protein["SelectedPrecursors"], errors="coerce")
    strata = {
        "1": selected.eq(1),
        "2": selected.eq(2),
        "3+": selected.ge(3),
    }
    for label, mask in strata.items():
        subset = protein[mask]
        metric = effect_metric(subset["RealizedConditionLog2FC"], subset["ObservedLog2FC"])
        protein_strata_rows.append({"SelectedPrecursors": label, **asdict(metric)})
    protein_strata = pd.DataFrame(protein_strata_rows)
    protein_strata.to_csv(out / "protein_effect_by_precursor_count.tsv", sep="\t", index=False)

    peptide_residual = residual_summary(peptide["EffectResidual"])
    protein_residual = residual_summary(protein["EffectResidual"])
    associations = {
        "peptide_abs_residual_vs_rt_error_spearman": safe_spearman(peptide["AbsEffectResidual"], peptide["MeanAbsRTErrorSeconds"]),
        "peptide_abs_residual_vs_im_error_spearman": safe_spearman(peptide["AbsEffectResidual"], peptide["MeanAbsIMError"]),
        "peptide_abs_residual_vs_signal_spearman": safe_spearman(peptide["AbsEffectResidual"], peptide["MedianLog2RealizedEvents"]),
    }

    summary = {
        "peptide_effect": asdict(peptide_metric),
        "peptide_effect_trimmed_top_1pct_residual": asdict(peptide_trimmed),
        "peptide_residual": peptide_residual,
        "protein_effect": asdict(protein_metric),
        "protein_effect_trimmed_top_1pct_residual": asdict(protein_trimmed),
        "protein_residual": protein_residual,
        "run_signal": {
            "mean_pearson_r": float(pd.to_numeric(run_signal["Log2SignalPearsonR"], errors="coerce").mean()),
            "min_pearson_r": float(pd.to_numeric(run_signal["Log2SignalPearsonR"], errors="coerce").min()),
            "mean_slope": float(pd.to_numeric(run_signal["Log2SignalSlope"], errors="coerce").mean()),
            "min_slope": float(pd.to_numeric(run_signal["Log2SignalSlope"], errors="coerce").min()),
            "max_slope": float(pd.to_numeric(run_signal["Log2SignalSlope"], errors="coerce").max()),
        },
        "associations": associations,
    }
    (out / "quantification_diagnostics.json").write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")

    lines = [
        "# Study quantification diagnostics",
        "",
        "This report separates per-run abundance tracking from condition-effect recovery and highlights heavy-tail failures.",
        "",
        "## Effect recovery",
        "",
        "| Level | n | Pearson r | slope | intercept | MAE | RMSE |",
        "|---|---:|---:|---:|---:|---:|---:|",
        f"| Peptide | {peptide_metric.n} | {peptide_metric.pearson_r:.4f} | {peptide_metric.slope:.4f} | {peptide_metric.intercept:.4f} | {peptide_metric.mae:.4f} | {peptide_metric.rmse:.4f} |",
        f"| Peptide, excluding worst 1% residuals | {peptide_trimmed.n} | {peptide_trimmed.pearson_r:.4f} | {peptide_trimmed.slope:.4f} | {peptide_trimmed.intercept:.4f} | {peptide_trimmed.mae:.4f} | {peptide_trimmed.rmse:.4f} |",
        f"| Protein | {protein_metric.n} | {protein_metric.pearson_r:.4f} | {protein_metric.slope:.4f} | {protein_metric.intercept:.4f} | {protein_metric.mae:.4f} | {protein_metric.rmse:.4f} |",
        f"| Protein, excluding worst 1% residuals | {protein_trimmed.n} | {protein_trimmed.pearson_r:.4f} | {protein_trimmed.slope:.4f} | {protein_trimmed.intercept:.4f} | {protein_trimmed.mae:.4f} | {protein_trimmed.rmse:.4f} |",
        "",
        "The regression is `ObservedLog2FC = intercept + slope × RealizedConditionLog2FC`; slope near 1 indicates unbiased effect amplitude.",
        "",
        "## Residual tails",
        "",
        f"Peptide absolute residual: median={peptide_residual.get('median_abs', math.nan):.4f}, p95={peptide_residual.get('p95_abs', math.nan):.4f}, p99={peptide_residual.get('p99_abs', math.nan):.4f}; "
        f">0.5={peptide_residual.get('gt_0_5', 0)}, >1.0={peptide_residual.get('gt_1_0', 0)}.",
        "",
        f"Protein absolute residual: median={protein_residual.get('median_abs', math.nan):.4f}, p95={protein_residual.get('p95_abs', math.nan):.4f}, p99={protein_residual.get('p99_abs', math.nan):.4f}; "
        f">0.5={protein_residual.get('gt_0_5', 0)}, >1.0={protein_residual.get('gt_1_0', 0)}.",
        "",
        "## Per-run abundance tracking",
        "",
        f"Mean log2 TimSim-input vs OpenDIA-intensity Pearson r={summary['run_signal']['mean_pearson_r']:.4f}; minimum run r={summary['run_signal']['min_pearson_r']:.4f}. ",
        f"Mean regression slope={summary['run_signal']['mean_slope']:.4f} (range {summary['run_signal']['min_slope']:.4f}–{summary['run_signal']['max_slope']:.4f}).",
        "",
        "See `run_signal_fidelity.tsv` for every run.",
        "",
        "## Residual associations",
        "",
        f"Spearman |effect residual| vs mean absolute RT error: {associations['peptide_abs_residual_vs_rt_error_spearman']:.4f}.",
        f"Spearman |effect residual| vs mean absolute IM error: {associations['peptide_abs_residual_vs_im_error_spearman']:.4f}.",
        f"Spearman |effect residual| vs median realized signal: {associations['peptide_abs_residual_vs_signal_spearman']:.4f}.",
        "",
        "## Protein precursor-count strata",
        "",
        "`protein_effect_by_precursor_count.tsv` tests whether protein errors are concentrated in proteins represented by only one selected precursor.",
        "",
        f"The top {max(1, args.top)} peptide and protein effect residuals are written to `peptide_quant_outliers.tsv` and `protein_quant_outliers.tsv`.",
        "",
    ]
    (out / "quantification_diagnostics.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"Quantification diagnostics: {out}")
    print(f"Peptide effect: r={peptide_metric.pearson_r:.4f}, slope={peptide_metric.slope:.4f}, MAE={peptide_metric.mae:.4f}")
    print(f"Peptide trimmed 1%: r={peptide_trimmed.pearson_r:.4f}, slope={peptide_trimmed.slope:.4f}, MAE={peptide_trimmed.mae:.4f}")
    print(f"Per-run signal fidelity: mean r={summary['run_signal']['mean_pearson_r']:.4f}, mean slope={summary['run_signal']['mean_slope']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
