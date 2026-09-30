#!/usr/bin/env python3
"""Benchmark raw TDF oracle signal between TimSim truth and OpenDIA quantification."""
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
class Metric:
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


def metric(x: pd.Series, y: pd.Series) -> Metric:
    xa, ya = finite_xy(x, y)
    if len(xa) == 0:
        return Metric(0, *(math.nan for _ in range(7)))
    error = ya - xa
    if len(xa) > 1 and not np.allclose(xa, xa[0]) and not np.allclose(ya, ya[0]):
        pearson = float(stats.pearsonr(xa, ya).statistic)
        spearman = float(stats.spearmanr(xa, ya).statistic)
        regression = stats.linregress(xa, ya)
        slope = float(regression.slope)
        intercept = float(regression.intercept)
    else:
        pearson = spearman = slope = intercept = math.nan
    return Metric(
        n=int(len(xa)),
        bias=float(np.mean(error)),
        mae=float(np.mean(np.abs(error))),
        rmse=float(np.sqrt(np.mean(error**2))),
        pearson_r=pearson,
        spearman_rho=spearman,
        slope=slope,
        intercept=intercept,
    )


def positive_log2(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    return np.log2(values.where(values > 0))


def safe_log2_ratio(treatment: float, control: float) -> float:
    if not math.isfinite(treatment) or not math.isfinite(control) or treatment <= 0 or control <= 0:
        return math.nan
    return math.log2(treatment / control)


def trimmed_metric(table: pd.DataFrame, x: str, y: str, fraction: float = 0.01) -> Metric:
    work = table[[x, y]].copy()
    work["residual"] = pd.to_numeric(work[y], errors="coerce") - pd.to_numeric(work[x], errors="coerce")
    work = work[np.isfinite(work["residual"])].copy()
    if len(work) < 2:
        return metric(work[x], work[y])
    keep = max(1, int(math.floor(len(work) * (1.0 - fraction))))
    work = work.assign(abs_residual=work["residual"].abs()).nsmallest(keep, "abs_residual")
    return metric(work[x], work[y])


def format_metric(value: Metric) -> str:
    return (
        f"n={value.n}, r={value.pearson_r:.4f}, slope={value.slope:.4f}, "
        f"MAE={value.mae:.4f}, RMSE={value.rmse:.4f}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oracle-dir", type=Path, required=True)
    parser.add_argument("--study-results-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()

    oracle = args.oracle_dir.resolve()
    study = args.study_results_dir.resolve()
    out = (args.out_dir or oracle).resolve()
    out.mkdir(parents=True, exist_ok=True)
    oracle_path = oracle / "raw_signal_oracle_measurements.tsv"
    measurement_path = study / "precursor_measurements.tsv"
    peptide_path = study / "peptide_effects.tsv"
    for path in (oracle_path, measurement_path, peptide_path):
        if not path.is_file():
            raise SystemExit(f"Missing benchmark input: {path}")

    raw = pd.read_csv(oracle_path, sep="\t", low_memory=False)
    observed = pd.read_csv(measurement_path, sep="\t", low_memory=False)
    peptide = pd.read_csv(peptide_path, sep="\t", low_memory=False)

    raw["sequence"] = raw["PeptideSequence"].astype(str).str.upper()
    raw["charge"] = pd.to_numeric(raw["PrecursorCharge"], errors="raise").astype(int)
    raw["run_name"] = raw["RunName"].astype(str)
    observed["sequence"] = observed["sequence"].astype(str).str.upper()
    observed["charge"] = pd.to_numeric(observed["charge"], errors="raise").astype(int)
    observed["run_name"] = observed["run_name"].astype(str)
    joined = raw.merge(
        observed[["run_name", "sequence", "charge", "condition", "intensity", "rt_error_seconds", "im_error"]],
        on=["run_name", "sequence", "charge"],
        how="left",
        validate="one_to_one",
    )
    joined.to_csv(out / "raw_oracle_openDIA_measurements.tsv", sep="\t", index=False)

    run_rows: list[dict[str, float | int | str]] = []
    for run_name, group in joined.groupby("run_name", sort=True):
        work = group.copy()
        input_metric = metric(positive_log2(work["RealizedInputEvents"]), positive_log2(work["RawOracleIntensity"]))
        proxy_metric = metric(positive_log2(work["RealizedEventProxy"]), positive_log2(work["RawOracleIntensity"]))
        raw_opendia_metric = metric(positive_log2(work["RawOracleIntensity"]), positive_log2(work["intensity"]))
        run_rows.append(
            {
                "RunName": run_name,
                "Condition": str(group["condition"].dropna().iloc[0]) if group["condition"].notna().any() else "unknown",
                "N": int(len(group)),
                "RawPositive": int(pd.to_numeric(group["RawOracleIntensity"], errors="coerce").gt(0).sum()),
                "InputVsRawPearsonR": input_metric.pearson_r,
                "InputVsRawSlope": input_metric.slope,
                "ProxyVsRawPearsonR": proxy_metric.pearson_r,
                "ProxyVsRawSlope": proxy_metric.slope,
                "RawVsOpenDIAPearsonR": raw_opendia_metric.pearson_r,
                "RawVsOpenDIASlope": raw_opendia_metric.slope,
            }
        )
    run_table = pd.DataFrame(run_rows)
    run_table.to_csv(out / "raw_oracle_run_metrics.tsv", sep="\t", index=False)

    raw_effect_rows = []
    for (sequence, charge), group in joined.groupby(["sequence", "charge"], sort=True):
        control = group[group["condition"].eq("control")]
        treatment = group[group["condition"].eq("treatment")]
        control_raw = pd.to_numeric(control["RawOracleIntensity"], errors="coerce")
        treatment_raw = pd.to_numeric(treatment["RawOracleIntensity"], errors="coerce")
        control_raw = control_raw[control_raw > 0]
        treatment_raw = treatment_raw[treatment_raw > 0]
        raw_log2fc = (
            safe_log2_ratio(float(treatment_raw.mean()), float(control_raw.mean()))
            if len(control_raw) and len(treatment_raw)
            else math.nan
        )
        raw_effect_rows.append(
            {
                "PeptideSequence": sequence,
                "PrecursorCharge": int(charge),
                "RawOracleLog2FC": raw_log2fc,
                "RawControlRuns": int(len(control_raw)),
                "RawTreatmentRuns": int(len(treatment_raw)),
                "RawMedianIntensity": float(pd.to_numeric(group["RawOracleIntensity"], errors="coerce").median()),
                "MatchedTransitionsMedian": float(pd.to_numeric(group["MatchedTransitions"], errors="coerce").median()),
            }
        )
    raw_effects = pd.DataFrame(raw_effect_rows)
    effects = peptide.merge(raw_effects, on=["PeptideSequence", "PrecursorCharge"], how="left", validate="one_to_one")
    effects.to_csv(out / "raw_oracle_peptide_effects.tsv", sep="\t", index=False)

    realized_raw = metric(effects["RealizedConditionLog2FC"], effects["RawOracleLog2FC"])
    raw_opendia_effect = metric(effects["RawOracleLog2FC"], effects["ObservedLog2FC"])
    realized_raw_trimmed = trimmed_metric(effects, "RealizedConditionLog2FC", "RawOracleLog2FC")
    raw_opendia_trimmed = trimmed_metric(effects, "RawOracleLog2FC", "ObservedLog2FC")

    mean_run = {
        "input_vs_raw_pearson_r": float(pd.to_numeric(run_table["InputVsRawPearsonR"], errors="coerce").mean()),
        "input_vs_raw_slope": float(pd.to_numeric(run_table["InputVsRawSlope"], errors="coerce").mean()),
        "proxy_vs_raw_pearson_r": float(pd.to_numeric(run_table["ProxyVsRawPearsonR"], errors="coerce").mean()),
        "proxy_vs_raw_slope": float(pd.to_numeric(run_table["ProxyVsRawSlope"], errors="coerce").mean()),
        "raw_vs_opendia_pearson_r": float(pd.to_numeric(run_table["RawVsOpenDIAPearsonR"], errors="coerce").mean()),
        "raw_vs_opendia_slope": float(pd.to_numeric(run_table["RawVsOpenDIASlope"], errors="coerce").mean()),
    }
    summary = {
        "runs": int(len(run_table)),
        "target_run_rows": int(len(joined)),
        "raw_positive_fraction": float(pd.to_numeric(joined["RawOracleIntensity"], errors="coerce").gt(0).mean()),
        "mean_run_metrics": mean_run,
        "condition_effects": {
            "realized_vs_raw": asdict(realized_raw),
            "realized_vs_raw_trimmed_top_1pct": asdict(realized_raw_trimmed),
            "raw_vs_opendia": asdict(raw_opendia_effect),
            "raw_vs_opendia_trimmed_top_1pct": asdict(raw_opendia_trimmed),
        },
    }
    (out / "raw_signal_oracle_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# TimSim raw-signal oracle benchmark",
        "",
        "The raw oracle integrates the frozen target library directly from generated DIA-PASEF raw data at simulator-realized RT/IM coordinates. OpenDIA coordinates, scores, and peak boundaries are not used during extraction.",
        "",
        "## Per-run absolute signal scaling",
        "",
        "| Run | Condition | Input→raw r | Input→raw slope | Event-proxy→raw r | Event-proxy→raw slope | Raw→OpenDIA r | Raw→OpenDIA slope |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in run_table.itertuples(index=False):
        lines.append(
            f"| {row.RunName} | {row.Condition} | {row.InputVsRawPearsonR:.4f} | {row.InputVsRawSlope:.4f} | "
            f"{row.ProxyVsRawPearsonR:.4f} | {row.ProxyVsRawSlope:.4f} | "
            f"{row.RawVsOpenDIAPearsonR:.4f} | {row.RawVsOpenDIASlope:.4f} |"
        )
    lines.extend(
        [
            "",
            f"Mean input→raw slope: **{mean_run['input_vs_raw_slope']:.4f}**; mean event-proxy→raw slope: **{mean_run['proxy_vs_raw_slope']:.4f}**; mean raw→OpenDIA slope: **{mean_run['raw_vs_opendia_slope']:.4f}**.",
            "",
            "`RealizedEventProxy` is the preferred charge-state-specific absolute-signal comparator because it includes peptide events, chromatographic frame mass, ion relative abundance, and mobility scan mass.",
            "",
            "## Condition-effect decomposition",
            "",
            "| Comparison | n | Pearson r | slope | MAE | RMSE |",
            "|---|---:|---:|---:|---:|---:|",
            f"| Realized input effect → raw oracle effect | {realized_raw.n} | {realized_raw.pearson_r:.4f} | {realized_raw.slope:.4f} | {realized_raw.mae:.4f} | {realized_raw.rmse:.4f} |",
            f"| Realized → raw, excluding worst 1% residuals | {realized_raw_trimmed.n} | {realized_raw_trimmed.pearson_r:.4f} | {realized_raw_trimmed.slope:.4f} | {realized_raw_trimmed.mae:.4f} | {realized_raw_trimmed.rmse:.4f} |",
            f"| Raw oracle effect → OpenDIA effect | {raw_opendia_effect.n} | {raw_opendia_effect.pearson_r:.4f} | {raw_opendia_effect.slope:.4f} | {raw_opendia_effect.mae:.4f} | {raw_opendia_effect.rmse:.4f} |",
            f"| Raw → OpenDIA, excluding worst 1% residuals | {raw_opendia_trimmed.n} | {raw_opendia_trimmed.pearson_r:.4f} | {raw_opendia_trimmed.slope:.4f} | {raw_opendia_trimmed.mae:.4f} | {raw_opendia_trimmed.rmse:.4f} |",
            "",
            "## Interpretation",
            "",
            "- If realized→raw has slope near 1 while raw→OpenDIA is compressed, the response compression enters during OpenDIA extraction/quantification.",
            "- If realized→raw is already compressed, raw realized signal is the appropriate measurement-level truth and the simulator/raw response should be characterized before changing OpenDIA.",
            "- If both slopes are near 1 but a heavy residual tail remains, inspect localization and interference for those specific precursors.",
            "",
            "This oracle is diagnostic only and must not be used to select targets, construct the assay library, or tune OpenDIA on the canonical benchmark.",
            "",
            "## Outputs",
            "",
            "- `raw_signal_oracle_measurements.tsv`: truth-centered raw signal joined to target/run identity.",
            "- `raw_signal_oracle_transitions.tsv`: per-transition raw integrated signal.",
            "- `raw_oracle_openDIA_measurements.tsv`: raw-oracle measurements joined to OpenDIA target-only measurements.",
            "- `raw_oracle_run_metrics.tsv`: per-run input/raw/OpenDIA regressions.",
            "- `raw_oracle_peptide_effects.tsv`: realized, raw-oracle, and OpenDIA condition effects.",
            "- `raw_signal_oracle_summary.json`: machine-readable metrics.",
            "",
        ]
    )
    (out / "raw_signal_oracle_report.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"Raw-signal oracle benchmark: {out}")
    print(f"Realized -> raw effect: {format_metric(realized_raw)}")
    print(f"Raw -> OpenDIA effect: {format_metric(raw_opendia_effect)}")
    print(
        "Mean per-run slopes: "
        f"event-proxy->raw={mean_run['proxy_vs_raw_slope']:.4f}, "
        f"raw->OpenDIA={mean_run['raw_vs_opendia_slope']:.4f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
