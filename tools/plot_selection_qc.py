#!/usr/bin/env python3
"""Plot the frozen high-signal selection distributions."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def save_hist(frame: pd.DataFrame, column: str, xlabel: str, title: str, path: Path, *, log10: bool = False) -> None:
    values = pd.to_numeric(frame[column], errors="coerce").dropna().to_numpy(dtype=float)
    if log10:
        values = np.log10(np.clip(values, 1.0, None))
        xlabel = f"log10({xlabel})"
    plt.figure(figsize=(7.2, 4.6))
    plt.hist(values, bins=30)
    plt.xlabel(xlabel)
    plt.ylabel("Selected precursor groups")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def save_bar(counts: pd.Series, xlabel: str, title: str, path: Path) -> None:
    plt.figure(figsize=(8.0, 4.6))
    x = np.arange(len(counts))
    plt.bar(x, counts.to_numpy())
    plt.xticks(x, [str(value) for value in counts.index], rotation=90 if len(counts) > 12 else 0)
    plt.xlabel(xlabel)
    plt.ylabel("Selected precursor groups")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("build_dir", type=Path)
    args = parser.parse_args()
    root = args.build_dir.resolve()
    selected_path = root / "OpenSwathTimSim.high_signal_selection.tsv"
    selected = pd.read_csv(selected_path, sep="\t", low_memory=False)
    if selected.empty:
        raise SystemExit(f"Selection TSV is empty: {selected_path}")

    out = root / "selection_qc" / "plots"
    out.mkdir(parents=True, exist_ok=True)
    save_hist(selected, "min_realized_event_proxy", "minimum realized event proxy", "Minimum realized precursor signal across three runs", out / "01_realized_signal.png", log10=True)
    save_hist(selected, "min_frame_abundance_sum", "minimum frame abundance sum", "Retained chromatographic mass", out / "02_frame_abundance_sum.png")
    save_hist(selected, "min_scan_abundance_sum", "minimum scan abundance sum", "Retained mobility mass", out / "03_scan_abundance_sum.png")
    save_hist(selected, "min_ion_relative_abundance", "minimum charge-state relative abundance", "Selected charge-state abundance", out / "04_charge_fraction.png")
    save_hist(selected, "assay_rt", "assay RT (s)", "RT coverage of frozen correctness set", out / "05_rt.png")
    save_hist(selected, "precursor_mz", "precursor m/z", "Precursor m/z coverage of frozen correctness set", out / "06_precursor_mz.png")
    save_hist(selected, "assay_im", "assay inverse mobility (1/K0)", "Ion-mobility coverage of frozen correctness set", out / "07_precursor_im.png")

    charge_counts = selected["charge"].astype(int).value_counts().sort_index()
    save_bar(charge_counts, "precursor charge", "Charge-state coverage", out / "08_charge.png")
    map_counts = selected["swath_window_index"].astype(int).value_counts().sort_index()
    save_bar(map_counts, "diaPASEF SWATH map", "SWATH-map coverage", out / "09_swath_map.png")

    print(f"Wrote selection QC plots under {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
