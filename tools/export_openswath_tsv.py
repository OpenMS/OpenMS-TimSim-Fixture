#!/usr/bin/env python3
"""Export a compact OpenSWATH assay-generator TSV from TimSim synthetic_data.db.

The script intentionally supports only unmodified canonical peptides. This keeps
m/z calculations transparent and appropriate for a small OpenMS test fixture.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

PROTON = 1.007276466621
WATER = 18.0105646837
AA_MASS = {
    "A": 71.037113805, "R": 156.101111050, "N": 114.042927470,
    "D": 115.026943065, "C": 103.009184505, "E": 129.042593135,
    "Q": 128.058577540, "G": 57.021463735, "H": 137.058911875,
    "I": 113.084063975, "L": 113.084063975, "K": 128.094963015,
    "M": 131.040484645, "F": 147.068413945, "P": 97.052763875,
    "S": 87.032028435, "T": 101.047678505, "W": 186.079312980,
    "Y": 163.063328575, "V": 99.068413945,
}

OUTPUT_COLUMNS = [
    "PrecursorMz", "ProductMz", "LibraryIntensity", "NormalizedRetentionTime",
    "PrecursorIonMobility", "ProteinId", "PeptideSequence",
    "ModifiedPeptideSequence", "PeptideGroupLabel", "PrecursorCharge",
    "ProductCharge", "FragmentType", "FragmentSeriesNumber",
    "TransitionGroupId", "TransitionId", "Decoy", "DetectingTransition",
    "IdentifyingTransition", "QuantifyingTransition",
]


@dataclass(frozen=True)
class Window:
    mz_lower: float
    mz_upper: float
    im_lower: float
    im_upper: float


@dataclass
class Candidate:
    ion_id: int
    peptide_id: int
    sequence: str
    protein: str
    precursor_mz: float
    precursor_charge: int
    precursor_im: float
    rt_seconds: float
    fragments: list[dict[str, Any]]

    @property
    def total_fragment_intensity(self) -> float:
        return sum(float(fragment["intensity"]) for fragment in self.fragments)


def quote_identifier(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def table_names(connection: sqlite3.Connection) -> set[str]:
    return {row[0] for row in connection.execute(
        "SELECT name FROM sqlite_master WHERE type='table'"
    )}


def columns(connection: sqlite3.Connection, table: str) -> list[str]:
    return [row[1] for row in connection.execute(
        f"PRAGMA table_info({quote_identifier(table)})"
    )]


def first_present(options: Iterable[str], available: Iterable[str], *, label: str) -> str:
    available_set = set(available)
    for option in options:
        if option in available_set:
            return option
    raise RuntimeError(f"Could not identify {label}; available columns: {sorted(available_set)}")


def optional_present(options: Iterable[str], available: Iterable[str]) -> str | None:
    available_set = set(available)
    return next((option for option in options if option in available_set), None)


def parse_array(value: Any) -> list[float]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return list(value)
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    text = str(value).strip()
    if not text:
        return []
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        # Older SQLite exports may use Python/NumPy-like bracketed strings.
        text = text.strip("[]")
        if not text.strip():
            return []
        return [float(token) for token in text.replace(",", " ").split()]
    if not isinstance(parsed, list):
        raise ValueError(f"Expected a list, received: {type(parsed).__name__}")
    return list(parsed)


def decode_prosit_index(index: int) -> tuple[str, int, int]:
    """Decode TimSim's 3 x (29 y + 29 b) flattened intensity layout."""
    if not 0 <= index < 174:
        raise ValueError(f"Unexpected Prosit index: {index}")
    product_charge = index // 58 + 1
    within_charge = index % 58
    if within_charge < 29:
        return "y", within_charge + 1, product_charge
    return "b", within_charge - 29 + 1, product_charge


def fragment_mz(sequence: str, ion_type: str, ordinal: int, charge: int) -> float:
    if ordinal < 1 or ordinal >= len(sequence):
        raise ValueError("Fragment ordinal must be between 1 and peptide_length - 1")
    residues = sequence[:ordinal] if ion_type == "b" else sequence[-ordinal:]
    neutral_mass = sum(AA_MASS[aa] for aa in residues)
    if ion_type == "y":
        neutral_mass += WATER
    return (neutral_mass + charge * PROTON) / charge


def load_windows(connection: sqlite3.Connection) -> list[Window]:
    names = table_names(connection)
    if "dia_ms_ms_windows" not in names or "scans" not in names:
        print("WARNING: DIA window or scan metadata is absent; skipping window-margin filtering.", file=sys.stderr)
        return []

    window_cols = columns(connection, "dia_ms_ms_windows")
    scan_cols = columns(connection, "scans")
    start_col = first_present(("scan_start", "ScanNumBegin"), window_cols, label="window scan start")
    end_col = first_present(("scan_end", "ScanNumEnd"), window_cols, label="window scan end")
    center_col = first_present(("isolation_mz", "IsolationMz"), window_cols, label="window isolation m/z")
    width_col = first_present(("isolation_width", "IsolationWidth"), window_cols, label="window isolation width")
    scan_id_col = first_present(("scan", "scan_id", "Scan"), scan_cols, label="scan index")
    mobility_col = first_present(
        ("mobility", "inv_mobility", "one_over_k0", "inverse_mobility"),
        scan_cols,
        label="scan mobility",
    )

    scan_to_im = {
        int(row[0]): float(row[1])
        for row in connection.execute(
            f"SELECT {quote_identifier(scan_id_col)}, {quote_identifier(mobility_col)} "
            f"FROM {quote_identifier('scans')}"
        )
    }
    if not scan_to_im:
        return []

    windows: list[Window] = []
    query = (
        f"SELECT {quote_identifier(start_col)}, {quote_identifier(end_col)}, "
        f"{quote_identifier(center_col)}, {quote_identifier(width_col)} "
        f"FROM {quote_identifier('dia_ms_ms_windows')}"
    )
    available_scans = sorted(scan_to_im)
    for scan_start, scan_end, center, width in connection.execute(query):
        start = int(scan_start)
        end = int(scan_end)
        # TimSim scan_end may be exclusive. Use the nearest available boundary scans.
        start_key = min(available_scans, key=lambda value: abs(value - start))
        end_key = min(available_scans, key=lambda value: abs(value - max(start, end - 1)))
        im_a, im_b = scan_to_im[start_key], scan_to_im[end_key]
        center = float(center)
        width = float(width)
        windows.append(Window(
            center - width / 2.0,
            center + width / 2.0,
            min(im_a, im_b),
            max(im_a, im_b),
        ))
    return windows


def inside_window(candidate: Candidate, windows: list[Window], mz_margin: float, im_margin: float) -> bool:
    if not windows:
        return True
    return any(
        window.mz_lower + mz_margin < candidate.precursor_mz < window.mz_upper - mz_margin
        and window.im_lower + im_margin < candidate.precursor_im < window.im_upper - im_margin
        for window in windows
    )


def read_candidates(
    connection: sqlite3.Connection,
    min_product_mz: float,
    max_product_mz: float,
    minimum_fragments: int,
) -> list[Candidate]:
    names = table_names(connection)
    required = {"peptides", "ions", "fragment_ions"}
    missing = required - names
    if missing:
        raise RuntimeError(f"Missing TimSim tables: {sorted(missing)}; found {sorted(names)}")

    peptide_cols = columns(connection, "peptides")
    ion_cols = columns(connection, "ions")
    fragment_cols = columns(connection, "fragment_ions")

    peptide_id_col = first_present(("peptide_id", "id"), peptide_cols, label="peptide ID")
    sequence_col = first_present(("sequence", "peptide"), peptide_cols, label="peptide sequence")
    protein_col = optional_present(("protein", "protein_id", "protein_name"), peptide_cols)
    rt_col = first_present(
        ("retention_time_gru_predictor", "retention_time", "rt", "predicted_rt"),
        peptide_cols,
        label="retention time",
    )
    decoy_col = optional_present(("decoy", "is_decoy"), peptide_cols)

    peptide_select = [peptide_id_col, sequence_col, rt_col]
    if protein_col:
        peptide_select.append(protein_col)
    if decoy_col:
        peptide_select.append(decoy_col)
    peptide_query = "SELECT " + ", ".join(quote_identifier(value) for value in peptide_select) + " FROM peptides"

    peptides: dict[int, dict[str, Any]] = {}
    for row in connection.execute(peptide_query):
        data = dict(zip(peptide_select, row))
        peptide_id = int(data[peptide_id_col])
        if decoy_col and bool(data.get(decoy_col)):
            continue
        sequence = str(data[sequence_col]).strip().upper()
        if not sequence or any(aa not in AA_MASS for aa in sequence):
            continue
        peptides[peptide_id] = {
            "sequence": sequence,
            "rt": float(data[rt_col]),
            "protein": str(data.get(protein_col) or f"TIMSIM_PROTEIN_{peptide_id}"),
        }

    ion_id_col = first_present(("ion_id", "id"), ion_cols, label="ion ID")
    ion_peptide_col = first_present(("peptide_id",), ion_cols, label="ion peptide ID")
    charge_col = first_present(("charge", "precursor_charge"), ion_cols, label="precursor charge")
    precursor_mz_col = first_present(("mz", "precursor_mz"), ion_cols, label="precursor m/z")
    im_col = first_present(
        ("inv_mobility_gru_predictor", "inverse_mobility", "ion_mobility", "mobility"),
        ion_cols,
        label="precursor ion mobility",
    )

    ions: dict[int, dict[str, Any]] = {}
    ion_query = (
        f"SELECT {quote_identifier(ion_id_col)}, {quote_identifier(ion_peptide_col)}, "
        f"{quote_identifier(charge_col)}, {quote_identifier(precursor_mz_col)}, {quote_identifier(im_col)} "
        f"FROM ions"
    )
    for ion_id, peptide_id, charge, precursor_mz, mobility in connection.execute(ion_query):
        peptide_id = int(peptide_id)
        if peptide_id not in peptides:
            continue
        ions[int(ion_id)] = {
            "peptide_id": peptide_id,
            "charge": int(charge),
            "mz": float(precursor_mz),
            "im": float(mobility),
        }

    fragment_ion_id_col = first_present(("ion_id",), fragment_cols, label="fragment ion ID")
    indices_col = first_present(("indices", "sparse_indices"), fragment_cols, label="fragment indices")
    values_col = first_present(("values", "sparse_values"), fragment_cols, label="fragment values")

    # One ion may have predictions at multiple collision energies. Keep the row
    # with the greatest total predicted intensity for a deterministic library.
    best_fragment_row: dict[int, tuple[float, list[int], list[float]]] = {}
    fragment_query = (
        f"SELECT {quote_identifier(fragment_ion_id_col)}, {quote_identifier(indices_col)}, "
        f"{quote_identifier(values_col)} FROM fragment_ions"
    )
    for ion_id, indices_raw, values_raw in connection.execute(fragment_query):
        ion_id = int(ion_id)
        if ion_id not in ions:
            continue
        indices = [int(value) for value in parse_array(indices_raw)]
        values = [float(value) for value in parse_array(values_raw)]
        if len(indices) != len(values):
            raise RuntimeError(f"Fragment index/value length mismatch for ion_id={ion_id}")
        score = sum(max(0.0, value) for value in values)
        previous = best_fragment_row.get(ion_id)
        if previous is None or score > previous[0]:
            best_fragment_row[ion_id] = (score, indices, values)

    candidates: list[Candidate] = []
    for ion_id, ion in ions.items():
        if ion_id not in best_fragment_row:
            continue
        peptide = peptides[ion["peptide_id"]]
        sequence = peptide["sequence"]
        fragments: list[dict[str, Any]] = []
        _, indices, values = best_fragment_row[ion_id]
        for index, intensity in zip(indices, values):
            if intensity <= 0:
                continue
            ion_type, ordinal, product_charge = decode_prosit_index(index)
            if ordinal >= len(sequence) or product_charge > ion["charge"]:
                continue
            mz = fragment_mz(sequence, ion_type, ordinal, product_charge)
            if min_product_mz <= mz <= max_product_mz:
                fragments.append({
                    "type": ion_type,
                    "ordinal": ordinal,
                    "charge": product_charge,
                    "mz": mz,
                    "intensity": intensity,
                })
        # Avoid repeated annotations and preserve only the strongest prediction.
        unique: dict[tuple[str, int, int], dict[str, Any]] = {}
        for fragment in fragments:
            key = (fragment["type"], fragment["ordinal"], fragment["charge"])
            if key not in unique or fragment["intensity"] > unique[key]["intensity"]:
                unique[key] = fragment
        fragments = sorted(unique.values(), key=lambda item: item["intensity"], reverse=True)
        if len(fragments) < minimum_fragments:
            continue
        candidates.append(Candidate(
            ion_id=ion_id,
            peptide_id=ion["peptide_id"],
            sequence=sequence,
            protein=peptide["protein"],
            precursor_mz=ion["mz"],
            precursor_charge=ion["charge"],
            precursor_im=ion["im"],
            rt_seconds=peptide["rt"],
            fragments=fragments,
        ))
    return candidates


def choose_candidates(candidates: list[Candidate], max_precursors: int) -> list[Candidate]:
    # Retain one strong charge state per sequence, then select across the full RT
    # range. A strength-only selection can cluster calibrants in one RT region
    # and make automatic iRT fitting fragile even when many candidates exist.
    best_by_sequence: dict[str, Candidate] = {}
    for candidate in candidates:
        previous = best_by_sequence.get(candidate.sequence)
        if previous is None or candidate.total_fragment_intensity > previous.total_fragment_intensity:
            best_by_sequence[candidate.sequence] = candidate

    ordered = sorted(best_by_sequence.values(), key=lambda item: (item.rt_seconds, item.sequence))
    if len(ordered) <= max_precursors:
        return ordered

    # Split the RT-ordered candidates into contiguous quantile blocks and keep
    # the strongest candidate in each block. This preserves RT coverage while
    # preferring transitions with robust predicted signal.
    chosen: list[Candidate] = []
    total = len(ordered)
    for block in range(max_precursors):
        start = block * total // max_precursors
        end = (block + 1) * total // max_precursors
        chunk = ordered[start:max(start + 1, end)]
        chosen.append(max(chunk, key=lambda item: item.total_fragment_intensity))
    return sorted(chosen, key=lambda item: (item.rt_seconds, item.sequence))


def write_outputs(
    candidates: list[Candidate],
    output_tsv: Path,
    truth_tsv: Path,
    transitions_per_precursor: int,
    gradient_length: float,
    rt_mode: str,
) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    truth_rows: list[dict[str, Any]] = []

    for candidate in candidates:
        group_id = f"TIMSIM_{candidate.sequence}_{candidate.precursor_charge}"
        selected = candidate.fragments[:transitions_per_precursor]
        max_intensity = max(fragment["intensity"] for fragment in selected)
        normalized_rt = (
            candidate.rt_seconds if rt_mode == "seconds"
            else 100.0 * candidate.rt_seconds / gradient_length
        )
        truth_rows.append({
            "TransitionGroupId": group_id,
            "PeptideSequence": candidate.sequence,
            "PrecursorCharge": candidate.precursor_charge,
            "PrecursorMz": f"{candidate.precursor_mz:.8f}",
            "PrecursorIonMobility": f"{candidate.precursor_im:.8f}",
            "RetentionTimeSeconds": f"{candidate.rt_seconds:.6f}",
            "NormalizedRetentionTime": f"{normalized_rt:.6f}",
            "NumberOfTransitions": len(selected),
        })
        for fragment in selected:
            annotation = f"{fragment['type']}{fragment['ordinal']}^{fragment['charge']}"
            transition_id = f"{group_id}_{annotation}"
            rows.append({
                "PrecursorMz": f"{candidate.precursor_mz:.8f}",
                "ProductMz": f"{fragment['mz']:.8f}",
                "LibraryIntensity": f"{10000.0 * fragment['intensity'] / max_intensity:.6f}",
                "NormalizedRetentionTime": f"{normalized_rt:.6f}",
                "PrecursorIonMobility": f"{candidate.precursor_im:.8f}",
                "ProteinId": candidate.protein,
                "PeptideSequence": candidate.sequence,
                "ModifiedPeptideSequence": candidate.sequence,
                "PeptideGroupLabel": f"TIMSIM_{candidate.sequence}",
                "PrecursorCharge": candidate.precursor_charge,
                "ProductCharge": fragment["charge"],
                "FragmentType": fragment["type"],
                "FragmentSeriesNumber": fragment["ordinal"],
                "TransitionGroupId": group_id,
                "TransitionId": transition_id,
                "Decoy": 0,
                "DetectingTransition": 1,
                "IdentifyingTransition": 0,
                "QuantifyingTransition": 1,
            })

    with output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    truth_columns = [
        "TransitionGroupId", "PeptideSequence", "PrecursorCharge", "PrecursorMz",
        "PrecursorIonMobility", "RetentionTimeSeconds", "NormalizedRetentionTime",
        "NumberOfTransitions",
    ]
    with truth_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=truth_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(truth_rows)



def load_selected_precursors(path: Path) -> list[tuple[str, int]]:
    """Load a frozen simulator-only precursor selection in selection-rank order."""
    if not path.is_file():
        raise FileNotFoundError(f"Selected precursor TSV does not exist: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"Selected precursor TSV has no rows: {path}")
    required = {"sequence", "charge"}
    missing = required - set(rows[0])
    if missing:
        raise RuntimeError(f"Selected precursor TSV is missing columns: {sorted(missing)}")
    if "selection_rank" in rows[0]:
        rows.sort(key=lambda row: int(row["selection_rank"]))
    keys = [(str(row["sequence"]).strip().upper(), int(row["charge"])) for row in rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("Selected precursor TSV contains duplicate sequence/charge keys")
    if len({sequence for sequence, _ in keys}) != len(keys):
        raise RuntimeError("Selected precursor TSV contains more than one charge state for a peptide sequence")
    return keys

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--truth-out", type=Path)
    parser.add_argument("--selected-precursors", type=Path,
                        help="Frozen simulator-only precursor selection TSV from select_high_signal_precursors.py")
    parser.add_argument("--max-precursors", type=int, default=500)
    parser.add_argument("--require-exact", action="store_true",
                        help="Fail unless exactly --max-precursors suitable groups are exported")
    parser.add_argument("--transitions-per-precursor", type=int, default=8)
    parser.add_argument("--minimum-fragments", type=int, default=6)
    parser.add_argument("--min-product-mz", type=float, default=350.0)
    parser.add_argument("--max-product-mz", type=float, default=2000.0)
    parser.add_argument("--mz-window-margin", type=float, default=1.0)
    parser.add_argument("--im-window-margin", type=float, default=0.005)
    parser.add_argument("--gradient-length", type=float, default=30.0)
    parser.add_argument("--rt-mode", choices=("percent", "seconds"), default="percent")
    args = parser.parse_args()

    if args.max_precursors < 1 or args.transitions_per_precursor < 1:
        parser.error("precursor and transition counts must be positive")
    if args.transitions_per_precursor < args.minimum_fragments:
        parser.error("transitions-per-precursor must be >= minimum-fragments")
    if args.gradient_length <= 0:
        parser.error("gradient-length must be positive")
    if not args.db.is_file():
        raise SystemExit(f"Database does not exist: {args.db}")

    truth_out = args.truth_out or args.out.with_name(args.out.stem + ".ground_truth.tsv")
    with sqlite3.connect(args.db) as connection:
        windows = load_windows(connection)
        candidates = read_candidates(
            connection,
            min_product_mz=args.min_product_mz,
            max_product_mz=args.max_product_mz,
            minimum_fragments=args.minimum_fragments,
        )

    filtered = [
        candidate for candidate in candidates
        if inside_window(candidate, windows, args.mz_window_margin, args.im_window_margin)
    ]
    unique_sequences = len({candidate.sequence for candidate in filtered})
    print(
        f"Candidate summary: {len(candidates)} ions with sufficient fragments; "
        f"{len(filtered)} ions / {unique_sequences} unique sequences inside DIA windows."
    )
    if args.selected_precursors:
        selected_keys = load_selected_precursors(args.selected_precursors)
        candidate_by_key = {
            (candidate.sequence, candidate.precursor_charge): candidate
            for candidate in filtered
        }
        missing_keys = [key for key in selected_keys if key not in candidate_by_key]
        if missing_keys:
            raise SystemExit(
                f"Frozen selection contains {len(missing_keys)} precursor groups that no longer pass "
                f"the export fragment/window checks; examples: {missing_keys[:10]}"
            )
        chosen = [candidate_by_key[key] for key in selected_keys]
        if args.require_exact and len(chosen) != args.max_precursors:
            raise SystemExit(
                f"Frozen selection contains {len(chosen)} precursor groups; "
                f"--max-precursors requires exactly {args.max_precursors}."
            )
        print(
            f"Using frozen simulator-only selection: {len(chosen)} precursor groups "
            f"from {args.selected_precursors}"
        )
    else:
        chosen = choose_candidates(filtered, args.max_precursors)
        if not chosen:
            raise SystemExit(
                "No suitable precursor was found. Inspect synthetic_data.db, reduce window margins, "
                "or increase --simulated-peptides when running generate_all.sh."
            )
        if len(chosen) < args.max_precursors:
            message = (
                f"Requested {args.max_precursors} precursor groups but only {len(chosen)} suitable "
                "unique peptide precursors were available. Increase --simulated-peptides (recommended), "
                "increase the generated FASTA pool, or reduce the m/z/IM window margins."
            )
            if args.require_exact:
                raise SystemExit(message)
            print(f"WARNING: {message}", file=sys.stderr)

    write_outputs(
        chosen,
        args.out,
        truth_out,
        args.transitions_per_precursor,
        args.gradient_length,
        args.rt_mode,
    )
    print(f"Wrote {args.out}")
    print(f"Wrote {truth_out}")
    print(f"Selected {len(chosen)} precursors and up to {args.transitions_per_precursor} transitions each.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
