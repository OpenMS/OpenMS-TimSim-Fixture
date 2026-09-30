#!/usr/bin/env python3
"""Append known-absent entrapment assays to a target OpenSWATH TSV.

Two deterministic entrapment constructions are supported:

``paired_hard``
    Sequence-shuffled isobaric assays retain the paired true target's RT/IM.
    This deliberately places an absent assay on top of a real high-signal target
    and is therefore an adversarial interference stress test.

``independent``
    De-novo absent peptide sequences are sampled from the target amino-acid and
    length distributions. Their precursor m/z and theoretical b/y fragments are
    computed from their own sequences. Charge and library-intensity templates
    preserve the target marginals, while RT/IM coordinates come from a
    *different* true target through a deterministic global one-to-one matching.
    Candidate sequence-derived m/z values and donor IM coordinates are matched
    jointly so every assigned m/z/IM pair remains inside a diaPASEF acquisition
    window.

In both modes the entrapment peptide sequence is absent from the TimSim peptide
DB, theoretical b/y fragment m/z values are computed from the entrapment
sequence, and the TimSim TDF files are never modified. Entrapments stay
``Decoy=0`` so OpenDIA's generated decoys remain a separate internal null.
"""
from __future__ import annotations

import argparse
import csv
import math
import random
import sqlite3
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from export_openswath_tsv import (
    AA_MASS,
    OUTPUT_COLUMNS,
    PROTON,
    WATER,
    Window,
    fragment_mz,
    load_windows,
)

MODES = ("paired_hard", "independent")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError(f"TSV has no rows: {path}")
    return rows


def _write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _load_simulated_sequences(db_path: Path) -> set[str]:
    if not db_path.is_file():
        raise FileNotFoundError(db_path)
    with sqlite3.connect(db_path) as connection:
        columns = [row[1] for row in connection.execute('PRAGMA table_info("peptides")')]
        sequence_col = next((c for c in ("sequence", "peptide") if c in columns), None)
        if sequence_col is None:
            raise RuntimeError(f"Could not find peptide sequence column in {db_path}; columns={columns}")
        return {
            str(row[0]).strip().upper()
            for row in connection.execute(f'SELECT "{sequence_col}" FROM "peptides"')
            if row[0]
        }


def _load_dia_windows(db_path: Path) -> list[Window]:
    with sqlite3.connect(db_path) as connection:
        return load_windows(connection)


def _group_rows(rows: list[dict[str, str]]) -> list[tuple[str, list[dict[str, str]]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    order: list[str] = []
    for row in rows:
        group = str(row["TransitionGroupId"])
        if group not in grouped:
            order.append(group)
        grouped[group].append(row)
    return [(group, grouped[group]) for group in order]


def _normalise_sequence(value: str) -> str:
    return str(value).strip().upper()


def _precursor_mz(sequence: str, charge: int) -> float:
    if charge < 1:
        raise ValueError("Precursor charge must be positive")
    neutral = sum(AA_MASS[aa] for aa in sequence) + WATER
    return (neutral + charge * PROTON) / charge


def _aa_population(grouped: list[tuple[str, list[dict[str, str]]]]) -> list[str]:
    population: list[str] = []
    for _, rows in grouped:
        sequence = _normalise_sequence(rows[0]["PeptideSequence"])
        body = sequence[:-1] if sequence and sequence[-1] in {"K", "R"} else sequence
        population.extend(aa for aa in body if aa in AA_MASS)
    if not population:
        raise RuntimeError("Could not derive amino-acid population from target library")
    return population


def _candidate_independent_sequence(
    source_sequence: str,
    aa_population: list[str],
    rng: random.Random,
) -> str:
    preserve_terminal = source_sequence[-1] in {"K", "R"}
    body_length = len(source_sequence) - (1 if preserve_terminal else 0)
    body = "".join(rng.choice(aa_population) for _ in range(body_length))
    terminal = source_sequence[-1] if preserve_terminal else ""
    return body + terminal


def _candidate_shuffle(sequence: str, rng: random.Random) -> str:
    """Shuffle peptide while preserving the C-terminal tryptic residue when present."""
    if len(sequence) < 4:
        chars = list(sequence)
        rng.shuffle(chars)
        return "".join(chars)
    preserve_terminal = sequence[-1] in {"K", "R"}
    body = list(sequence[:-1] if preserve_terminal else sequence)
    rng.shuffle(body)
    return "".join(body) + (sequence[-1] if preserve_terminal else "")


def _transition_specs(source_rows: list[dict[str, str]]) -> list[tuple[str, int, int, float, float]]:
    specs = []
    for row in source_rows:
        ion_type = str(row["FragmentType"]).strip().lower()
        ordinal = int(row["FragmentSeriesNumber"])
        charge = int(row["ProductCharge"])
        source_mz = float(row["ProductMz"])
        intensity = float(row["LibraryIntensity"])
        if ion_type not in {"b", "y"}:
            raise RuntimeError(f"Unsupported target fragment type for entrapment: {ion_type}")
        specs.append((ion_type, ordinal, charge, source_mz, intensity))
    return specs


def _theoretical_rows(
    sequence: str,
    specs: list[tuple[str, int, int, float, float]],
    *,
    precursor_charge: int,
    min_product_mz: float,
    max_product_mz: float,
) -> list[dict[str, Any]] | None:
    """Construct a matched transition set for a shuffled sequence."""
    pool: dict[tuple[str, int, int], float] = {}
    for ion_type in ("b", "y"):
        for ordinal in range(1, len(sequence)):
            for charge in range(1, min(max(1, precursor_charge), 3) + 1):
                mz = fragment_mz(sequence, ion_type, ordinal, charge)
                if min_product_mz <= mz <= max_product_mz:
                    pool[(ion_type, ordinal, charge)] = mz
    if len(pool) < len(specs):
        return None

    out: list[dict[str, Any]] = []
    used: set[tuple[str, int, int]] = set()
    for ion_type, ordinal, charge, source_mz, intensity in specs:
        desired = (ion_type, ordinal, charge)
        if desired in pool and desired not in used:
            chosen = desired
        else:
            available = [key for key in pool if key not in used]
            if not available:
                return None
            chosen = min(available, key=lambda key: abs(pool[key] - source_mz))
        used.add(chosen)
        out.append(
            {
                "type": chosen[0],
                "ordinal": chosen[1],
                "charge": chosen[2],
                "mz": pool[chosen],
                "source_mz": source_mz,
                "library_intensity": intensity,
            }
        )
    return out


def _collision_count(
    entrapment_fragments: list[dict[str, Any]],
    source_rows: list[dict[str, str]],
    ppm: float,
) -> int:
    """Count unchanged paired annotations after sequence shuffling."""
    count = 0
    for fragment, source_row in zip(entrapment_fragments, source_rows):
        mz = float(fragment["mz"])
        target = float(source_row["ProductMz"])
        if abs(mz - target) <= max(abs(target), 1.0) * ppm * 1e-6:
            count += 1
    return count


def _find_entrapment(
    sequence: str,
    source_rows: list[dict[str, str]],
    forbidden_sequences: set[str],
    *,
    seed: int,
    max_attempts: int,
    min_product_mz: float,
    max_product_mz: float,
    collision_ppm: float,
    max_fragment_collisions: int,
) -> tuple[str, list[dict[str, Any]], int] | None:
    specs = _transition_specs(source_rows)
    best: tuple[str, list[dict[str, Any]], int] | None = None
    rng = random.Random(seed)
    for _ in range(max_attempts):
        candidate = _candidate_shuffle(sequence, rng)
        if candidate == sequence or candidate in forbidden_sequences:
            continue
        fragments = _theoretical_rows(
            candidate,
            specs,
            precursor_charge=int(source_rows[0]["PrecursorCharge"]),
            min_product_mz=min_product_mz,
            max_product_mz=max_product_mz,
        )
        if fragments is None:
            continue
        collisions = _collision_count(fragments, source_rows, collision_ppm)
        if best is None or collisions < best[2]:
            best = (candidate, fragments, collisions)
        if collisions <= max_fragment_collisions:
            return candidate, fragments, collisions
    return best if best is not None and best[2] <= max_fragment_collisions else None


def _independent_candidate_options(
    source_group: str,
    source_rows: list[dict[str, str]],
    forbidden_sequences: set[str],
    aa_population: list[str],
    grouped: list[tuple[str, list[dict[str, str]]]],
    donor_mz_intervals: list[list[tuple[float, float]]],
    *,
    seed: int,
    max_attempts: int,
    min_precursor_mz: float,
    max_precursor_mz: float,
    min_product_mz: float,
    max_product_mz: float,
    collision_ppm: float,
    max_fragment_collisions: int,
    min_rt_separation: float,
    min_im_separation: float,
    desired_donor_options: int = 64,
) -> dict[int, tuple[str, float, list[dict[str, Any]], int]]:
    """Generate independent peptide candidates and map them to compatible RT/IM donors.

    A fixed donor can be impossible for a particular peptide length/charge because
    diaPASEF couples precursor m/z and ion mobility.  Instead of forcing a random
    sequence to fit one preassigned donor, sample sequence-derived precursor m/z
    values and retain every *different* donor whose IM makes that m/z fall inside
    a real acquisition window.  A global one-to-one matching is performed later.

    The returned mapping stores one valid sequence candidate per donor index.
    """
    source_sequence = _normalise_sequence(source_rows[0]["PeptideSequence"])
    source_charge = int(source_rows[0]["PrecursorCharge"])
    specs = _transition_specs(source_rows)
    rng = random.Random(seed + 104729)
    options: dict[int, tuple[str, float, list[dict[str, Any]], int]] = {}

    # A local set prevents repeated random candidates from wasting work.  The
    # caller still checks global uniqueness after the one-to-one donor match.
    local_sequences: set[str] = set()
    min_attempts_before_early_stop = min(50, max_attempts)

    for attempt in range(max_attempts):
        candidate = _candidate_independent_sequence(source_sequence, aa_population, rng)
        if (
            candidate == source_sequence
            or candidate in forbidden_sequences
            or candidate in local_sequences
        ):
            continue
        local_sequences.add(candidate)

        precursor_mz = _precursor_mz(candidate, source_charge)
        if not (min_precursor_mz <= precursor_mz <= max_precursor_mz):
            continue

        fragments = _theoretical_rows(
            candidate,
            specs,
            precursor_charge=source_charge,
            min_product_mz=min_product_mz,
            max_product_mz=max_product_mz,
        )
        if fragments is None:
            continue
        collisions = _collision_count(fragments, source_rows, collision_ppm)
        if collisions > max_fragment_collisions:
            continue

        for donor_index, (donor_group, donor_rows) in enumerate(grouped):
            if donor_index in options:
                continue
            if not _coordinate_separated(
                source_group,
                source_rows,
                donor_group,
                donor_rows,
                min_rt_separation=min_rt_separation,
                min_im_separation=min_im_separation,
            ):
                continue
            if any(lower < precursor_mz < upper for lower, upper in donor_mz_intervals[donor_index]):
                options[donor_index] = (candidate, precursor_mz, fragments, collisions)

        if (
            attempt + 1 >= min_attempts_before_early_stop
            and len(options) >= min(desired_donor_options, len(grouped))
        ):
            break

    return options


@dataclass(frozen=True)
class GeneratedEntrapment:
    source_group: str
    source_rows: list[dict[str, str]]
    source_sequence: str
    source_charge: int
    sequence: str
    precursor_mz: float
    fragments: list[dict[str, Any]]
    collision_count: int


def _inside_window(
    precursor_mz: float,
    precursor_im: float,
    windows: list[Window],
    mz_margin: float,
    im_margin: float,
) -> bool:
    if not windows:
        return True
    return any(
        window.mz_lower + mz_margin < precursor_mz < window.mz_upper - mz_margin
        and window.im_lower + im_margin < precursor_im < window.im_upper - im_margin
        for window in windows
    )


def _coordinate_separated(
    source_group: str,
    source_rows: list[dict[str, str]],
    donor_group: str,
    donor_rows: list[dict[str, str]],
    *,
    min_rt_separation: float,
    min_im_separation: float,
) -> bool:
    if donor_group == source_group:
        return False
    source = source_rows[0]
    donor = donor_rows[0]
    source_rt = float(source["NormalizedRetentionTime"])
    source_im = float(source["PrecursorIonMobility"])
    donor_rt = float(donor["NormalizedRetentionTime"])
    donor_im = float(donor["PrecursorIonMobility"])
    # Avoid recreating the paired/hard geometry. A donor is acceptable when it
    # differs materially in at least one coordinate. Acquisition compatibility
    # is enforced later while generating the independent precursor mass.
    return not (
        abs(donor_rt - source_rt) < min_rt_separation
        and abs(donor_im - source_im) < min_im_separation
    )


def _mz_intervals_for_im(
    precursor_im: float,
    windows: list[Window],
    mz_margin: float,
    im_margin: float,
) -> list[tuple[float, float]]:
    """Return diaPASEF precursor-m/z intervals usable at one IM coordinate."""
    if not windows:
        return [(-math.inf, math.inf)]
    intervals: list[tuple[float, float]] = []
    for window in windows:
        if window.im_lower + im_margin < precursor_im < window.im_upper - im_margin:
            lower = window.mz_lower + mz_margin
            upper = window.mz_upper - mz_margin
            if lower < upper:
                intervals.append((lower, upper))
    return intervals


def _plan_independent_entrapments(
    sources: list[tuple[str, list[dict[str, str]]]],
    grouped: list[tuple[str, list[dict[str, str]]]],
    forbidden_sequences: set[str],
    aa_population: list[str],
    windows: list[Window],
    *,
    seed: int,
    max_attempts: int,
    min_precursor_mz: float,
    max_precursor_mz: float,
    min_product_mz: float,
    max_product_mz: float,
    collision_ppm: float,
    max_fragment_collisions: int,
    min_rt_separation: float,
    min_im_separation: float,
    mz_window_margin: float,
    im_window_margin: float,
) -> tuple[list[GeneratedEntrapment], list[list[dict[str, str]]]]:
    """Jointly choose independent sequences and one-to-one RT/IM donors.

    This is a bipartite matching problem because diaPASEF acquisition geometry
    couples an independently generated precursor m/z to the donor IM.  Each
    source template first generates a pool of valid sequence candidates.  Every
    candidate contributes edges to all compatible RT/IM donors.  We then solve
    one global one-to-one source->donor matching and use the sequence candidate
    associated with each matched edge.

    For the default 500 sources and 500 donors, a complete matching uses every
    donor exactly once, preserving the true-target RT/IM marginal distribution.
    """
    if len(sources) > len(grouped):
        raise RuntimeError(
            "Independent mode requires --count <= number of true target groups so RT/IM donors can be one-to-one"
        )

    donor_mz_intervals = [
        _mz_intervals_for_im(
            float(rows[0]["PrecursorIonMobility"]),
            windows,
            mz_window_margin,
            im_window_margin,
        )
        for _, rows in grouped
    ]

    option_maps: list[dict[int, tuple[str, float, list[dict[str, Any]], int]]] = []
    # Reserve every candidate that contributes a donor edge so two source
    # templates cannot independently plan the same absent peptide sequence.
    planning_forbidden = set(forbidden_sequences)
    for source_index, (source_group, source_rows) in enumerate(sources, start=1):
        options = _independent_candidate_options(
            source_group,
            source_rows,
            planning_forbidden,
            aa_population,
            grouped,
            donor_mz_intervals,
            seed=seed + 1009 * source_index,
            max_attempts=max_attempts,
            min_precursor_mz=min_precursor_mz,
            max_precursor_mz=max_precursor_mz,
            min_product_mz=min_product_mz,
            max_product_mz=max_product_mz,
            collision_ppm=collision_ppm,
            max_fragment_collisions=max_fragment_collisions,
            min_rt_separation=min_rt_separation,
            min_im_separation=min_im_separation,
        )
        if options:
            planning_forbidden.update(candidate[0] for candidate in options.values())
        if not options:
            source_sequence = _normalise_sequence(source_rows[0]["PeptideSequence"])
            source_charge = int(source_rows[0]["PrecursorCharge"])
            raise RuntimeError(
                "Could not generate any independent peptide/donor combination for "
                f"source={source_group} ({source_sequence}/{source_charge}+) after {max_attempts} attempts. "
                "This is now a joint sequence+geometry search; increasing --max-attempts is appropriate only if this persists."
            )
        option_maps.append(options)

    # Match the most constrained sources first.  The augmenting-path matcher can
    # reassign previous records, so a complete solution is found whenever one
    # exists in the generated candidate graph.
    record_order = sorted(range(len(sources)), key=lambda i: len(option_maps[i]))
    donor_to_record: dict[int, int] = {}

    def augment(record_index: int, seen: set[int]) -> bool:
        # Deterministic ordering with a seed-derived tie break keeps fixture
        # generation reproducible while avoiding systematic low-index donors.
        donor_indices = list(option_maps[record_index])
        rng = random.Random(seed + 15485863 * (record_index + 1))
        rng.shuffle(donor_indices)
        for donor_index in donor_indices:
            if donor_index in seen:
                continue
            seen.add(donor_index)
            previous = donor_to_record.get(donor_index)
            if previous is None or augment(previous, seen):
                donor_to_record[donor_index] = record_index
                return True
        return False

    for record_index in record_order:
        if not augment(record_index, set()):
            constrained = sorted(
                ((len(option_maps[i]), sources[i][0]) for i in range(len(sources))),
                key=lambda item: item[0],
            )[:10]
            raise RuntimeError(
                "Could not construct a complete one-to-one independent sequence/RT/IM donor matching. "
                f"Most constrained source option counts={constrained}. Increase --max-attempts before relaxing separation or DIA margins."
            )

    record_to_donor = {record: donor for donor, record in donor_to_record.items()}
    if len(record_to_donor) != len(sources):
        raise AssertionError("Internal error: incomplete independent donor assignment")

    generated: list[GeneratedEntrapment] = []
    coordinate_donors: list[list[dict[str, str]]] = []
    selected_sequences: set[str] = set()
    for record_index, (source_group, source_rows) in enumerate(sources):
        donor_index = record_to_donor[record_index]
        candidate, precursor_mz, fragments, collisions = option_maps[record_index][donor_index]
        if candidate in forbidden_sequences or candidate in selected_sequences:
            # Collisions between independently planned sources are extraordinarily
            # unlikely.  Fail explicitly rather than silently duplicate an assay.
            raise RuntimeError(
                f"Independent planner selected a duplicate/forbidden sequence: {candidate}"
            )
        selected_sequences.add(candidate)
        source_sequence = _normalise_sequence(source_rows[0]["PeptideSequence"])
        source_charge = int(source_rows[0]["PrecursorCharge"])
        generated.append(
            GeneratedEntrapment(
                source_group=source_group,
                source_rows=source_rows,
                source_sequence=source_sequence,
                source_charge=source_charge,
                sequence=candidate,
                precursor_mz=precursor_mz,
                fragments=fragments,
                collision_count=collisions,
            )
        )
        coordinate_donors.append(grouped[donor_index][1])

    forbidden_sequences.update(selected_sequences)
    return generated, coordinate_donors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-transitions", type=Path, required=True)
    parser.add_argument("--baseline-db", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True, help="Combined true-target + entrapment transition TSV")
    parser.add_argument("--truth-out", type=Path, required=True, help="Entrapment ground-truth TSV")
    parser.add_argument("--mode", choices=MODES, default="paired_hard")
    parser.add_argument("--count", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--max-attempts", type=int, default=1000)
    parser.add_argument("--min-product-mz", type=float, default=350.0)
    parser.add_argument("--max-product-mz", type=float, default=2000.0)
    parser.add_argument(
        "--fragment-collision-ppm",
        type=float,
        default=25.0,
        help="Tolerance used to count entrapment fragments colliding with the paired true target.",
    )
    parser.add_argument(
        "--max-paired-fragment-collisions",
        type=int,
        default=1,
        help="Maximum allowed entrapment fragments within the collision tolerance of a paired true-target fragment.",
    )
    parser.add_argument(
        "--independent-min-rt-separation",
        type=float,
        default=2.0,
        help="Minimum assay-RT separation unless the IM separation criterion is met (NormalizedRetentionTime units).",
    )
    parser.add_argument(
        "--independent-min-im-separation",
        type=float,
        default=0.03,
        help="Minimum IM separation unless the RT separation criterion is met (1/K0).",
    )
    parser.add_argument("--mz-window-margin", type=float, default=1.0)
    parser.add_argument("--im-window-margin", type=float, default=0.005)
    args = parser.parse_args()

    if args.count < 0:
        parser.error("--count must be non-negative")
    if args.max_attempts < 1:
        parser.error("--max-attempts must be positive")
    if args.max_paired_fragment_collisions < 0:
        parser.error("--max-paired-fragment-collisions must be non-negative")
    if args.independent_min_rt_separation < 0 or args.independent_min_im_separation < 0:
        parser.error("Independent coordinate separations must be non-negative")

    target_rows = _read_tsv(args.target_transitions)
    missing = set(OUTPUT_COLUMNS) - set(target_rows[0])
    if missing:
        raise SystemExit(f"Target transition TSV is missing columns: {sorted(missing)}")
    if any(str(row.get("Decoy", "0")).strip().lower() not in {"0", "false"} for row in target_rows):
        raise SystemExit("--target-transitions must contain target rows only (Decoy=0)")

    grouped = _group_rows(target_rows)
    if args.count and not grouped:
        raise SystemExit("No target precursor groups were found")
    if args.mode == "independent" and args.count > len(grouped):
        raise SystemExit(
            f"Independent mode requested {args.count} entrapments but only {len(grouped)} true target groups are available"
        )

    simulated_sequences = _load_simulated_sequences(args.baseline_db)
    target_sequences = {_normalise_sequence(rows[0]["PeptideSequence"]) for _, rows in grouped}
    forbidden = set(simulated_sequences) | set(target_sequences)
    amino_acid_population = _aa_population(grouped)
    target_precursor_mz = [float(rows[0]["PrecursorMz"]) for _, rows in grouped]
    min_precursor_mz = min(target_precursor_mz)
    max_precursor_mz = max(target_precursor_mz)

    windows = _load_dia_windows(args.baseline_db)
    generated: list[GeneratedEntrapment] = []

    if args.mode == "paired_hard":
        source_index = 0
        cycle = 0
        while len(generated) < args.count:
            if source_index >= len(grouped):
                source_index = 0
                cycle += 1
            source_group, source_rows = grouped[source_index]
            source_index += 1
            source_sequence = _normalise_sequence(source_rows[0]["PeptideSequence"])
            source_charge = int(source_rows[0]["PrecursorCharge"])
            local_seed = args.seed + 1000003 * cycle + 1009 * source_index
            paired = _find_entrapment(
                source_sequence,
                source_rows,
                forbidden,
                seed=local_seed,
                max_attempts=args.max_attempts,
                min_product_mz=args.min_product_mz,
                max_product_mz=args.max_product_mz,
                collision_ppm=args.fragment_collision_ppm,
                max_fragment_collisions=args.max_paired_fragment_collisions,
            )
            if paired is None:
                if cycle > max(5, math.ceil(args.count / max(len(grouped), 1)) + 5):
                    raise SystemExit(
                        f"Could generate only {len(generated)} clean paired-hard entrapments after repeated passes. "
                        "Increase --max-attempts or relax --max-paired-fragment-collisions explicitly."
                    )
                continue
            entrapment_sequence, fragments, collision_count = paired
            entrapment_precursor_mz = float(source_rows[0]["PrecursorMz"])
            if entrapment_sequence in forbidden:
                raise AssertionError("Internal error: duplicate/forbidden entrapment sequence")
            forbidden.add(entrapment_sequence)
            generated.append(
                GeneratedEntrapment(
                    source_group=source_group,
                    source_rows=source_rows,
                    source_sequence=source_sequence,
                    source_charge=source_charge,
                    sequence=entrapment_sequence,
                    precursor_mz=entrapment_precursor_mz,
                    fragments=fragments,
                    collision_count=collision_count,
                )
            )
        coordinate_donors = [item.source_rows for item in generated]
    else:
        # Independent sequence generation and RT/IM assignment are solved jointly.
        # A random peptide mass is not forced onto a preassigned IM donor because
        # diaPASEF m/z and IM windows are coupled; instead, candidate sequence masses
        # define the set of compatible donors and a global one-to-one matching picks
        # a valid donor for every entrapment.
        independent_sources = grouped[: args.count]
        generated, coordinate_donors = _plan_independent_entrapments(
            independent_sources,
            grouped,
            forbidden,
            amino_acid_population,
            windows,
            seed=args.seed,
            max_attempts=args.max_attempts,
            min_precursor_mz=min_precursor_mz,
            max_precursor_mz=max_precursor_mz,
            min_product_mz=args.min_product_mz,
            max_product_mz=args.max_product_mz,
            collision_ppm=args.fragment_collision_ppm,
            max_fragment_collisions=args.max_paired_fragment_collisions,
            min_rt_separation=args.independent_min_rt_separation,
            min_im_separation=args.independent_min_im_separation,
            mz_window_margin=args.mz_window_margin,
            im_window_margin=args.im_window_margin,
        )

    entrapment_rows: list[dict[str, Any]] = []
    truth_rows: list[dict[str, Any]] = []
    for rank, (item, donor_rows) in enumerate(zip(generated, coordinate_donors), start=1):
        source_first = item.source_rows[0]
        donor_first = donor_rows[0]
        donor_group = donor_first["TransitionGroupId"]
        donor_sequence = _normalise_sequence(donor_first["PeptideSequence"])
        donor_charge = int(donor_first["PrecursorCharge"])
        assigned_rt = float(donor_first["NormalizedRetentionTime"])
        assigned_im = float(donor_first["PrecursorIonMobility"])
        source_rt = float(source_first["NormalizedRetentionTime"])
        source_im = float(source_first["PrecursorIonMobility"])

        group_id = f"ENTRAPMENT_{args.mode.upper()}_{rank:04d}_{item.sequence}_{item.source_charge}"
        protein_id = f"ENTRAPMENT_{args.mode.upper()}_PROTEIN_{rank:04d}"
        label = f"ENTRAPMENT_{args.mode.upper()}_{item.sequence}"

        for source_row, fragment in zip(item.source_rows, item.fragments):
            row = dict(source_row)
            annotation = f"{fragment['type']}{fragment['ordinal']}^{fragment['charge']}"
            row.update(
                {
                    "PrecursorMz": f"{item.precursor_mz:.8f}",
                    "ProductMz": f"{float(fragment['mz']):.8f}",
                    "LibraryIntensity": f"{float(fragment['library_intensity']):.6f}",
                    "NormalizedRetentionTime": f"{assigned_rt:.8f}",
                    "PrecursorIonMobility": f"{assigned_im:.8f}",
                    "ProteinId": protein_id,
                    "PeptideSequence": item.sequence,
                    "ModifiedPeptideSequence": item.sequence,
                    "PeptideGroupLabel": label,
                    "ProductCharge": int(fragment["charge"]),
                    "FragmentType": fragment["type"],
                    "FragmentSeriesNumber": int(fragment["ordinal"]),
                    "TransitionGroupId": group_id,
                    "TransitionId": f"{group_id}_{annotation}",
                    "Decoy": 0,
                    "DetectingTransition": 1,
                    "IdentifyingTransition": 0,
                    "QuantifyingTransition": 1,
                }
            )
            entrapment_rows.append(row)

        truth_rows.append(
            {
                "EntrapmentRank": rank,
                "EntrapmentMode": args.mode,
                "EntrapmentTransitionGroupId": group_id,
                "EntrapmentSequence": item.sequence,
                "PrecursorCharge": item.source_charge,
                "PrecursorMz": f"{item.precursor_mz:.8f}",
                "NormalizedRetentionTime": f"{assigned_rt:.8f}",
                "PrecursorIonMobility": f"{assigned_im:.8f}",
                "ProteinId": protein_id,
                "SourceTargetTransitionGroupId": item.source_group,
                "SourceTargetSequence": item.source_sequence,
                "SourceTargetCharge": item.source_charge,
                "SourceTargetPrecursorMz": source_first["PrecursorMz"],
                "SourceTargetNormalizedRetentionTime": source_first["NormalizedRetentionTime"],
                "SourceTargetPrecursorIonMobility": source_first["PrecursorIonMobility"],
                "CoordinateDonorTransitionGroupId": donor_group,
                "CoordinateDonorSequence": donor_sequence,
                "CoordinateDonorCharge": donor_charge,
                "CoordinateDonorNormalizedRetentionTime": donor_first["NormalizedRetentionTime"],
                "CoordinateDonorPrecursorIonMobility": donor_first["PrecursorIonMobility"],
                "CoordinateRtSeparation": abs(assigned_rt - source_rt),
                "CoordinateImSeparation": abs(assigned_im - source_im),
                "IndependentCoordinatePasefCompatible": _inside_window(
                    item.precursor_mz,
                    assigned_im,
                    windows,
                    args.mz_window_margin,
                    args.im_window_margin,
                ),
                "NumberOfTransitions": len(item.source_rows),
                "PairedFragmentCollisions": item.collision_count,
                "FragmentCollisionPpm": args.fragment_collision_ppm,
                "SequenceAbsentFromTimSim": True,
            }
        )

    combined = [*target_rows, *entrapment_rows]
    transition_ids = [str(row["TransitionId"]) for row in combined]
    if len(transition_ids) != len(set(transition_ids)):
        raise SystemExit("Combined target/entrapment TransitionId values are not unique")

    _write_tsv(args.out, combined, OUTPUT_COLUMNS)
    truth_columns = [
        "EntrapmentRank",
        "EntrapmentMode",
        "EntrapmentTransitionGroupId",
        "EntrapmentSequence",
        "PrecursorCharge",
        "PrecursorMz",
        "NormalizedRetentionTime",
        "PrecursorIonMobility",
        "ProteinId",
        "SourceTargetTransitionGroupId",
        "SourceTargetSequence",
        "SourceTargetCharge",
        "SourceTargetPrecursorMz",
        "SourceTargetNormalizedRetentionTime",
        "SourceTargetPrecursorIonMobility",
        "CoordinateDonorTransitionGroupId",
        "CoordinateDonorSequence",
        "CoordinateDonorCharge",
        "CoordinateDonorNormalizedRetentionTime",
        "CoordinateDonorPrecursorIonMobility",
        "CoordinateRtSeparation",
        "CoordinateImSeparation",
        "IndependentCoordinatePasefCompatible",
        "NumberOfTransitions",
        "PairedFragmentCollisions",
        "FragmentCollisionPpm",
        "SequenceAbsentFromTimSim",
    ]
    _write_tsv(args.truth_out, truth_rows, truth_columns)

    print(
        f"Entrapment library OK ({args.mode}): {len(grouped)} true target groups + "
        f"{len(truth_rows)} known-absent entrapment groups; wrote {len(combined)} transitions"
    )
    if args.mode == "independent" and truth_rows:
        rt_sep = [float(row["CoordinateRtSeparation"]) for row in truth_rows]
        im_sep = [float(row["CoordinateImSeparation"]) for row in truth_rows]
        print(
            "Independent coordinate assignment: "
            f"median |delta RT|={sorted(rt_sep)[len(rt_sep)//2]:.4g}, "
            f"median |delta IM|={sorted(im_sep)[len(im_sep)//2]:.4g}; "
            "all assigned m/z/IM coordinates pass diaPASEF window compatibility"
        )
    print(f"Combined library:  {args.out}")
    print(f"Entrapment truth:  {args.truth_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
