#!/usr/bin/env python3
"""Generate a deterministic synthetic tryptic proteome for TimSim fixtures."""
from __future__ import annotations

import argparse
import random
from pathlib import Path

# No internal K/R so every designed peptide is a clean tryptic product.
BODY_ALPHABET = "ACDEFGHILMNPQSTVWY"
FIRST_ALPHABET = BODY_ALPHABET.replace("P", "")


def make_peptide(rng: random.Random, index: int, used: set[str]) -> str:
    """Create one unique, canonical, fully tryptic peptide."""
    for _ in range(10000):
        total_length = rng.randint(9, 18)
        body_length = total_length - 1
        body = rng.choice(FIRST_ALPHABET) + "".join(
            rng.choice(BODY_ALPHABET) for _ in range(body_length - 1)
        )
        terminal = "K" if index % 2 == 0 else "R"
        peptide = body + terminal
        if peptide not in used:
            used.add(peptide)
            return peptide
    raise RuntimeError("Unable to generate another unique peptide")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--peptides", type=int, required=True)
    parser.add_argument("--peptides-per-protein", type=int, default=20)
    parser.add_argument("--seed", type=int, default=1729)
    args = parser.parse_args()

    if args.peptides < 1:
        parser.error("--peptides must be positive")
    if args.peptides_per_protein < 1:
        parser.error("--peptides-per-protein must be positive")

    rng = random.Random(args.seed)
    used: set[str] = set()
    peptides = [make_peptide(rng, index, used) for index in range(args.peptides)]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as handle:
        for protein_index, start in enumerate(
            range(0, len(peptides), args.peptides_per_protein), start=1
        ):
            chunk = peptides[start:start + args.peptides_per_protein]
            handle.write(f">TIMSIM_PROTEIN_{protein_index:05d}\n")
            sequence = "".join(chunk)
            for line_start in range(0, len(sequence), 80):
                handle.write(sequence[line_start:line_start + 80] + "\n")

    protein_count = (len(peptides) + args.peptides_per_protein - 1) // args.peptides_per_protein
    print(f"Wrote {args.out}")
    print(f"Designed peptides: {len(peptides)}")
    print(f"Synthetic proteins: {protein_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
