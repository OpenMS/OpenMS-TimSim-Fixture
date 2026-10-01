#!/usr/bin/env python3
"""Render one TimSim from-existing config for a lazily materialized study run."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def toml_escape(value: str | Path) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def read_run(manifest: Path, ordinal: int) -> dict[str, str]:
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if int(row["RunOrdinal"]) == ordinal]
    if len(rows) != 1:
        raise SystemExit(f"Expected one manifest row for ordinal {ordinal}; found {len(rows)}")
    return rows[0]


def render(template: Path, destination: Path, replacements: dict[str, str]) -> None:
    text = template.read_text(encoding="utf-8")
    for token, value in replacements.items():
        text = text.replace(token, value)
    unresolved = sorted({piece for piece in text.replace("=", " ").split() if piece.startswith("@") and piece.endswith("@")})
    if unresolved:
        raise SystemExit(f"Unresolved placeholders in {template}: {unresolved}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--study-root", type=Path, required=True)
    parser.add_argument("--study-manifest", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--run-ordinal", type=int, required=True)
    parser.add_argument("--existing-path", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--timsim-threads", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--frame-batch-size", type=int, default=100)
    parser.add_argument("--use-gpu", action="store_true")
    args = parser.parse_args()

    for path in (args.study_manifest, args.plan, args.reference, args.fasta):
        if not path.exists():
            raise SystemExit(f"Required input does not exist: {path}")
    if not (args.existing_path / "synthetic_data.db").is_file():
        raise SystemExit(f"Materialized source DB missing: {args.existing_path / 'synthetic_data.db'}")

    row = read_run(args.study_manifest, args.run_ordinal)
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    gradient = float(plan.get("gradient_length_seconds", 30.0))
    root = args.repo_root.resolve()
    template = root / "configs" / "study_from_existing.toml.in"
    modifications = root / "data" / "unmodified_modifications.toml"

    render(
        template,
        args.out,
        {
            "@OUTPUT_ROOT@": toml_escape(args.study_root.resolve()),
            "@REFERENCE_D@": toml_escape(args.reference.resolve()),
            "@FASTA@": toml_escape(args.fasta.resolve()),
            "@MODIFICATIONS@": toml_escape(modifications.resolve()),
            "@GRADIENT_LENGTH@": f"{gradient:.6f}",
            "@TIMSIM_THREADS@": str(args.timsim_threads),
            "@BATCH_SIZE@": str(args.batch_size),
            "@FRAME_BATCH_SIZE@": str(args.frame_batch_size),
            "@USE_GPU@": "true" if args.use_gpu else "false",
            "@EXPERIMENT_NAME@": row["RunName"],
            "@SAMPLE_SEED@": str(int(row["TimSimSampleSeed"])),
            "@EXISTING_PATH@": toml_escape(args.existing_path.resolve()),
        },
    )
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
