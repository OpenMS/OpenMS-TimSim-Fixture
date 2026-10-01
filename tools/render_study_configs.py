#!/usr/bin/env python3
"""Render TimSim TOML files for a blueprint plus arbitrary study runs."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def toml_escape(value: str | Path) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


def render_template(template: Path, destination: Path, replacements: dict[str, str]) -> None:
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
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--rendered-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--study-manifest", type=Path)
    parser.add_argument("--blueprint-name", default="OpenSwathTimSim_blueprint")
    parser.add_argument("--blueprint-sample-seed", type=int, required=True)
    parser.add_argument("--n-proteins", type=int, required=True)
    parser.add_argument("--num-peptides-total", type=int, required=True)
    parser.add_argument("--num-sample-peptides", type=int, required=True)
    parser.add_argument("--gradient-length", type=float, default=30.0)
    parser.add_argument("--timsim-threads", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--frame-batch-size", type=int, default=100)
    parser.add_argument("--use-gpu", action="store_true", help="Render TimSim configs with performance.use_gpu=true")
    parser.add_argument("--blueprint-only", action="store_true")
    args = parser.parse_args()

    root = args.repo_root.resolve()
    output_root = args.output_root.resolve()
    rendered = args.rendered_dir.resolve()
    reference = args.reference.resolve()
    fasta = args.fasta.resolve()
    modifications = (root / "data" / "unmodified_modifications.toml").resolve()
    blueprint_template = root / "configs" / "study_blueprint.toml.in"
    existing_template = root / "configs" / "study_from_existing.toml.in"
    for path in (blueprint_template, existing_template, modifications, fasta):
        if not path.exists():
            raise SystemExit(f"Required input does not exist: {path}")

    rendered.mkdir(parents=True, exist_ok=True)
    common = {
        "@OUTPUT_ROOT@": toml_escape(output_root),
        "@REFERENCE_D@": toml_escape(reference),
        "@FASTA@": toml_escape(fasta),
        "@MODIFICATIONS@": toml_escape(modifications),
        "@GRADIENT_LENGTH@": f"{args.gradient_length:.6f}",
        "@TIMSIM_THREADS@": str(args.timsim_threads),
        "@BATCH_SIZE@": str(args.batch_size),
        "@FRAME_BATCH_SIZE@": str(args.frame_batch_size),
        "@USE_GPU@": "true" if args.use_gpu else "false",
    }
    blueprint_path = rendered / "000_blueprint.toml"
    render_template(
        blueprint_template,
        blueprint_path,
        {
            **common,
            "@EXPERIMENT_NAME@": args.blueprint_name,
            "@SAMPLE_SEED@": str(args.blueprint_sample_seed),
            "@N_PROTEINS@": str(args.n_proteins),
            "@NUM_PEPTIDES_TOTAL@": str(args.num_peptides_total),
            "@NUM_SAMPLE_PEPTIDES@": str(args.num_sample_peptides),
        },
    )
    print(blueprint_path)

    if args.blueprint_only:
        return 0
    if args.study_manifest is None or not args.study_manifest.is_file():
        raise SystemExit("--study-manifest is required unless --blueprint-only is used")

    with args.study_manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise SystemExit("Study manifest contains no runs")
    for row in rows:
        ordinal = int(row["RunOrdinal"])
        run_name = row["RunName"]
        existing_path = Path(row["InputDirectory"]).resolve()
        if not (existing_path / "synthetic_data.db").is_file():
            raise SystemExit(f"Missing study source DB: {existing_path / 'synthetic_data.db'}")
        config_path = rendered / f"{ordinal:03d}_{row['RunId']}.toml"
        render_template(
            existing_template,
            config_path,
            {
                **common,
                "@EXPERIMENT_NAME@": run_name,
                "@SAMPLE_SEED@": str(int(row["TimSimSampleSeed"])),
                "@EXISTING_PATH@": toml_escape(existing_path),
            },
        )
        print(config_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
