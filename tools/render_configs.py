#!/usr/bin/env python3
"""Render TimSim TOML templates with absolute paths and fixture parameters."""
from __future__ import annotations

import argparse
from pathlib import Path


def toml_string(path: Path) -> str:
    return str(path.resolve()).replace("\\", "\\\\").replace('"', '\\"')


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--rendered-dir", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--replicate-existing-path", type=Path, required=True)
    parser.add_argument("--treatment-existing-path", type=Path, required=True)
    parser.add_argument("--n-proteins", type=int, required=True)
    parser.add_argument("--num-peptides-total", type=int, required=True)
    parser.add_argument("--num-sample-peptides", type=int, required=True)
    parser.add_argument("--gradient-length", type=float, default=30.0)
    parser.add_argument("--sample-seed", type=int, default=41)
    parser.add_argument("--timsim-threads", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--frame-batch-size", type=int, default=100)
    args = parser.parse_args()

    root = args.root.resolve()
    reference = args.reference.resolve()
    output_root = args.output_root.resolve()
    rendered = args.rendered_dir.resolve()
    fasta = args.fasta.resolve()

    if not reference.is_dir():
        raise SystemExit(f"Reference .d directory does not exist: {reference}")
    for required in ("analysis.tdf", "analysis.tdf_bin"):
        if not (reference / required).is_file():
            raise SystemExit(f"Reference is missing {required}: {reference}")
    if not fasta.is_file():
        raise SystemExit(f"Generated FASTA does not exist: {fasta}")

    modifications = root / "data" / "unmodified_modifications.toml"
    rendered.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)

    replacements = {
        "@OUTPUT_ROOT@": toml_string(output_root),
        "@REFERENCE_D@": toml_string(reference),
        "@FASTA@": toml_string(fasta),
        "@MODIFICATIONS@": toml_string(modifications),
        "@REPLICATE_EXISTING_PATH@": toml_string(args.replicate_existing_path),
        "@TREATMENT_EXISTING_PATH@": toml_string(args.treatment_existing_path),
        "@N_PROTEINS@": str(args.n_proteins),
        "@NUM_PEPTIDES_TOTAL@": str(args.num_peptides_total),
        "@NUM_SAMPLE_PEPTIDES@": str(args.num_sample_peptides),
        "@GRADIENT_LENGTH@": f"{args.gradient_length:.6f}",
        "@SAMPLE_SEED@": str(args.sample_seed),
        "@TIMSIM_THREADS@": str(args.timsim_threads),
        "@BATCH_SIZE@": str(args.batch_size),
        "@FRAME_BATCH_SIZE@": str(args.frame_batch_size),
    }

    # Render only the supported experiment templates and clear any stale rendered configs.
    template_names = (
        "run_01_baseline.toml.in",
        "run_02_biological_replicate.toml.in",
        "run_03_treatment.toml.in",
    )

    rendered.mkdir(parents=True, exist_ok=True)
    for stale in rendered.glob("*.toml"):
        stale.unlink()

    for template_name in template_names:
        template = root / "configs" / template_name
        if not template.is_file():
            raise SystemExit(f"Required config template is missing: {template}")
        text = template.read_text(encoding="utf-8")
        for token, value in replacements.items():
            text = text.replace(token, value)
        unresolved = sorted({token for token in text.split() if token.startswith("@") and token.endswith("@")})
        if unresolved:
            raise SystemExit(f"Unresolved placeholders in {template}: {unresolved}")
        destination = rendered / template.name.removesuffix(".in")
        destination.write_text(text, encoding="utf-8")
        print(destination)

    # Validate the rendered configs before the expensive simulation starts.
    try:
        import tomllib
    except ModuleNotFoundError:  # pragma: no cover - Python 3.11+ is required
        raise SystemExit("Python 3.11 or newer is required")

    for config_path in sorted(rendered.glob("*.toml")):
        with config_path.open("rb") as handle:
            config = tomllib.load(handle)
        flat = {}
        for value in config.values():
            if isinstance(value, dict):
                flat.update(value)
        required = ("save_path", "reference_path", "fasta_path")
        missing = [key for key in required if not str(flat.get(key, "")).strip()]
        if flat.get("from_existing"):
            existing_value = str(flat.get("existing_path", "")).strip()
            if not existing_value:
                missing.append("existing_path")
            elif Path(existing_value).suffix == ".db":
                raise SystemExit(
                    f"Rendered TimSim config {config_path} sets existing_path to a database file. "
                    "Use the directory containing synthetic_data.db."
                )
        if missing:
            raise SystemExit(
                f"Rendered TimSim config {config_path} is missing required options: "
                + ", ".join(missing)
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
