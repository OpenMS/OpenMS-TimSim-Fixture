# Reproducibility

## Upstream baseline

The cleaned repository uses `imspy-simulation==0.4.2` directly. It does not build a private rustims connector and it does not modify installed TimSim Python files.

This supersedes the old project's local fixes for deterministic peptide/isotope calculations and the `add_real_data_noise=false` gate.

## Seed roles

The generator records the seed plan under `fixture_manifest.json -> randomization`. The three distinct seed roles are:

- `fasta_seed`: fixes the synthetic peptide/protein universe;
- `timsim_sample_seed`: controls TimSim sampling for a synthetic experiment;
- `condition_seed`: controls replicate and treatment abundance perturbations.

The user-facing `--simulation-seed` deterministically derives the latter two while the FASTA seed remains independently configurable.

Keeping these roles separate allows repeated simulated experiments to use the same synthetic proteome without conflating peptide-universe changes with acquisition/sampling changes.

## Fixture manifest

`fixture_manifest.json` records:

- fixture kind and run names;
- Python, `imspy-simulation`, and connector versions;
- reference `.d` path;
- all seed values;
- target/entrapment counts;
- condition effect parameters;
- high-signal selection thresholds;
- key generated artifact names;
- SHA-256 hashes for the synthetic FASTA, target/combined libraries, entrapment truth, and major truth files.

## Replay expectation

For a pinned software environment, identical inputs and seeds should reproduce the scientific fixture truth and key artifact hashes.

The repository intentionally tests this property rather than forcing deterministic behavior with local simulator patches.

A full integration replay should be performed whenever the pinned TimSim version changes:

```text
same reference + same configuration + same seeds
           |
           +--> fixture A
           +--> fixture B

compare manifest parameters and key truth/library SHA-256 hashes
```

The raw `.d` byte stream may be tested separately if exact binary replay is required, but scientific truth replay is the minimum acceptance requirement for this project.

## Reference data

The reference `.d` is required for acquisition geometry. Canonical configs explicitly disable real-data superposition/noise. The reference dataset itself should therefore be tracked by an external immutable identifier or checksum in any published/CI fixture deployment.

## Multi-run study seed plan

Multi-run studies use seed-plan version 3. The synthetic FASTA seed, blueprint TimSim
seed, study biological-variation seed, entrapment seed, and per-run TimSim child seeds
are separate. `OpenSwathTimSim.study_manifest.tsv` records each run's TimSim and
abundance seed, and `fixture_manifest.json` records the study-level seeds and hashes of
all major truth/library artifacts.

Reusing the same seeds and pinned simulator version is expected to reproduce the same
study design and source abundance databases. Reproducibility of generated TDF signal is
validated end-to-end rather than maintained through local simulator patches.
