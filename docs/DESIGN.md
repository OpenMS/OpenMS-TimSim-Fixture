# OpenMS-TimSim-Fixture design

## Purpose

The repository exists to answer three regression questions about OpenDIA with one frozen synthetic study:

1. **Identification:** does OpenDIA recover precursor groups that TimSim says are present and measurable?
2. **Quantification:** does OpenDIA recover the known relative abundance structure across control and treatment conditions?
3. **Error control:** how many known-absent entrapment targets pass the same statistical filters as the present targets?

The fixture is deliberately a correctness/regression fixture, not an attempt to simulate every source of experimental complexity.

## Scientific separation

The pipeline is intentionally one-way:

```text
synthetic proteome
      |
      v
TimSim simulation + explicit condition truth
      |
      v
simulator-only high-signal target selection
      |
      +----> frozen target truth
      |
      +----> known-absent entrapments
      |
      v
OpenDIA input library
      |
      v
OpenDIA
      |
      v
benchmark against truth fixed before OpenDIA ran
```

No OpenDIA q-value, score, feature, peak, or identification result may influence fixture target selection.

## Three runs

### Baseline

Defines the shared peptide identities, precursor charge states, RT, m/z, ion mobility, and baseline abundances.

### Biological replicate

Starts from the baseline molecular database and applies modest deterministic peptide-level abundance variation. Molecular coordinates are preserved.

### Treatment

Starts from the same baseline molecular database and applies coherent protein-level up/down fold changes plus modest sample-level variation. Molecular coordinates are preserved.

This means differences in measured OpenDIA quantities can be compared to explicit abundance truth without confounding them with a different synthetic peptide population.

## High-signal target selection

The current selector requires each target to be measurable in all three runs and considers:

- realized event proxy;
- retained RT/frame abundance mass;
- retained IM/scan abundance mass;
- selected charge-state abundance fraction;
- minimum usable fragment count;
- safe m/z and IM placement within diaPASEF acquisition geometry;
- exclusion from RT acquisition edges;
- diversity across RT, precursor m/z, IM, charge, protein, and SWATH map.

Only one realized charge state per peptide sequence is retained. Selection fails if the requested target count cannot be satisfied under the fixed criteria.

For multi-run production studies, eligibility is evaluated on the blueprint only. The canonical production contract is **250 proteins × 4 precursors/protein = 1,000 targets**. A protein must have at least four eligible unique peptide precursors. Within each protein, selection is two-stage: the strongest 2× signal shortlist is formed first, then the final four are chosen by simulator-only fragment-collision burden. Protein selection is balanced across RT, precursor m/z, ion mobility, charge, and SWATH geometry.

Collision specificity is computed from the complete TimSim blueprint fragment universe using precursor window group, RT/IM proximity, fragment m/z proximity, and predicted fragment intensity. It does not use observed raw intensity, the raw-oracle sweep result, or any OpenDIA output. The production selector fails closed if fewer than the declared number of proteins satisfy the multiplicity contract; the remedy is to enlarge the synthetic candidate universe, not relax high-signal thresholds.

## Transition library

For each selected precursor, the generator exports a fixed number of OpenSWATH/OpenDIA input transitions from the baseline TimSim truth.

The target-only library is kept separately from the active target+entrapment library so the benchmark can always distinguish the true analytes from known-absent negatives.

## Entrapments

Entrapments are target-like precursor groups that are guaranteed absent from the TimSim peptide population. They remain `Decoy=0` because they are external negative controls, not OpenDIA decoys.

The default `independent` construction is the canonical external null for q-value/FDR calibration. It generates absent peptide sequences with their own precursor m/z and theoretical fragment ions, then assigns target-like RT/IM coordinates through a deterministic one-to-one donor matching while keeping the m/z/IM pair inside the diaPASEF acquisition geometry.

`paired_hard` is retained as a separate adversarial interference stress test. It deliberately keeps the paired true target's precursor m/z, RT, and IM while changing the absent sequence and fragment annotations. Because that null is intentionally colocated with real high-signal analyte evidence, its acceptance fraction must not be interpreted as calibration of OpenDIA's nominal q-values.

The fixture checker verifies that entrapment peptide sequences do not occur in the baseline TimSim database.

## Decoys

The fixture does not provide OpenDIA decoys. OpenDIA should run its normal assay-preparation/decoy-generation path. This keeps external entrapment error estimation independent from the model's own target-decoy mechanism.


## Benchmark lane isolation

OpenDIA scoring is learned from the library/search population. External target-labelled negatives can therefore alter Percolator training and which target peak group is ultimately selected, even when the synthetic `.d` files and target truth are unchanged.

For that reason, the canonical benchmark uses separate OpenDIA runs:

- **target-only lane:** `OpenSwathTimSim.target_only.transitions.tsv`; authoritative for identification, localization, correctness, and quantification;
- **independent-entrapment lane:** `OpenSwathTimSim.transitions.tsv`; authoritative for external-null q-value/FDR calibration;
- **paired-hard lane:** optional adversarial interference challenge, interpreted separately from nominal FDR calibration.

Quantification metrics must not be compared across different entrapment constructions when they were produced from a combined target+entrapment scoring population. The frozen TimSim truth is shared, but each scientific question gets its own OpenDIA scoring lane.

## Noise policy

`add_real_data_noise=false` and `superimpose_on_reference=false` are part of the canonical fixture contract.

The Bruker reference run contributes acquisition geometry/layout. It must not contribute hidden real analyte signal to the correctness fixture.

## Version policy

The fixture pins:

```text
imspy-simulation==0.4.2
```

The pinned upstream simulator is used directly. If reproducibility regresses in a later release, the repository should expose that regression through tests and either retain a known-good version or address the issue upstream.

## Multi-run study design

The multi-run benchmark freezes target identities from a blueprint simulation before
any condition effect is assigned. Condition runs reuse that molecular blueprint and
change abundance only. This ensures treatment-induced missingness and fold changes are
outcomes of the benchmark rather than criteria used to choose targets.

The canonical study is 25 controls plus 25 treatments with 1,000 frozen targets and
1,000 independent entrapments. Protein treatment effects and replicate variation are
seeded independently from the synthetic-proteome and TimSim blueprint seeds.

Quantification has separate design, realized-input, and post-simulation realized truth.
OpenDIA is evaluated against realized finite-cohort abundance; design-vs-realized
statistics describe biological sampling variation. See `STUDY_BENCHMARK.md`.

## Raw-signal oracle non-circularity

The optional raw-signal oracle is a post hoc diagnostic between simulator truth and
OpenDIA quantification. It integrates only frozen target library fragment m/z values at
simulator-realized RT/IM coordinates from the generated Bruker raw data. It does not use
OpenDIA feature coordinates, scores, q-values, peak boundaries, or intensities during
extraction.

Oracle outputs are not permitted to change the frozen target set, transition library,
entrapment construction, or canonical OpenDIA parameters. This keeps the oracle useful
for locating quantitative response compression without turning raw truth into an
analysis-time oracle.

The interference-sensitivity sweep preserves the same rule. It evaluates a fixed extraction
grid on already generated raw files and ranks transition subsets from collision geometry in
the complete TimSim blueprint fragment universe. Transition ranks are determined from
precursor window group, RT/IM proximity, fragment m/z proximity, and predicted fragment
intensity only. Raw observed intensity and OpenDIA output are excluded from specificity
ranking. Sweep results are diagnostic evidence about cofragmentation and must not be used
to tune canonical OpenDIA parameters.

