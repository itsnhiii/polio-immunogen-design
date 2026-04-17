# Polio Immunogen Design (BoltzGen)

Exploring paratope-conditioned generative design of protein immunogens for poliovirus.

## Overview

This project uses BoltzGen to design de novo protein scaffolds that mimic antibody-recognized epitopes on poliovirus. The focus is on:

* Avoiding Ig-domain collapse
* Minimal anchor-based binding constraints
* Geometry-aware scaffold topology (e.g., B → A → C)
* Designing compact, manufacturable immunogens

## Structure

```
configs/        # YAML configs for BoltzGen experiments
notes/          # experiment logs and observations
scripts/        # helper scripts (optional)
```

## Key Ideas

* Use **minimal anchor residues** (not full epitope)
* Control structural bias via **structure_groups**
* Break antibody priors to avoid **Ig-like folds**
* Use **loop placement** to connect spatial binding clusters

## Running Experiments

Configs in `configs/` are run with BoltzGen. Each config represents a specific hypothesis:

* anchor-only vs patch visibility
* topology (B→A→C vs alternatives)
* length constraints
* structure group settings

## Experiment Tracking

Each run should be logged in `notes/` with:

* config used
* key settings (anchors, visibility, length)
* outcome (Ig collapse vs non-Ig)
* binding behavior (qualitative)

Example:

```
Run 001:
- config: b_a_c_topology.yaml
- structure_groups: anchors only
- result: no Ig collapse
- binding: weak
- notes: cluster A not fully engaged
```

## Current Focus

* Breaking Ig-domain collapse
* Achieving stable canyon-binding topology
* Improving binding without reintroducing structural bias

## Future Directions

* Multi-objective design (binding + developability)
* Guided diffusion with surrogate models
* Experimental validation via yeast surface display
