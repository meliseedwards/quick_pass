# GP2 NULISA proteomics pipeline

Quality control and differential abundance analysis for the GP2 NULISA pilot,
run on the Alamar NULISA platform (Neuro220 panel, ~220 protein targets) in
plasma.

The pipeline compares Parkinson's disease cases to controls within each cohort
and produces per-target effect sizes, standard errors and p-values alongside
per-target quality annotations, so that summary statistics from different
cohorts can be combined in a meta-analysis.

## Approach

Differential abundance runs first; quality control annotates the results
afterwards. Every target that can be modelled gets a p-value, and QC flags are
carried alongside it. This is deliberate: a target that is poorly measured is
flagged rather than silently dropped, so it is visible whether a finding rests
on a well-measured protein or a marginal one.

Only data that cannot enter a model is removed. For example, samples with no
resolvable donor ID, samples without a case/control label or complete covariates, 
and targets with no measurements at all. Everything else is retained and annotated.

## Scripts

Run in order. All cohort-specific settings live in `00_config_v2.R`; changing the
cohort there is the only edit needed to run the chain on a different dataset.

| Script | Purpose |
|---|---|
| `00_config_v2.R` | Cohort settings, file paths, thresholds, target categories |
| `01_load_filter_v2.R` | Load NPQ, resolve donor IDs, compute detectability, build a SummarizedExperiment |
| `02a_merge_v2.R` | Join GP2 release metadata, harmonise phenotype and sex |
| `02b_SC_drift_v2.R` | Intra- and inter-plate CV per target from sample control wells |
| `03_DA_v2.R` | Differential abundance, sensitivity refits, QC annotation |
| `04_QC_v2.R` | Plate ICC, covariate screen, zero report, QC figures |

## Input data

- **NPQ**:  NULISA Protein Quantification values as delivered by Alamar,
  already internal-control and inter-plate-control normalised and log2-scaled.
  No additional normalisation is applied.
- **GP2 release metadata**: phenotype, age, sex and ancestry, joined by donor
  ID. Where a donor is not yet in the release, these fields are taken from the
  cohort manifest and the source is recorded per sample.
- **Cohort manifests**: map site-local sample identifiers to GP2 donor IDs.

Input data are controlled-access and are not included in this repository.

## Model

```
NPQ ~ phenotype + age + sex
```

Fitted per target with `limma::lmFit` without empirical Bayes moderation, so
that ordinary standard errors and degrees of freedom are available for
meta-analysis. Plate can be included as a covariate, but only where phenotype
is reasonably balanced across plates; where it is not, adjusting for plate
removes much of the case/control contrast rather than a batch effect. Each
cohort's plate balance is tested and reported.

Two sensitivity refits accompany every result: one excluding samples with an
assay QC warning, one excluding PCA outliers. A hit is marked robust only if it
survives both.

## Quality annotations

Carried on every differential abundance result:

- **Detectability**: proportion of samples above the target's limit of detection
- **Intra- and inter-plate CV** — measured on sample control wells, computed on
  unlogged values
- **Plate ICC** with an F-test and Tukey pairwise comparison; a target is called
  batch-affected only when all three criteria are met
- **Zero imbalance**: the difference in exact-zero rate between cases and
  controls, which distinguishes a genuine effect from a detection-floor artifact
- **Vendor flags**: cross-reactivity, low vendor detectability and high vendor
  CV, from the Alamar Neuro220 data sheet

## Environment

R 4.5.2. See `environment.yml` for the full package set. Key dependencies are
SummarizedExperiment, limma, lme4 and tidyverse.

## Status

This pipeline is under active development as additional cohorts are processed.
Results produced with it are preliminary.
