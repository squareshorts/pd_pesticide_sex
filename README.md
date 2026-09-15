# Sex-stratified occupational pesticide exposure and Parkinson's disease MR

This repository contains the analysis code and derived results for a sex-stratified two-sample Mendelian randomization study of occupational pesticide exposure and Parkinson's disease (PD).

## Study design

Exposure associations were obtained from sex-specific Neale Lab UK Biobank round-2 GWAS for occupational phenotypes. The primary exposure was UK Biobank Data-Field 22614, category 1: "Worked with pesticides: Sometimes". Secondary phenotypes were pesticides rarely/never (22614_0), pesticides often (22614_2), and workplace chemical or other fumes often (22610_2).

Sex-specific PD outcome summary statistics were from the International Parkinson's Disease Genomics Consortium dataset described by Blauwendraat et al. (2021), using versions with all UK Biobank participants removed.

Instruments were selected at p < 1e-5, excluding low-confidence variants, requiring MAF >= 0.01 and F > 10, followed by PLINK clumping at r2 < 0.001 within 10 Mb using a European 1000 Genomes LD reference.

The primary estimator was multiplicative random-effects inverse-variance weighted MR. Sensitivity analyses included weighted median, MR-Egger, MR-RAPS, MR-PRESSO, Cochran's Q, radial residual screening, leave-one-out analysis, and weak-instrument-robust inference.

## Primary result

For the primary 22614_1 phenotype, 21 female and 16 male instruments remained after harmonization.

- Female IVW: beta = -1.9652, SE = 1.4671, p = 0.1804.
- Male IVW: beta = 0.1352, SE = 0.9539, p = 0.8873.
- Female-minus-male contrast: delta beta = -2.1004, SE = 1.7499, p = 0.2300, 95% CI [-5.5302, 1.3293].

The primary analysis therefore provided no evidence that sex modifies the MR association. The secondary occupational-exposure phenotypes likewise showed no significant sex interaction.

These null interaction results should not be interpreted as evidence of equivalence between females and males, because the interaction confidence intervals remain broad. The exposure GWAS are linear-probability models of broad occupational questionnaire categories; MR coefficients should not be interpreted as effects per unit of chemical dose.

## Repository structure

- `scripts/`: acquisition, preparation, validation, analysis, and summary scripts.
- `results/primary_results.csv`: primary sex-specific MR estimates.
- `results/sex_interactions.csv`: formal female-minus-male interaction tests.
- `results/sensitivity_summary.csv`: sensitivity-analysis summary.
- `results/*_analysis.json`: complete phenotype-by-sex analysis outputs.
- `results/*_harmonized.csv`: harmonized instrument datasets.
- `results/*_loo.csv`: leave-one-out results.
- `results/*_radial.csv`: radial-MR diagnostic outputs.

## Software

Analyses were run in R 4.5.2. MR-PRESSO 1.0 and mr.raps 0.4.3 were used for sensitivity analyses. LD clumping and retained-instrument LD checks used PLINK v1.9.0-b.7.15.

## Data availability

The repository does not redistribute the source GWAS files or the 1000 Genomes LD reference. Exposure GWAS are available from the Neale Lab UK Biobank resource. PD summary statistics are available from the relevant IPDGC data release subject to its access conditions.

## Manuscript

The associated manuscript evaluates whether sex modifies the association between genetically proxied occupational pesticide exposure and PD. The principal conclusion is that sex-matched Mendelian randomization provides no evidence of such effect modification for the broad UK Biobank occupational exposure phenotypes analyzed here.
