# Sex-specific pesticide–Parkinson's disease MR: primary analysis

Exposure: UK Biobank field 22614, category “Sometimes”, sex-specific Neale Lab v3/v2 summary statistics.
Outcome: Blauwendraat et al. sex-stratified Parkinson's disease GWAS, NO_UKB_AT_ALL.
Instrument selection: p < 1e-5; low-confidence variants excluded; MAF >= 0.01; F > 10; 1000 Genomes European LD clumping r2 < 0.001 within 10,000 kb.
Harmonization: allele-matched; palindromic variants with ambiguous MAF removed; variants absent from the outcome removed.

Female: 16 harmonized instruments. IVW-MRE beta=-1.2758, SE=1.8729, p=0.4958, OR=0.279 (95% CI 0.007–10.970); Q p=0.204. MR-Egger intercept p=0.463. Weighted-median p=0.650.

Male: 15 harmonized instruments. IVW-MRE beta=-1.6708, SE=0.8691, p=0.05454, OR=0.188 (95% CI 0.034–1.033); Q p=0.683. MR-Egger intercept p=0.631. Weighted-median p=0.077.

Formal sex heterogeneity: beta_female - beta_male = 0.3950, SE=2.0648, z=0.191, p=0.8483; 95% CI [-3.6519, 4.4419].

Interpretation: no evidence that sex modifies the MR estimate. This is not evidence of equivalence because the interaction confidence interval is wide. The sex-specific analyses do not reproduce the TCC's nominal female-positive IVW signal. The male estimate trends in the opposite (protective) direction but misses conventional significance in standard IVW, Egger, weighted-median, and weak-instrument-aware profile-likelihood analyses. Therefore, the data do not support a female-specific pesticide susceptibility claim.

Important scale caveat: the exposure GWAS is a binary occupational questionnaire category, so the MR coefficient is per unit change in genetically predicted liability/probability for reporting “Sometimes” exposure; the large OR scale should not be interpreted as a toxicological dose-response effect.
