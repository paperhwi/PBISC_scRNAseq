# Results draft: PBISC-enabled cohort expansion

## PBISC-enabled cohort expansion reveals clinical associations of rare HSPCs

To extend cell-type-resolved clinical association profiling beyond patients with measured single-cell data, we used PBISC to include 29 additional bulk-only patients alongside 37 Real patients. The resulting 66-patient cohort contained 13 deaths. Patients were the independent analytical units, repeated draws were aggregated within patients, and overlapping patients were counted once. All arms were evaluated using the same 43-population testing family within each endpoint. Rare populations were defined by a mean and median Real patient abundance below 1%; this designation did not reduce the multiplicity denominator.

HSPC abundance emerged as a correlate of peak WHO severity in the expanded cohort. This association did not meet nominal significance in the Real-only analysis (Spearman ρ=0.095, P=0.5742), but was pronounced in the additional cohort (ρ=0.659) and met the fixed 43-population FDR criterion in the combined analysis (ρ=0.354, q=0.0168). Source-restricted permutation testing supported the combined association (q=0.01816). HSPCs had already been observed in 35 of 37 Real patients; the newly detected result therefore concerns their clinical association rather than the first detection of the population itself.

The expanded cohort also strengthened statistical support for outcome-associated populations. The association between HSPC abundance and mortality showed similar AUCs in the Real, additional, and combined cohorts (0.810, 0.804, and 0.808, respectively), reaching q=0.01356 in the combined analysis. IgG and IgA antibody-secreting cell abundances showed concordant severity associations across sources and met the combined FDR criterion (ρ=0.391, q=0.007077 and ρ=0.384, q=0.007818, respectively). These results represent FDR rescue of associations already nominally present in Real data, whereas the HSPC–severity association met our stricter Real-missed definition of Real P≥.05 and combined q<.05. Across all 43 populations, severity discoveries changed from 10 to 12, with six gains and four losses; mortality discoveries changed from one to two.

## Cell-type-resolved profiling identifies a candidate HSPC transcriptomic correlate of severity

We next examined all 5,762 shared genes across the nine cell labels meeting coverage criteria, including seven labels satisfying the primary rarity definition. HSPC analyses included 33 Real and 29 additional patients. KIAA0391 expression within the HSPC label fell short of nominal significance in Real-only data (ρ=0.336, P=0.05593), showed a strong association in the additional cohort (ρ=0.803), and met the HSPC-wide gene-family FDR criterion in the combined cohort (ρ=0.571, q=0.007249). Analytic BH adjustment across all nine eligible cell types within the severity endpoint also yielded q=0.04408; adjustment jointly across both endpoints yielded q=0.08816.

A fixed high-precision source-restricted permutation analysis supported the combined KIAA0391 association (9,999,999 resamples; plus-one P=2.0e-06). Applying a Bonferroni factor for the full HSPC family of 5,762 genes across both clinical endpoints retained the association (11,524 tests; adjusted P=0.023048). The Monte Carlo 95% upper bound, multiplied by the same family size, also remained below .05. This sensitivity analysis retained the full HSPC testing denominator rather than restricting correction to selected candidates; it does not constitute permutation-based global BH across all cell types.

Additional transcriptomic FDR rescues included ADK within Tγδ cells for severity (combined ρ=-0.567, within-cell-type q=0.01959; endpoint-global q=0.04408) and MXRA7 within the same label for mortality (AUC=0.041, within-cell-type q=0.03266). Lower MXRA7 values were associated with death; this candidate did not meet the endpoint-global criterion (q=0.266) and involved only three deaths in the additional cohort, warranting supplementary presentation.

Together, these findings support PBISC-enabled inclusion of bulk-only patients for rare-cell and transcriptomic association discovery. They identify candidates for direct follow-up rather than establish patient-level reconstruction fidelity or clinical utility. HSPC severity effects differed in magnitude across sources, and the shared-patient KIAA0391 correlation (ρ=-0.181, n=24) did not support high expression-reconstruction fidelity. Additional-cohort support is part of the discovery analysis, not independent replication. Threshold transitions do not themselves test augmentation benefit, and mortality associations do not establish severity-independent prognosis.

## Optional Discussion paragraph: bulk context

In the same 29 additional patients, the single-gene bulk KIAA0391 score did not show a detectable severity association (ρ=0.167, P=0.3867). This motivates further investigation of cell-type-resolved context. The difference between bulk and inferred-HSPC associations was not formally tested, and the inferred values lack independent measured HSPC ground truth in these patients; the comparison therefore does not establish cell specificity or superiority over bulk profiling.

## Editorial notes

- All q values are retrospective, explicitly defined family-specific results; the analysis contract is not prospective preregistration.
- Rho values refer to the canonical within-source rank-normalized patient scores. Abundance and expression aggregations differ and are detailed in Methods.
- Mortality AUC describes an univariable association, not externally validated prediction.
- The source-specific HSPC severity contrast is descriptive: heterogeneity P=0.007643;43-population-adjusted q=0.06515.
- Shared-patient KIAA0391 bridge correlation: ρ=-0.181, P=0.3961, n=24. Noncontemporaneous sampling limits its interpretation.
- This is a manuscript text proposal, not a report that newly rendered figures were created.
