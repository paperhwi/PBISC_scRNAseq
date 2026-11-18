# Methods and interpretation

## Fixed comparison

Composition: 37 Real patients/7 deaths and 29 genuinely additional patients/6 deaths;66 unique patients,217 draws. Shared patients are represented once using Real measurements. RBC and platelet labels are excluded before closure to100% leukocytes. Draws are averaged equally within each patient; each feature is rank inverse-normal transformed within source before combined association analysis. Severity is patient-level peak WHO severity; mortality is a binary endpoint.

The cohort uses a recoded four-level WHO-derived severity variable: 0=Mild,1=Moderate,2=Severe,3=Critical. These codes are not the original WHO 0–8 scale. Figures with numeric0–3 axes refer to this recoding.

All 43 fine populations remain in every composition BH family, including constant/nonestimable features with P=1. Primary q<.05 is applied separately within each endpoint. Joint86-test and historical death q<.10 sensitivity results are distinct. Rarity is assigned after testing using both Real patient mean and median below 1%; rarity does not reduce the multiplicity denominator.

RNA: 14 labels are shared by the available rare-subset expression matrices. Primary eligibility requires at least 3 cells per patient, at least 10 patients in each source, and at least 3 deaths and 3 survivors per source. Nine labels qualify, of which seven satisfy the primary rarity rule. All 5,762 shared unique symbols are tested within each eligible label. The available panel is not the whole human transcriptome or all 43 composition labels.

Per-cell shared-gene CP10K/log1p values are averaged with cell weighting within patient, then rank inverse-normal transformed within source. Equal-visit aggregation, a minimum10-cell threshold, and duplicate-symbol summation are separate sensitivities. PBISC expression remains inferred expression even when stored as integer-valued counts. Primary duplicate mapping uses the first occurrence; the duplicate-sum sensitivity is retained.

Severity uses two-sided Spearman association. Mortality uses a two-sided Mann–Whitney U test with featurewise exact/asymptotic selection compatible with scalar SciPy `method='auto'`; AUC is U/(deaths×survivors), with deaths as the first group. Larger values in deaths give AUC>.5. Values are never silently flipped into a risk direction. Constant features have rho=0, AUC=.5 and P=1.

## Multiplicity and uncertainty

- Composition: BH43 per endpoint; jointBH86 is separate.
- RNA: BH5,762 within label/endpoint; joint11,524 within label; global51,858 per endpoint; global-joint103,716 across both endpoints. These are separate families and cannot be substituted without naming their scope.
- Modules: Cycling and Cycle4 are exact gene-weight aliases and count once, giving16 unique hypotheses; the legacy17-name result is retained separately. The two-endpoint module family contains 32 tests.
- Historical marker panel: BH20 per endpoint and joint40 are retrospective targeted analyses, not genome-wide discoveries.
- Intervals:50,000 fixed-score percentile patient bootstrap resamples; severity rho is source-stratified and mortality AUC source×outcome-stratified. Maps, generator and labels are held fixed. These intervals omit generator, annotation and candidate-selection uncertainty.

Composition permutation tests use19,999 source-restricted resamples. KIAA0391 high-precision checks use a fixed 9,999,999 resamples per cohort and plus-one P=(k+1)/(B+1). Full-HSPC Bonferroni sensitivities use5,762 or11,524 as the multiplier, including the Monte Carlo95% upper bound. They do not verify permutation-based BH across all nine labels. Raw P, analytic BH q and Bonferroni-adjusted P are labelled separately.

## What the results establish

The findings support cohort extension for clinical association discovery and identify concrete rare-cell and transcriptomic candidates. HSPC abundance–severity and HSPC KIAA0391–severity meet the strict Real-missed definition; IgG/IgA abundance–severity, HSPC abundance–mortality and Tγδ ADK/MXRA7 associations are FDR rescues of nominal Real associations. The New cohort is part of the selection/discovery analysis.

HSPCs already occur in Real data. KIAA0391 Real P=.0559 is close to the nominal threshold. HSPC severity effects are stronger in the added source; source-restricted permutation does not remove all case-mix or timing confounding. Significance in Combined and nonsignificance in Real do not test a difference in effect or prove PBISC-specific augmentation gain.

KIAA0391 shared-patient correlation is negative and nonsignificant (rho≈−.181,n=24); this does not establish a significant inverse relationship, but it does not support high individual expression-reconstruction fidelity. Sampling is not guaranteed contemporaneous. The same 29-patient single-gene bulk comparison motivates further cell-context analysis; a nonsignificant bulk association and significant inferred-cell association do not prove their difference or cell specificity.

Mortality associations are not established as severity-independent. Peak WHO severity is not a baseline prognostic variable. The sparse additional-cohort Tγδ mortality group contains three deaths. No result here establishes an externally validated biomarker, a causal mechanism, a new cell identity, or clinical decision benefit.

## Next decisive experiment

Use patient-held-out masking among patients with paired bulk and observed single-cell data. Keep generator training, reference construction and tuning inside the training fold. Compare reference-only, reference+inferred masked patients, reference+measured masked patients, and bulk/deconvolution controls. Repeat masks and evaluate effect recovery, false positives and candidate recovery with the complete selection procedure inside the evaluation. Same-draw direct HSPC abundance and sorted-HSPC KIAA0391 measurement would test the cell-level interpretation.
