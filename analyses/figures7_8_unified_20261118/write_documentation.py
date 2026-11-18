"""Generate publication-facing documentation from the frozen numerical tables."""
from pathlib import Path
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'

def main():
    primary=pd.read_csv(DATA/'PRIMARY_CANDIDATES.csv')
    rna=pd.read_csv(DATA/'transcriptome/fixed_primary_candidates_with_rarity_and_composition.csv')
    comp=pd.read_csv(DATA/'composition/candidate_evidence_all_rescues.csv')
    def pick(t,**kw):
        for key,value in kw.items():t=t[t[key].eq(value)]
        assert len(t)==1
        return t.iloc[0]
    hs=pick(primary,level='composition',feature='HSPC',endpoint='severity')
    hd=pick(primary,level='composition',feature='HSPC',endpoint='death')
    k=pick(rna,celltype='HSPC',gene='KIAA0391',endpoint='severity')
    kp=pick(comp,feature='HSPC',endpoint='severity',arm='Combined')
    fmt=lambda x:f'{x:.4g}'
    readme=f'''# Figure 7–8: unified rare-cell clinical association analysis

**PBISC-enabled inclusion of bulk-only patients reveals rare-HSPC and candidate transcriptomic associations with COVID-19 outcomes.**

This release implements the approved five-panel Figure 7, five-panel Figure 8, and five supplementary figures using the canonical `25-unified-v1.0` results. It is the current unified-analysis figure package; the repository's earlier Figure 7 versions remain historical analyses with different rules and scope.

## Read and view

- [Results and complete legends](docs/Results_and_Legends.md)
- [All seven figures in one PDF](figures/Figures_7_8_and_S1_S5.pdf)
- [한국어 서술·패널 구성안](docs/Positive_Narrative_KO.md), [English Results](docs/Results_Draft_EN.md)
- [Analysis methods and claim boundaries](docs/Methods_and_Interpretation.md)
- [Independent release checks](provenance/RELEASE_VALIDATION.json)

### Figure 7: rare-cell clinical associations

![Figure 7](figures/Figure_7.png)

[PDF](figures/Figure_7.pdf) · [editable SVG](figures/Figure_7.svg) · [panel source data](source_data/)

### Figure 8: a candidate HSPC transcriptomic correlate of severity

![Figure 8](figures/Figure_8.png)

[PDF](figures/Figure_8.pdf) · [editable SVG](figures/Figure_8.svg) · [panel source data](source_data/)

### Supplementary figures

| Figure | Content | Files |
|---|---|---|
| S1 | All43 composition results, rarity and gain/loss context | [PNG](figures/Supplementary_Figure_S1.png) · [PDF](figures/Supplementary_Figure_S1.pdf) · [SVG](figures/Supplementary_Figure_S1.svg) |
| S2 | All14 RNA-label coverage and nine-label shared-gene screen | [PNG](figures/Supplementary_Figure_S2.png) · [PDF](figures/Supplementary_Figure_S2.pdf) · [SVG](figures/Supplementary_Figure_S2.svg) |
| S3 | Tγδ ADK–severity and MXRA7–mortality FDR rescues | [PNG](figures/Supplementary_Figure_S3.png) · [PDF](figures/Supplementary_Figure_S3.pdf) · [SVG](figures/Supplementary_Figure_S3.svg) |
| S4 | Same29-patient bulk KIAA0391 context and shared-patient comparison | [PNG](figures/Supplementary_Figure_S4.png) · [PDF](figures/Supplementary_Figure_S4.pdf) · [SVG](figures/Supplementary_Figure_S4.svg) |
| S5 | All16 unique HSPC modules and historical20-gene panel | [PNG](figures/Supplementary_Figure_S5.png) · [PDF](figures/Supplementary_Figure_S5.pdf) · [SVG](figures/Supplementary_Figure_S5.svg) |

S1–S5 identify supplements within this release; they do not overwrite the supplements of earlier figure versions.

## Main evidence

| Finding | Real-only | Combined | Interpretation |
|---|---|---|---|
| Rare HSPC abundance–severity | ρ={hs.Real_effect:.3f}, P={fmt(hs.Real_P)} | ρ={hs.Combined_effect:.3f}, q43={fmt(hs.Combined_q)} | Strict Real-missed association; source-permutation q43={fmt(kp.permutation_q43)} |
| HSPC KIAA0391–severity | ρ={k.Real_effect:.3f}, P={fmt(k.Real_P)} | ρ={k.Combined_effect:.3f}, within-HSPC q={fmt(k.Combined_q_gene)} | Shared-panel transcriptomic candidate |
| HSPC abundance–mortality | AUC={hd.Real_effect:.3f}, q43={fmt(hd.Real_q)} | AUC={hd.Combined_effect:.3f}, q43={fmt(hd.Combined_q)} | FDR rescue with concordant effects |

Composition includes37 Real+29 additional=66 patients/13 deaths/217 draws. HSPC RNA includes33 Real+29 additional=62 patients. Generated cells and repeated draws are not independent sample-size units. HSPCs were already observed in35/37 Real patients: the new result concerns the clinical association. Rarity means both the Real mean and median leukocyte abundance are below1%.

KIAA0391 passes analytic endpoint-global BH across9×5,762 tests (q={fmt(k.Combined_q_global)}), but not the two-endpoint global-joint family (q={fmt(k.Combined_q_global_joint)}). Fixed10M-scale source-restricted permutation gives Combined P=2.0×10⁻⁶ and full-HSPC two-endpoint Bonferroni P=.023048. This is not a full-screen permutation-BH calculation. The shared24-patient expression correlation is ρ={k.bridge_rho:.3f}; high individual expression-reconstruction fidelity is not established.

## Reproduce the figures

```bash
python -m pip install -r requirements.txt
python verify_release.py
python render_all.py
python verify_release.py --figures
```

All paths used for rendering are relative to this package. Original H5AD, NPZ patient matrices, generator weights, private runtime directories, and the original patient-ID mapping are not needed for figure regeneration and are not bundled. Plotting vectors use stable release-local pseudonyms; age and sex are excluded. The full primary association table contains155,574 cohort×cell-type×gene records, including null results.

`verify_release.py` independently recomputes all43 composition associations from the bundled patient scores, verifies RNA candidate statistics and the shared bridge, checks120 RNA BH column/family comparisons, and checks module and permutation correction formulas. It does not rerun all raw-expression calculations or the generator. Historical raw-data audits and the fixed10M calculation are documented in `provenance/`.

The optional `stage_source_data.py --canonical PATH` exports the frozen source release again. It needs the original canonical directory; it is not part of ordinary figure regeneration. `write_documentation.py` regenerates publication documentation from the bundled tables.

## Interpretation

Strict Real-missed: Real P≥.05 and Combined primary q<.05. FDR rescue: Real q≥.05 and Combined primary q<.05; these definitions overlap. These labels describe observed threshold transitions, not a formal test of PBISC-specific augmentation benefit. The additional cohort contributes to discovery and is not independent replication. HSPC severity effects differ in magnitude across sources. Mortality AUCs are univariable associations and do not establish severity-independent prognosis or clinical utility.
'''
    (ROOT/'README.md').write_text(readme,encoding='utf-8')

    methods='''# Methods and interpretation

## Fixed comparison

Composition:37 Real patients/7 deaths and29 genuinely additional patients/6 deaths;66 unique patients,217 draws. Shared patients are represented once using Real measurements. RBC and platelet labels are excluded before closure to100% leukocytes. Draws are averaged equally within each patient; each feature is rank inverse-normal transformed within source before combined association analysis. Severity is patient-level peak WHO severity; mortality is a binary endpoint.

The cohort uses a recoded four-level WHO-derived severity variable:0=Mild,1=Moderate,2=Severe,3=Critical. These codes are not the original WHO0–8 scale. Figures with numeric0–3 axes refer to this recoding.

All43 fine populations remain in every composition BH family, including constant/nonestimable features with P=1. Primary q<.05 is applied separately within each endpoint. Joint86-test and historical death q<.10 sensitivity results are distinct. Rarity is assigned after testing using both Real patient mean and median below1%; rarity does not reduce the multiplicity denominator.

RNA:14 labels are shared by the available rare-subset expression matrices. Primary eligibility requires at least3 cells per patient, at least10 patients in each source, and at least3 deaths and3 survivors per source. Nine labels qualify, of which seven satisfy the primary rarity rule. All5,762 shared unique symbols are tested within each eligible label. The available panel is not the whole human transcriptome or all43 composition labels.

Per-cell shared-gene CP10K/log1p values are averaged with cell weighting within patient, then rank inverse-normal transformed within source. Equal-visit aggregation, a minimum10-cell threshold, and duplicate-symbol summation are separate sensitivities. PBISC expression remains inferred expression even when stored as integer-valued counts. Primary duplicate mapping uses the first occurrence; the duplicate-sum sensitivity is retained.

Severity uses two-sided Spearman association. Mortality uses a two-sided Mann–Whitney U test with featurewise exact/asymptotic selection compatible with scalar SciPy `method='auto'`; AUC is U/(deaths×survivors), with deaths as the first group. Larger values in deaths give AUC>.5. Values are never silently flipped into a risk direction. Constant features have rho=0, AUC=.5 and P=1.

## Multiplicity and uncertainty

- Composition: BH43 per endpoint; jointBH86 is separate.
- RNA: BH5,762 within label/endpoint; joint11,524 within label; global51,858 per endpoint; global-joint103,716 across both endpoints. These are separate families and cannot be substituted without naming their scope.
- Modules: Cycling and Cycle4 are exact gene-weight aliases and count once, giving16 unique hypotheses; the legacy17-name result is retained separately. The two-endpoint module family contains32 tests.
- Historical marker panel: BH20 per endpoint and joint40 are retrospective targeted analyses, not genome-wide discoveries.
- Intervals:50,000 fixed-score percentile patient bootstrap resamples; severity rho is source-stratified and mortality AUC source×outcome-stratified. Maps, generator and labels are held fixed. These intervals omit generator, annotation and candidate-selection uncertainty.

Composition permutation tests use19,999 source-restricted resamples. KIAA0391 high-precision checks use a fixed9,999,999 resamples per cohort and plus-one P=(k+1)/(B+1). Full-HSPC Bonferroni sensitivities use5,762 or11,524 as the multiplier, including the Monte Carlo95% upper bound. They do not verify permutation-based BH across all nine labels. Raw P, analytic BH q and Bonferroni-adjusted P are labelled separately.

## What the results establish

The findings support cohort extension for clinical association discovery and identify concrete rare-cell and transcriptomic candidates. HSPC abundance–severity and HSPC KIAA0391–severity meet the strict Real-missed definition; IgG/IgA abundance–severity, HSPC abundance–mortality and Tγδ ADK/MXRA7 associations are FDR rescues of nominal Real associations. The New cohort is part of the selection/discovery analysis.

HSPCs already occur in Real data. KIAA0391 Real P=.0559 is close to the nominal threshold. HSPC severity effects are stronger in the added source; source-restricted permutation does not remove all case-mix or timing confounding. Significance in Combined and nonsignificance in Real do not test a difference in effect or prove PBISC-specific augmentation gain.

KIAA0391 shared-patient correlation is negative and nonsignificant (rho≈−.181,n=24); this does not establish a significant inverse relationship, but it does not support high individual expression-reconstruction fidelity. Sampling is not guaranteed contemporaneous. The same29-patient single-gene bulk comparison motivates further cell-context analysis; a nonsignificant bulk association and significant inferred-cell association do not prove their difference or cell specificity.

Mortality associations are not established as severity-independent. Peak WHO severity is not a baseline prognostic variable. The sparse additional-cohort Tγδ mortality group contains three deaths. No result here establishes an externally validated biomarker, a causal mechanism, a new cell identity, or clinical decision benefit.

## Next decisive experiment

Use patient-held-out masking among patients with paired bulk and observed single-cell data. Keep generator training, reference construction and tuning inside the training fold. Compare reference-only, reference+inferred masked patients, reference+measured masked patients, and bulk/deconvolution controls. Repeat masks and evaluate effect recovery, false positives and candidate recovery with the complete selection procedure inside the evaluation. Same-draw direct HSPC abundance and sorted-HSPC KIAA0391 measurement would test the cell-level interpretation.
'''
    (ROOT/'docs/Methods_and_Interpretation.md').write_text(methods,encoding='utf-8')

    results=(ROOT/'docs/Results_Draft_EN.md').read_text(encoding='utf-8')
    results=results.split('## Editorial notes')[0]
    captions='''
## Figure legends

### Figure7. PBISC-enabled cohort expansion reveals rare-cell clinical associations

**a**, Cohort design:37 Real patients and29 additional bulk-only patients,66 unique patients and13 deaths. Shared patients are counted once; repeated draws and generated cells are not independent observations. **b**, Patient-level HSPC abundance across peak WHO severity, with source distinguished by marker and color. Abundance is leukocyte percentage after RBC/platelet exclusion and equal-draw patient averaging; the association tests use within-source rank-normalized values. **c**, HSPC–severity Spearman rho in Real, New and Combined, with the corresponding P, BH43 q and source-restricted permutation support. **d**, HSPC–mortality AUC in each cohort. **e**, IgG and IgA ASC–severity associations. Error bars are95% fixed-score percentile intervals from50,000 patient bootstrap resamples, stratified by source for rho and source×outcome for AUC. Generator, annotation and feature-selection uncertainty are not included. AUC>.5 indicates larger values in deceased patients and is an univariable association, not validated prediction. All43 populations are retained in endpoint-specific BH families. Strict Real-missed means Real P≥.05 and Combined q<.05; FDR rescue means Real q≥.05 and Combined q<.05. These threshold transitions are not a formal test of augmentation benefit. The complete screen, including lost and null associations, is shown in Supplementary FigureS1 and the bundled tables.

### Figure8. Shared-gene profiling identifies an HSPC transcriptomic correlate of severity

**a**, RNA coverage and testing scope:14 shared labels, nine primary-eligible labels including seven rare labels, and5,762 shared genes per eligible label. HSPC analyses include33 Real and29 additional patients. **b**, KIAA0391–severity rho across Real, New and Combined, with95% fixed-score source-stratified50,000-resample bootstrap intervals. **c**, Patient-level KIAA0391 within the HSPC label and peak WHO severity in the additional29 patients. **d**, Distinct multiplicity scopes for the analytic screen and high-precision permutation sensitivity. Analytic within-HSPC BH, nine-label endpoint-global BH, global-joint BH, and permutation-based full-HSPC two-endpoint Bonferroni are reported separately. Fixed9,999,999 source-restricted resamples per cohort use a plus-one P value. Monte Carlo uncertainty concerns simulation precision, not the population-effect confidence interval. No full-screen permutation-global BH was performed. **e**, Real and PBISC HSPC KIAA0391 patient means in24 shared patients with at least three cells per source. Sampling is not guaranteed contemporaneous. This comparison does not establish high individual reconstruction fidelity. The additional cohort participates in discovery and is not independent replication.

### Supplementary FigureS1. Full composition screen and rarity context

All43 populations are retained, including gains, losses and null associations. Effects are Spearman rho for peak WHO severity and AUC for mortality, with Real/New/Combined in fixed order. Rarity is determined by both Real mean and median below1%, independently of the testing denominator. Significance labels refer to endpoint-specific BH43 q<.05. Full numerical P, q, raw proportions,50,000-resample fixed-score intervals and source-permutation results are bundled under `data/composition/`. Severity discoveries change10→12, with six gains and four losses; mortality discoveries change1→2. These counts do not measure discovery accuracy.

### Supplementary FigureS2. Coverage and complete shared-panel RNA screen

Coverage is reported for all14 common labels at the primary three-cell-per-patient threshold, including labels ineligible for outcome testing. Primary eligibility requires at least10 patients, three deaths and three survivors per source. Counts of combined discoveries are shown for all nine eligible labels with within-label5,762-gene and endpoint-global51,858-test corrections distinguished. Seven labels satisfy the primary rarity rule; CD14mono_IFN and NK CD56bright do not. Ineligible coverage is not evidence of a biological null. The complete primary table includes155,574 cohort×cell-type×gene records.

### Supplementary FigureS3. Additional rare-cell transcriptomic FDR rescues

ADK–severity and MXRA7–mortality effects in Tγδ cells are shown by cohort with95% fixed-score patient bootstrap intervals. Tγδ primary RNA includes34 Real and24 additional patients; the latter include three deaths. ADK has a negative severity association. MXRA7 AUC is retained in its original direction: values below.5 indicate lower expression in deceased patients. Both genes had Real nominal P<.05. ADK passes endpoint-global analytic BH; MXRA7 passes only the within-cell-type gene family. Neither result is a new Real-missed nominal association or independent replication.

### Supplementary FigureS4. Same-patient bulk context and shared-patient expression comparison

The same29 additional patients are used for bulk KIAA0391 and PBISC-HSPC comparisons. Raw HTSeq gene counts exclude technical rows; version-stripped Ensembl mapping identifies ENSG00000100890. The primary bulk score is the equal-draw patient mean of log2(CPM+1), using all60,605 available raw Ensembl rows in the library denominator after the documented mapping/aggregation. The source tables also retain shared-panel normalization and HSPC-cell-count-weighted sensitivities. Shared24-patient Real/PBISC HSPC values describe expression fidelity without asserting identical sampling times. Differences in significance between bulk and inferred HSPC do not test a difference in association or establish cell specificity.

### Supplementary FigureS5. Historical HSPC modules and targeted marker panel

All16 unique historical HSPC module hypotheses are shown for both endpoints; Cycling and Cycle4 are one exact gene-weight alias. Primary and joint-endpoint families are kept distinct. The historical20-gene panel is presented as a retrospective targeted analysis with its own BH20 denominator; the same genes' full-panel gene-family q values retain the5,762 denominator. Targeted-panel positives are not full-panel discoveries. Ery5–severity passes primary q16 but not joint32 or equal-visit sensitivity; Ery5–mortality does not meet primary q16<.05. The complete module and marker tables retain all configurations and null results.
'''
    authored_supplement=ROOT/'supplementary_specs/captions.md'
    if authored_supplement.exists():
        supplement=authored_supplement.read_text(encoding='utf-8')
        supplement=re.sub(r'^## Supplementary_Figure_S([1-5])',r'### Supplementary FigureS\1',supplement,flags=re.M)
        captions=captions.split('### Supplementary FigureS1.')[0]+supplement
    # Readable prose spacing, preserving gene symbols and scientific notation.
    captions=re.sub(r'(?<=[a-z])(?=\d)', ' ', captions)
    captions=captions.replace('below.5','below .5')
    (ROOT/'docs/Results_and_Legends.md').write_text(results+captions,encoding='utf-8')

    # Archival context keeps its text; unavailable historical links become plain
    # canonical-path references instead of broken links in this portable release.
    for name in ['Canonical_Results_KO.md','KIAA0391_10M_Interpretation_KO.md','Bulk_KIAA0391_Context.md']:
        path=ROOT/'docs'/name
        text=path.read_text(encoding='utf-8')
        def replace_link(match):
            label,target=match.groups()
            if target.startswith(('http:','https:','#')):return match.group(0)
            candidates=[ROOT/'data'/target, ROOT/'provenance'/target, ROOT/target]
            found=next((c for c in candidates if c.exists()),None)
            if found:
                return f'[{label}](../{found.resolve().relative_to(ROOT).as_posix()})'
            return f'{label} (canonical source: `{target}`)'
        path.write_text(re.sub(r'\[([^\]]+)\]\(([^)]+)\)',replace_link,text),encoding='utf-8')
    # Typeset prose quantities without changing gene symbols or q-family names.
    prose_words='all|All|and|least|below|across|contains|includes|of|than|from|for|within|with|the|same|uses|Real|New|Combined|WHO|both|unique|shared|original|fixed|nine'
    for path in [ROOT/'README.md',ROOT/'docs/Methods_and_Interpretation.md',ROOT/'docs/Results_and_Legends.md']:
        text=path.read_text(encoding='utf-8')
        text=re.sub(r'\b('+prose_words+r')(?=\d)',r'\1 ',text)
        text=text.replace('FigureS','Figure S').replace('Figure7.','Figure 7.').replace('Figure8.','Figure 8.')
        text=text.replace('Composition:37','Composition: 37').replace('RNA:14','RNA: 14').replace('variable:0=','variable: 0=')
        path.write_text(text,encoding='utf-8')
    print('README, Results/legends and Methods written from frozen release tables')

if __name__=='__main__':main()
