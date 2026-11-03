import nbformat as nbf
nb = nbf.v4.new_notebook(); C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Qualifying PBISC cell proportions (analyses B–F), 2026-10-13
Question: is mapping generated cells into the clustering-ready real reference a sound way to qualify PBISC's cell proportions,
and what do the proportions look like when that mapping is checked?

| | Analysis | Files |
|---|---|---|
| B | Deconvolution of the same sepsis bulk (GSE279448) with the same reference (NNLS, weighted NNLS, CIBERSORT-style SVR, WNNLS with a 5-donor bulk→scRNA gene correction, BayesPrism) and a no-information null | `results/sepsis_deconv_props.csv.gz`, `results/B1_paired_error.csv`, `code/b_deconv.py`, `code/b2_*`, `code/b3_bp_collect.py` |
| C | Proportions by four mapping routes with a reject option | `results/*_mapping_counts.csv.gz`, `results/*_mapping.h5ad`, `code/c_mapping.py` |
| D | Confusion matrix from real cells passed through a PBISC-like read-out; recoverable labels | `results/*_confusion_recall.csv`, `results/D_recoverability.csv` |
| E | Label-free neighbourhoods (Milo-style) | `results/*_nhoods*.csv.gz`, `results/*_nhood_summary.json`, `results/*_k1000_nhood_embedding.h5ad`, `code/e_nhood.py` |
| F | Per-sample compositional statistics (CLR), effect sizes with CIs, attenuation slope | `results/compstats_*.{csv.gz,json}`, `results/F_summary.csv`, `code/f_compstats.py` |

COVID bulk (Gisby 2022) is controlled-access (EGA EGAS00001006778) and is not on F:, so B is sepsis only. A (ground truth for PBISC) needs new PBISC runs on held-out pseudobulk.""")
code("""import json, numpy as np, pandas as pd
from IPython.display import Image, display
pd.set_option('display.width', 220); pd.set_option('display.max_rows', 200); pd.set_option('display.max_columns', 30)
R = 'results'""")
md("""## B. Deconvolution baselines (sepsis)
**B1 – paired samples** (bulk with a single-cell partner, adult; 2-fold by sample, signature and calibration donors from the other fold).
Median per-sample L1 and Aitchison distance to the sample's own scRNA composition; `r_across_samples` = median over abundant types of the
across-sample correlation (does the method follow individuals?). `null` = pooled adult composition for every sample.""")
code("""B1 = pd.read_csv(f'{R}/B1_paired_error.csv', keep_default_na=False, na_values=['']); display(B1.pivot_table(index='method', columns='level', values=['L1', 'aitchison', 'r_across_samples']).round(3))""")
md("**B2 – bulk-only samples** (the 122 PBISC was generated from): group composition and disease effect against the real scRNA groups, PBISC vs deconvolution vs null (rows from F). WNNLS+cal uses each LoRA draw's own 5 training donors (five versions, averaged). BayesPrism was run on the paired folds only (uncalibrated and calibrated, one Gibbs pass); the bulk-only BayesPrism run was not done.")
code("""FS = pd.read_csv(f'{R}/F_summary.csv', keep_default_na=False, na_values=[''])
keep = ['null', 'deconv_WNNLS+cal', 'deconv_BayesPrism', 'k1000', 'k4000', 'k1000_R3', 'k4000_R3']
t = FS[(FS.cohort == 'sepsis') & (FS.contrast == 'adult_sepsis_vs_adult_HC') & FS.source.isin(keep)]
display(t.set_index(['level', 'source'])[['L1_ctrl', 'L1_case', 'aitchison_ctrl', 'aitchison_case', 'slope', 'slope_lo', 'slope_hi', 'spearman', 'sign_agree', 'spread_ratio_case']].round(3))
display(Image('figs/prop_B_baselines.png'))""")
md("""## C. Mapping robustness
Routes: **R1** pairwise Harmony + kNN (the report's current method), **R2** kNN in the reference PCA without Harmony,
**R3** logistic regression trained on reference cells passed through the PBISC-like read-out, **R4** marker gating (six lineages).
Reject: distance (R1, R2) or confidence (R3) beyond the value that rejects 5% of PBISC-like real cells (R4 has no usable reject threshold).
`realT2` = held-out real cells with the PBISC-like read-out and known labels (positive control).""")
code("""for coh in ['sepsis', 'covid']:
    L = pd.read_csv(f'{R}/{coh}_mapping_counts.csv.gz')
    print(coh, 'share of cells rejected (lin6 level)')
    display(L[(L.level == 'lin6') & (L.route != 'truth')].groupby(['route', 'source']).apply(lambda d: d.loc[d.rejected, 'n'].sum() / d.n.sum()).unstack().round(3))
    m = json.load(open(f'{R}/{coh}_mapping_meta.json'))['accuracy']
    print(coh, 'cell-level accuracy on held-out real cells'); display(pd.DataFrame(m).T.round(3))""")
code("""def pooled_L1(L, lv, rt):
    tr = L[(L.source == 'realT2') & (L.route == 'truth') & (L.level == lv)].groupby('label').n.sum(); tr = tr / tr.sum()
    p = L[(L.source == 'realT2') & (L.route == rt) & (L.level == lv)].groupby('label').n.sum(); p = p / p.sum(); u = tr.index.union(p.index)
    return float((p.reindex(u, fill_value=0) - tr.reindex(u, fill_value=0)).abs().sum())
rows = []
for coh in ['sepsis', 'covid']:
    L = pd.read_csv(f'{R}/{coh}_mapping_counts.csv.gz')
    for lv in ['lin6', 'major', 'fine']:
        for rt in ['R1', 'R2', 'R3', 'R4']:
            if rt == 'R4' and lv != 'lin6': continue
            rows.append({'cohort': coh, 'level': lv, 'route': rt, 'L1_realT2_pred_vs_truth': pooled_L1(L, lv, rt)})
print('Composition error the mapping alone introduces (real cells, PBISC-like read-out):')
display(pd.DataFrame(rows).pivot_table(index=['cohort', 'level'], columns='route', values='L1_realT2_pred_vs_truth').round(3))
Cr = pd.read_csv(f'{R}/C_routes_major.csv')
display(Cr.pivot_table(index=['cohort', 'label'], columns='route', values='log2_ratio').round(2))
display(Image('figs/prop_C_routes.png'))""")
md("""## D. Which fine labels can be judged at all
Recall and precision of each label on real cells with the PBISC-like read-out. A label is *recoverable* for a route when both are ≥ 0.5;
only for those labels does the generated/real ratio mean something. The NNLS confusion-matrix correction is kept for completeness
but is unstable for fine labels (near-singular confusion matrix) and is not used for conclusions.""")
code("""D = pd.read_csv(f'{R}/D_recoverability.csv')
for (c, r), d in D.groupby(['cohort', 'route']):
    rec = d[d.recoverable]
    print(f'{c} {r}: {len(rec)} of {len(d)} labels (real share >= 0.3%) recoverable; within 1.5x: {(rec.log2_ratio.abs() <= 0.585).sum()}, within 2x: {(rec.log2_ratio.abs() <= 1).sum()}')
display(D.sort_values(['cohort', 'route', 'real_share'], ascending=[True, True, False]).round(3))
display(Image('figs/prop_D_recover.png'))""")
code("""for coh in ['sepsis', 'covid']:
    cp = pd.read_csv(f'{R}/{coh}_corrected_props.csv')
    print(coh, 'R1 major, generated K1000, observed vs NNLS-corrected (case group)')
    g = 'adult_sepsis' if coh == 'sepsis' else 'covid_positive'
    display(cp[(cp.route == 'R1') & (cp.level == 'major') & (cp.source == 'k1000') & (cp.group == g)].set_index('label')[['observed', 'corrected']].round(3).T)""")
md("""## E. Label-free neighbourhoods
Joint Harmony space (real + one K), neighbourhoods of 300 cells around refined index cells; unit = sample.
**E2** (disease): case/control log2 ratio per neighbourhood computed within real samples and within generated samples.
**E1** (abundance, balanced subset): generated cells dominate most neighbourhoods because they are more compact in this space
(within-source 30-NN distance 3.8 vs 5.5, total variance 159 vs 227 for sepsis), so E1 mainly re-states the residual gap.""")
code("""for coh in ['sepsis', 'covid']:
    print(coh); display(pd.DataFrame(json.load(open(f'{R}/{coh}_nhood_summary.json'))).round(3))""")
md("## F. Per-sample compositional statistics\nSlope of generated effect on real effect across abundant cell types (1 = same size, < 1 = attenuated), with bootstrap CI over samples; `spread_ratio_case` = SD across samples of CLR, generated / real.")
code("""cols = ['slope', 'slope_lo', 'slope_hi', 'spearman', 'sign_agree', 'L1_case', 'L1_ctrl', 'spread_ratio_case', 'n_types']
display(FS[FS.source.isin(['k1000', 'k4000', 'k1000_R3', 'k4000_R3'])].set_index(['cohort', 'level', 'contrast', 'source'])[cols].round(3))
S = json.load(open(f'{R}/compstats_summary.json'))
dr = pd.DataFrame([{'key': k.rsplit('|', 1)[0], 'draw': k.rsplit('|', 1)[1], 'slope': v['slope']} for k, v in S.items() if '_draw' in k])
print('slope per draw (each draw = 100 cells per sample)'); display(dr.pivot_table(index='key', columns='draw', values='slope').round(2))
display(Image('figs/prop_EF_effects.png'))""")
md("## Hypotheses: verdicts\nSee the last cell of the analysis summary (`results/verdicts.md`).")
code("print(open(f'{R}/verdicts.md').read())")
nbf.write(nb, "/mnt/user-data/outputs/integration_report_20261012/proportion_qualification/proportion_qualification.ipynb")
