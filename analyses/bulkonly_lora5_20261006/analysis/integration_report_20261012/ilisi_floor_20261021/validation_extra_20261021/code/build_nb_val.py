import nbformat as nbf
nb = nbf.v4.new_notebook(); C = nb.cells; md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Floor read-out: extra validation runnable from local data (2026-10-21)

Every analysis compares the floor read-out (draw0) with the previous K = 4,000 cells and with real controls. Code in `code/`; this notebook only reads results.

| # | analysis | script | data |
|---|---|---|---|
| A | donor fingerprint (COVID, 32 donors with real scRNA) | `pb_analysis.py` | per-sample x six-type pseudobulks of all cells (`data/pb`, built on the PC by `pb_build.py`) |
| B | M6: disease effect per cell type (pseudobulk DE) | `pb_analysis.py` | same |
| C | classifier two-sample test, density / coverage | `c2st_dc.py` | ~60k-cell bundles + joint Harmony embeddings |
| D | within-lineage co-expression | `coexpr.py` | bundles |
| E | ligand–receptor recovery (LIANA consensus) | `lr.py` | bundles |
| F | paper claims, same classifier on every source | `claims_rare.py` | bundles + kNN-transferred parents |
| G | rare states without label transfer | `claims_rare.py` | bundles |
| H | integration methods (Harmony, Scanorama, scVI, BBKNN) | `integ.py` | 20k + 20k cells per pair |

Not possible here: the read-out-only comparison (no depth256 Poisson control on the PC), per-sample truth for sepsis (the 122 bulk-only samples have no scRNA), paired benchmark (server).""")
code("""import os, numpy as np, pandas as pd
from IPython.display import Image, display
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
R = '/mnt/user-data/outputs/val/results'; F = '/mnt/user-data/outputs/val/figs'
rd = lambda n: pd.read_csv(f'{R}/{n}')""")
md("## A. Donor fingerprint (COVID)\nDonor-stable genes picked on real samples only (ICC across donors with ≥ 2 samples); generated and real donor pseudobulks z-scored per gene within source; chance top-1 = 1/32.")
code("rd('covid_donor_fingerprint.csv').drop(columns=['y_cpm_log_real_median_male_female','y_cpm_log_gen_median_male_female'], errors='ignore').round(3)")
md("## B. Disease effect (M6)\nPer six-type, per-sample pseudobulk log2 CPM, Welch t. Real = all scRNA samples (sepsis adults: 93 sepsis / 7 HC / 8 ICU control; COVID 149 positive / 38 negative); generated = bulk-only samples (sepsis 68 / 36 / 17; COVID positive / negative without the 5 adapter-training samples).")
code("""M = rd('m6_disease_effect.csv')
M.groupby(['cohort', 'contrast', 'genes', 'source'])[['logFC_r', 'slope_on_DEG', 'sign_agree_DEG', 'top100_overlap', 'n_real_DEG']].median().round(3)""")
code("""P = rd('m6_precision_top100.csv'); P.pivot_table(index=['cohort', 'contrast', 'six_type'], columns='source', values='precision_top100_gen').round(2)""")
md("`precision_top100_gen` = share of the generated top-100 genes (by |t|) that are real FDR < 0.05 with the same sign. Sepsis vs ICU control has only 32 real DEGs (7–8 control samples), so it is not interpretable.")
code("display(Image(f'{F}/val_m6_donor.png', width=1100))")
md("## C. Classifier two-sample test and density / coverage")
code("""pd.concat([rd(f'{c}_c2st.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'lineage'], columns='source', values='c2st_auc').round(3)""")
code("""pd.concat([rd(f'{c}_density_coverage.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'space', 'lineage'], columns='source', values=['coverage', 'density']).round(3)""")
code("""T = pd.concat([rd(f'{c}_c2st_top_genes.csv') for c in ('sepsis', 'covid')]); T = T[T.source == 'floor']
T.groupby(['cohort', 'lineage', 'direction']).gene.apply(lambda g: ', '.join(g)).unstack()""")
code("display(Image(f'{F}/val_c2st_dc.png', width=1100))")
md("## D–E. Co-expression and ligand–receptor")
code("""pd.concat([rd(f'{c}_coexpression.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'genes', 'lineage'], columns='source', values=['corr_of_corr', 'mean_abs_corr_src']).round(3)""")
code("""display(pd.concat([rd(f'{c}_lr_recovery.csv') for c in ('sepsis', 'covid')]).round(3)); display(pd.concat([rd(f'{c}_lr_disease.csv') for c in ('sepsis', 'covid')]).round(3))""")
code("display(Image(f'{F}/val_coexpr_lr.png', width=1100))")
md("## F–G. Paper claims and rare states without label transfer")
code("""CL = pd.concat([rd(f'{c}_claims.csv') for c in ('sepsis', 'covid')])
CL.pivot_table(index=['claim', 'auc_realA'], columns='source', values='delta_pos_frac').round(3)""")
md("Real half A has a single adult HC sample, so real half B (4 HC samples in sepsis) is the real comparator. SEP-4A (AUC 0.49) and IFN in B cells (no IFN label) are not usable.")
code("""RR = pd.concat([rd(f'{c}_rare_transfer_free.csv') for c in ('sepsis', 'covid')])
RR['ratio_vs_realB_same_classifier'] = RR.gen_estimate_pct / RR.realB_estimate_pct; RR['calibration_realB'] = RR.realB_estimate_pct / RR.realB_truth_pct
RR[['cohort', 'subset', 'parent', 'model_depth', 'cv_auc_realA', 'realB_truth_pct', 'realB_estimate_pct', 'gen_estimate_pct', 'ratio_vs_realB_same_classifier', 'calibration_realB']].round(3)""")
code("""ok = RR[RR.cv_auc_realA >= 0.8].copy(); ok['within_0.5_2'] = ok.ratio_vs_realB_same_classifier.between(0.5, 2)
ok.groupby(['cohort', 'model_depth'])['within_0.5_2'].agg(['size', 'sum'])""")
code("display(Image(f'{F}/val_claims_rare.png', width=1100))")
md("## H. Integration methods")
code("""fn = [f'{R}/{c}_integration_methods.csv' for c in ('sepsis', 'covid') if os.path.exists(f'{R}/{c}_integration_methods.csv')]
if fn:
    I = pd.concat([pd.read_csv(f) for f in fn]); display(I.pivot_table(index=['cohort', 'method'], columns='pair', values=['iLISI', 'label_ASW']).round(3))
    print('median cLISI (six-type) = 1.00 for every method and pair: saturated, not informative, so label ASW is the biology metric.')
if os.path.exists(f'{F}/val_integration.png'): display(Image(f'{F}/val_integration.png', width=1000))""")
code("""S = rd('scorecard.csv') if os.path.exists(f'{R}/scorecard.csv') else None; S""")
nbf.write(nb, '/mnt/user-data/outputs/val/floor_validation_extra.ipynb')
