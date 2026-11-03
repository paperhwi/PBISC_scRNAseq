import nbformat as nbf
nb = nbf.v4.new_notebook(); C = nb.cells; md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Three-model comparison: scSemiProfiler vs PBISC prev K4000 vs PBISC floor (2026-10-22)

Runs everything that local data allow (the 5 analyses ①–⑤ proposed on 2026-10-21). Code in `code/`; this notebook only reads results.

| # | analysis | script |
|---|---|---|
| ① | the full validation battery on scSemiProfiler cells (integration, C2ST, density/coverage, co-expression, LR, M6, claims, transfer-free rare, 4 integration methods) | `scsemi_prep.py`, `scsemi_pb.py`, `big_joint.py`, `c2st_dc.py`, `coexpr.py`, `lr.py`, `claims_rare.py`, `integ.py`, `m6_s3.py` |
| ② | bulk-broadcast baseline for M6 (sepsis + COVID) and partial r of each model given broadcast | `m6_s3.py` |
| ③ | group-level composition and disease shifts (six-type and transferred fine labels) | `comp_s3.py` |
| ④ | patient-bootstrap CIs: M6, composition, donor identity, per-sample iLISI | `m6_s3.py`, `comp_s3.py`, `boot_s3.py` |
| ⑤ | myeloid gene-level diagnosis (depth-controlled residuals, T/NK co-detection, identity and MHC-II programme) | `myeloid_s3.py` |

**scSemiProfiler cells.** strict_5pairs, draw0 export (5 paired support samples; 1,500 rows per test sample; values = expected counts). The real
held-out panel *is* the scSemi test set, so scSemi cells were restricted to samples whose patient is not in real half A (the reference) and, for sepsis,
to adults: sepsis 34 samples (51,000 cells), COVID 102 samples (60,000 cells). One seeded Poisson draw per cell; a floor-rule copy is kept for a
sensitivity check. **The export is sparse by construction**: ~93 genes with non-zero rate per cell (sepsis; native output zero fraction 0.984;
smallest positive rate 2.04), so scSemi cells carry ~100 detected genes against ~500–900 in real cells.""")
code("""import os, json, numpy as np, pandas as pd
from IPython.display import Image, display
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_colwidth', 90)
R = '/mnt/user-data/outputs/s3/results'; F = '/mnt/user-data/outputs/s3/figs'; TR = '/mnt/user-data/outputs/ilisi_floor_20261021/results/big'
rd = lambda n: pd.read_csv(f'{R}/{n}')""")
md("## Scorecard (sepsis / COVID)")
code("S = rd('scorecard_three_models.csv'); S")
md("## Sparsity of the scSemi export")
code("display(Image(f'{F}/s3_genes_per_cell.png', width=1000))")
md("## ① Integration and cluster matching (pairwise Harmony with real half A)")
code("""J = pd.DataFrame({(c, s): json.load(open(f'{TR}/{c}_{s}_summary.json')) for c in ('sepsis', 'covid') for s in ('realU', 'scsemi', 'prev4k', 'floor')}).T
J[['n_source', 'iLISI_pre', 'iLISI_post', 'own_share_post', 'six_agree_with_own_label', 'src_cells_in_mixed_clusters', 'src_cells_in_src_only_clusters', 'cluster_freq_pearson', 'cluster_freq_JSD']]""")
code("""I = pd.concat([rd(f'{c}_integration_methods.csv') for c in ('sepsis', 'covid')]).drop_duplicates(['cohort', 'pair', 'method'], keep='last')
I.pivot_table(index=['cohort', 'method'], columns='pair', values=['iLISI', 'label_ASW']).round(3)""")
code("display(Image(f'{F}/s3_integration.png', width=1100)); display(Image(f'{F}/s3_umap.png', width=800))")
md("## ① Realism: C2ST, density / coverage")
code("""pd.concat([rd(f'{c}_c2st.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'lineage'], columns='source', values='c2st_auc').round(3)""")
code("""pd.concat([rd(f'{c}_density_coverage.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'space', 'lineage'], columns='source', values=['coverage', 'density']).round(3)""")
code("display(Image(f'{F}/s3_realism.png', width=1100))")
md("## ① Co-expression and ligand–receptor")
code("""pd.concat([rd(f'{c}_coexpression.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'genes', 'lineage'], columns='source', values='corr_of_corr').round(3)""")
code("""display(pd.concat([rd(f'{c}_lr_recovery.csv') for c in ('sepsis', 'covid')]).round(3)); display(pd.concat([rd(f'{c}_lr_disease.csv') for c in ('sepsis', 'covid')]).round(3))""")
md("## ① + ② + ④ Disease effect (M6), bulk broadcast, bootstrap\n`scsemi` = non-A adult test samples; `scsemi_all` = all adult test samples (the same samples the truth is computed on). `broadcast` = bulk RNA-seq logFC of the bulk-only samples the floor model was given (sepsis GSE279448, 121 samples; COVID Imperial ESKD, Gisby 2022 Zenodo 6497251, 164 samples), identical for every cell type.")
code("""M = rd('m6_three_models.csv')
M.groupby(['cohort', 'contrast', 'genes', 'source'])[['logFC_r', 'slope_on_DEG', 'sign_agree_DEG', 'top100_overlap']].median().round(3)""")
code("""PR = rd('m6_broadcast_partial.csv'); PR[PR.contrast == 'case vs control'].groupby(['genes', 'source'])[['r_broadcast_real', 'r_model_real', 'r_model_broadcast', 'partial_r_model_given_broadcast', 'R2_broadcast', 'R2_broadcast_plus_model']].median().round(3)""")
code("""CI = rd('m6_bootstrap_ci.csv'); CI[CI.comparison.str.contains(' - ')].round(3)""")
code("display(Image(f'{F}/s3_m6_bootstrap.png', width=1100)); display(Image(f'{F}/s3_broadcast_partial.png', width=1000))")
md("## ③ + ④ Composition and disease shifts")
code("""CP = rd('composition_shift.csv'); CP.round(3)""")
code("""CC = rd('composition_bootstrap_ci.csv'); CC.round(3)""")
code("display(Image(f'{F}/s3_composition_shift.png', width=1000))")
md("## ① Paper claims and transfer-free rare states")
code("""CL = pd.concat([rd(f'{c}_claims.csv') for c in ('sepsis', 'covid')])
CL.pivot_table(index=['claim', 'auc_realA'], columns='source', values='delta_pos_frac').round(3)""")
code("""RR = pd.concat([rd(f'{c}_rare_transfer_free.csv') for c in ('sepsis', 'covid')])
RR['ratio_vs_realB_same_classifier'] = RR.gen_estimate_pct / RR.realB_estimate_pct
ok = RR[RR.cv_auc_realA >= 0.8].copy(); ok['within_0.5_2'] = ok.ratio_vs_realB_same_classifier.between(0.5, 2)
display(ok.pivot_table(index=['cohort', 'subset'], columns='model_depth', values='ratio_vs_realB_same_classifier').round(2)); ok.groupby(['cohort', 'model_depth'])['within_0.5_2'].agg(['size', 'sum'])""")
md("## ④ Patient identity and per-sample mixing")
code("rd('bootstrap_identity_mixing.csv').round(3)")
md("## ⑤ Myeloid diagnosis")
code("""display(pd.concat([rd(f'{c}_myeloid_model_specific_summary.csv') for c in ('sepsis', 'covid')]).round(3))
display(pd.concat([rd(f'{c}_myeloid_identity.csv') for c in ('sepsis', 'covid')]).round(3))""")
code("""G = rd('sepsis_myeloid_model_specific_genes.csv'); G.groupby('model').head(20).round(3)""")
code("""pd.concat([rd(f'{c}_myeloid_m6_missed_by_category.csv') for c in ('sepsis', 'covid')]).pivot_table(index=['cohort', 'category'], columns='source', values=['n_DEG', 'miss_rate']).round(2)""")
code("display(Image(f'{F}/s3_myeloid.png', width=1100))")
nbf.write(nb, '/mnt/user-data/outputs/s3/three_model_comparison.ipynb')
