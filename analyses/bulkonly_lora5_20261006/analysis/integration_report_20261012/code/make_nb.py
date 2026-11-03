import nbformat as nbf
nb = nbf.v4.new_notebook(); C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# PBISC lab-meeting report — batch integration and clustering, K = 1,000 vs K = 4,000 vs public references

Run from the `integration_report_20261012` folder (all paths relative). COVID real data restricted to the published Table S1 roster (187 samples / 61 patients). Results were computed in the cloud workspace by the
scripts in `code/` (A `a_metrics.py`, C `c_joint.py`, D `d_composition.py`, figures `figs.py`); this notebook reloads their
outputs, prints every number shown on the slides and redraws the key panels.

| analysis | input | what |
|---|---|---|
| A. generated integration | `*_harmony_embedding.h5ad` (500 cells/sample/draw), X_pca vs X_pca_harmony | kBET (k=30, 10k queries, permutation null), iLISI(draw), cLISI(L1), Leiden on 10 Harmony PCs (r 0.5/1.0) → ARI/NMI vs predicted L1 and AIFI L2 |
| A'. real references | recap h5ads (sepsis Harmony by sample; COVID per-compartment `X_rep`) | iLISI(sample), cLISI(label) |
| C. joint co-embedding | 100 cells/sample/draw generated + 400 (sepsis) / 300 (COVID) real cells per sample | shared 6k-panel genes, log CP10k on panel totals, PCA fit on real, generated projected, Harmony batch = source; iLISI; kNN (k=30, distance-weighted) label transfer |
| D. composition | transferred labels vs full real reference | per-group composition r / L1, disease-shift (log2 case/control) Spearman and sign agreement, per-draw r, rare subsets |""")
code("""import json, numpy as np, pandas as pd, matplotlib.pyplot as plt
from IPython.display import Image, display
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)
plt.rcParams.update({'font.family': ['Arial', 'Liberation Sans'], 'figure.facecolor': 'white'})
A = json.load(open('results/metrics_A.json'))""")
md("## 1. Generated cells — integration across the five draws (LoRA adapters)")
code("""keys = ['n_cells','n_samples','kBET_pre','kBET_post','kBET_null_post','iLISI_draw_pre','iLISI_draw_post','cLISI_l1_post']
gen = pd.DataFrame({k: A[k] for k in A if '_k' in k}).T
gen[keys].astype(float).round(3)""")
code("display(Image('figs/batch_metrics.png'))\ndisplay(Image('figs/sepsis_gen_umap.png'))\ndisplay(Image('figs/covid_gen_umap.png'))")
md("## 2. Generated cells — Leiden clustering vs labels")
code("""ck = [c for c in gen.columns if c.startswith(('clusters','ARI','NMI','draw_dominated'))]
gen[ck].astype(float).round(3).T""")
code("display(Image('figs/clustering_metrics.png'))\ndisplay(Image('figs/sepsis_cluster_draw.png'))\ndisplay(Image('figs/covid_cluster_draw.png'))")
md("## 3. Real references (recapitulated public data)")
code("""pd.DataFrame({k: A[k] for k in A if 'real' in k}).T.round(3)""")
code("display(Image('figs/real_ref_umap.png'))")
md("## 4. Joint embedding of real and generated cells")
code("""J = {c: json.load(open(f'results/{c}_joint.json'))['mix'] for c in ['sepsis','covid']}
pd.DataFrame(J).round(3)""")
code("display(Image('figs/joint_mixing.png'))\ndisplay(Image('figs/sepsis_joint_umap.png'))\ndisplay(Image('figs/covid_joint_umap.png'))")
md("## 5. Composition after label transfer")
code("""for coh in ['sepsis','covid']:
    Cj = json.load(open(f'results/{coh}_composition.json'))
    rows = {}
    for K in ['k1000','k4000']:
        for col in ['major','label']:
            rows[(K, col)] = {g: Cj[f'{K}_{col}_{g}_r'] for g in Cj['groups']} | {f'L1 {g}': Cj[f'{K}_{col}_{g}_L1'] for g in Cj['groups']} | \\
                {'shift_spearman': Cj[f'{K}_{col}_shift_spearman'], 'shift_sign_agree': Cj[f'{K}_{col}_shift_sign_agree'], 'n_subsets': Cj[f'{K}_{col}_shift_n']}
    print(coh); display(pd.DataFrame(rows).T.round(3))""")
code("""# redraw: major composition per group, real vs K1000 vs K4000
fig, axs = plt.subplots(1, 2, figsize=(14, 4))
for ax, coh in zip(axs, ['sepsis','covid']):
    T = pd.concat({s: pd.read_csv(f'results/{coh}_comp_{s}_major.csv', index_col=0) for s in ['real','k1000','k4000']}, axis=1).fillna(0)
    T.T.plot.bar(stacked=True, ax=ax, legend=False, width=0.85, colormap='tab20'); ax.set_title(coh); ax.set_ylabel('fraction')
plt.tight_layout()""")
code("display(Image('figs/sepsis_composition.png'))\ndisplay(Image('figs/covid_composition.png'))")
md("### Per-draw composition of the case group (draw = LoRA adapter trained on a different 5 donors)")
code("""for coh, case in [('sepsis','adult_sepsis'), ('covid','covid_positive')]:
    Cj = json.load(open(f'results/{coh}_composition.json'))
    print(coh, {K: Cj[f'{K}_major_{case}_r_per_draw'] for K in ['k1000','k4000']})
    t = pd.read_csv(f'results/{coh}_k1000_transfer.csv.gz', keep_default_na=False); t = t[t.group == case]
    display(pd.crosstab(t.major, t.draw, normalize='columns').round(3))""")
code("display(Image('figs/sepsis_draw_comp.png'))\ndisplay(Image('figs/covid_draw_comp.png'))")
md("## 6. Disease-shift concordance (log2 case/control of subset fractions)")
code("""fig, axs = plt.subplots(2, 2, figsize=(10, 9))
for i, coh in enumerate(['sepsis','covid']):
    for j, K in enumerate(['k1000','k4000']):
        d = pd.read_csv(f'results/{coh}_{K}_shift.csv', index_col=0); ax = axs[i, j]
        ax.scatter(d.real_log2FC, d[f'{K}_log2FC'], s=12); ax.axhline(0, c='k', lw=.4); ax.axvline(0, c='k', lw=.4)
        from scipy.stats import spearmanr
        ax.set_title(f'{coh} {K}: rho = {spearmanr(d.real_log2FC, d[f"{K}_log2FC"])[0]:.2f}'); ax.set_xlabel('real'); ax.set_ylabel('generated')
plt.tight_layout()""")
md("## 7. Rare and disease-associated subsets")
code("""for coh in ['sepsis','covid']:
    print(coh); display(pd.read_csv(f'results/{coh}_rare.csv').set_index('subset').round(3))""")
code("display(Image('figs/sepsis_rare.png'))\ndisplay(Image('figs/covid_rare.png'))")
md("""## 8. Calibration checks (added 2026-10-12): real-vs-real baseline and real cells through a PBISC-like read-out

`code/e_checks.py`. Reference = real cells from half of the samples (600 / sample sepsis, 400 / sample COVID). Compared sources, each integrated
with the reference by pairwise Harmony in the same projected PCA space: held-out real samples (ceiling), the same held-out cells depth-matched
(binomial thinning to generated UMI), the same cells with a PBISC-like read-out (kNN-smoothed profile on the cell's own detected genes, matched
depth, one Poisson draw), real cells of the other study (cross-study baseline), and generated K = 1,000 / 4,000 cells. Panel = 5,675 genes present
in the generated panel and in both real datasets (duplicate gene names in the generated panel de-duplicated).""")
code("""C = pd.read_csv('results/check/calibration_summary.csv')
cols = ['cohort','source','umi_median','detected_median','frac_ones','iLISI_pre','iLISI_post','own_share_pre','own_share_post',
        'acc_coarse','acc_major','acc_fine','conf_major_median','conf_fine_median']
C[cols].round(3)""")
code("""for coh in ['sepsis','covid']:
    r = json.load(open(f'results/check/{coh}_checks.json'))
    srcs = [k for k in r if isinstance(r[k], dict) and 'iLISI_post' in r[k]]
    print(coh); display(pd.DataFrame({k: {l: r[k].get(f'iLISI_post_{l}') for l in ['T','NK','B','Myeloid']} for k in srcs}).round(3))
    print({k: v for k, v in r.items() if '3source' in k})""")
code("display(Image('figs/calib_mixing.png'))\ndisplay(Image('figs/calib_umap.png'))\ndisplay(Image('figs/calib_lineage.png'))")
md("""## 9. Paper-claim and rare populations within their parent lineage (added 2026-10-12, after review)

`code/f_rare_parent.py`. Target cells / parent-lineage cells, pooled per group, plus per-sample medians (samples with >= 20 parent cells) and the number
of samples carrying the population. Real = full recapitulated reference (sepsis in_paper_scheme; COVID published 187-sample roster);
generated = kNN-transferred labels in the 100-cells-per-sample-per-draw joint subset. COVID results in this notebook use the 187-sample roster;
the earlier 198-sample outputs are kept in `results/covid198_backup/`.""")
code("""P = pd.read_csv('results/rare_parent_denominators.csv')
cols = ['cohort','subset','claim','parent','group','real_pooled_pct','real_median_sample_pct','real_samples_with_target','real_samples_evaluable',
        'k1000_pooled_pct','k1000_n_target','k1000_n_parent','k1000_samples_with_target','k1000_samples_evaluable',
        'k4000_pooled_pct','k4000_n_target','k4000_n_parent']
P[cols].round(3)""")
code("display(Image('figs/claims_parent.png'))\ndisplay(Image('figs/rare_parent.png'))")
nb["cells"] = C
md("""## 10. Sensitivity to gene selection and PCA space (added 2026-10-13)
Question: is the generated-vs-real gap created by choosing the 2,000 HVGs on the real reference and fitting PCA on real cells only?
Four variants, each followed by the same Harmony (batch = source) and LISI as section 8:
- **V0** reference HVG + reference-fitted PCA (the original setting)
- **V1** joint HVG (seurat_v3, batch_key = source, on both sources) + reference PCA
- **V2** reference HVG + joint PCA (pooled mean/sd, PCA fitted on a balanced 30,000 + 30,000 cell pool, all cells projected)
- **V3** joint HVG + joint PCA

Real controls (held-out, depth-matched, PBISC-like read-out, other study) are processed identically in every variant. No model-based integration.
Code: `code/i_hvg_pca.py`, `code/j_hvgpca_fig.py`. V0 uses a randomized PCA solver, so PBISC-like real (realT2) in COVID differs slightly from section 8 (1.37 vs 1.44).""")
code("""S = pd.read_csv('results/check/hvgpca_sensitivity_summary.csv'); S['variant'] = S.variant.str[:2]
for m in ['iLISI_pre', 'iLISI_post', 'own_share_post', 'iLISI_post_Myeloid', 'hvg_overlap_with_ref_hvg']:
    print(m); display(S.pivot_table(index=['cohort', 'source'], columns='variant', values=m).round(3))""")
code("""# gap to the other real study after Harmony, per variant
for coh, oth in [('sepsis', 'covid'), ('covid', 'sepsis')]:
    t = S[S.cohort == coh].pivot_table(index='source', columns='variant', values='iLISI_post')
    print(coh); display((t.loc[oth] - t.loc[['k1000', 'k4000']]).round(3))""")
code("display(Image('figs/hvgpca_sensitivity.png'))")
nbf.write(nb, "/mnt/user-data/outputs/integration_report_20261012/integration_report.ipynb")
