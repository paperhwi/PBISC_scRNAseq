import nbformat as nbf
nb = nbf.v4.new_notebook(); C = nb.cells; md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Floor read-out validation for the deck (2026-10-21)

Checks behind `PBISC floor read-out 검증 보고서.pptx`. Larger bundles (~60,000 cells per source, `data/big/`, `code/big_prep.py`):
real held-out cells split by patient into half A (reference) and half B (control) exactly as `ilisi_prep.py`; floor read-out (draw0) and previous K = 4,000 cells drawn equally per sample.

* `code/big_txstats.py` — gene count curve, per-cell shape, pseudobulk r per six-type, marker panel (adds real half B thinned to floor / prev depth).
* `code/big_joint.py` — pairwise Harmony with half A (same V0 pipeline as `ilisi_run.py`), iLISI re-run, kNN (k = 30) label transfer, joint Leiden, UMAP. Joint embeddings: `data/big_joint/<coh>_joint_<source>.h5ad` (obsm `X_pca`, `X_pca_harmony`; obs incl. transferred labels and joint Leiden).
* `code/big_leiden_res.py` — cluster matching at Leiden resolution 0.5 / 1 / 2.5.
* `code/big_rare.py` — rare subsets within parent lineage (definitions of `integration_report_20261012/code/f_rare_parent.py`), marker AUROC, marker-called shares.
* Deck figures: `code/fig_*.py` → `figs/deck/`; deck: `code/build_deck_ko.py`.""")
code("""import json, numpy as np, pandas as pd
from IPython.display import Image, display
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
B = '/mnt/user-data/outputs/ilisi_floor_20261021'; RB = B + '/results/big'; FIG = B + '/figs/deck'""")
md("## 1. Gene count curve and cell shape")
code("""rows = []
for c in ('sepsis', 'covid'):
    G = pd.read_csv(f'{RB}/{c}_tx_gene.csv.gz'); ref = G[G.source == 'ref'].set_index('gene').mean_cp10k
    for s in ('realT1f', 'floor', 'prev4k'):
        d = G[G.source == s].set_index('gene'); r = d.mean_cp10k / ref.reindex(d.index)
        m = r.groupby(pd.cut(d.real_rank, [0, 500, 1500, 2000, 3000, 4000, 10**5], labels=['1-500', '501-1500', '1501-2000', '2001-3000', '3001-4000', '>4000']), observed=True).median()
        rows.append({'cohort': c, 'source': s, **m.round(3).to_dict(), 'genes_never_generated': int((d.mean_cp10k == 0).sum()), 'real_A_zero_genes': int((ref == 0).sum())})
pd.DataFrame(rows)""")
md("Median per-gene CP10k ratio to real half A by real rank bin. `realT1f` = real half B thinned to floor depth (the read-out control).")
code("""cell = pd.concat([pd.read_csv(f'{RB}/{c}_tx_cell.csv.gz').assign(cohort=c) for c in ('sepsis', 'covid')])
cell.groupby(['cohort', 'source'])[['umi', 'n_genes']].median().unstack(0).round(0)""")
code("display(Image(f'{FIG}/gene_count_curve.png', width=1100)); display(Image(f'{FIG}/cell_shape.png', width=1100))")
md("## 2. Transcriptome characteristics")
code("""pb = pd.concat([pd.read_csv(f'{RB}/{c}_tx_pb.csv').assign(cohort=c) for c in ('sepsis', 'covid')])
pb.pivot_table(index=['cohort', 'gene_set', 'source'], columns='six_type', values='pearson').round(3)""")
code("display(Image(f'{FIG}/pseudobulk_r.png', width=1100)); display(Image(f'{FIG}/marker_dot.png', width=1100))")
md("## 3. Harmony integration, iLISI, cluster matching, label transfer")
code("""S = pd.read_csv(f'{RB}/big_joint_summary.csv'); S.set_index(['cohort', 'source']).round(3).T""")
code("""LR = pd.concat([pd.read_csv(f'{RB}/{c}_{s}_leiden_res.csv') for c in ('sepsis', 'covid') for s in ('realU', 'floor', 'prev4k')])
LR.round(3)""")
code("""for c in ('sepsis', 'covid'): display(Image(f'{FIG}/umap_{c}.png', width=1100))
display(Image(f'{FIG}/ilisi_bar_deck.png', width=1100)); display(Image(f'{FIG}/cluster_matching.png', width=1000))""")
md("The 20,000-cell iLISI (main numbers) is in `results/ilisi_all.csv` and `ilisi_floor.ipynb`; the ~60,000-cell values above are the second, independent subsample.")
md("## 4. Rare populations")
code("""R = pd.read_csv(f'{RB}/rare_summary.csv')
R[['cohort', 'subset', 'parent', 'real B (author)', 'real B (transferred)', 'floor (transferred)', 'prev K4000 (transferred)', 'floor_vs_realB_author', 'prev_vs_realB_author',
   'floor_vs_realB_tr', 'prev_vs_realB_tr', 'n_floor', 'n_prev', 'auroc_floor', 'auroc_realB_author']].round(3)""")
md("""Shares are % of the parent lineage (pooled). Real half B is shown with author labels (truth) and with kNN-transferred labels (same route as generated cells = ceiling of the route).
Verdict on the slide: ratio to real B (author) within 0.5–2 = recovered, < 0.5 = under-produced, > 2 = over-produced, 0 = absent.""")
code("""cat = lambda v: 'absent' if v == 0 or np.isnan(v) else ('low' if v < 0.5 else ('over' if v > 2 else 'ok'))
pd.DataFrame({k: R.groupby('cohort')[f'{k}_vs_realB_author'].apply(lambda x: x.map(cat).value_counts()).unstack(fill_value=0).stack() for k in ('floor', 'prev')}).unstack(0)""")
code("""MC = pd.concat([pd.read_csv(f'{RB}/{c}_rare_marker_called.csv') for c in ('sepsis', 'covid')])
MC.pivot_table(index=['cohort', 'subset'], columns='source', values='pct_of_parent').round(2).join(MC.drop_duplicates(['cohort', 'subset']).set_index(['cohort', 'subset']).F1_realA.round(2))""")
md("""Marker-called shares (label-free inside the parent; threshold fitted on real A author labels) are shown for completeness only: F1 on real A is below 0.5 for most proliferating subsets, so they are not used for the verdict.""")
code("display(Image(f'{FIG}/rare_parent_share.png', width=1100))")
nbf.write(nb, '/mnt/user-data/outputs/ilisi_floor_20261021/floor_validation_deck.ipynb')
