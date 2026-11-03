import sys, nbformat as nbf
W = sys.argv[1]; nb = nbf.v4.new_notebook(); C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Floor read-out vs real: iLISI after Harmony (2026-10-21)

**Question.** After batch correction against real cells, how well do the new floor-read-out cells (`depth256-K4000-floor1e-4`, draw0) mix with real cells, compared with the previous K4000 cells and with real-vs-real controls?

**Pipeline** (same as `integration_report_20261012/code/i_hvg_pca.py`, `code/ilisi_lib.py` here): the real held-out cells (`cell_transcriptome_validation_20261014`; sepsis adults 40 patients, COVID 35 patients) are split **by patient** into a reference half A and a control half B. 2,000 `seurat_v3` HVGs and a 30-PC PCA are fitted on the reference; every source is projected into this PCA space; **pairwise Harmony** (batch = source) aligns reference and source; iLISI = inverse Simpson index of the source label in a perplexity-30 Gaussian neighbourhood (k = 90), median over 10,000 random cells of the pooled set (mean of three query draws; per-lineage values on ≤ 5,000 cells). Higher = better mixing, maximum 2.

| source | what it is |
|---|---|
| `realU` | real half B, untouched — upper control (different patients from the reference) |
| `realT1` | real half B thinned to the depth of floor cells of the same six-type |
| `floor` | new floor read-out file, draw0 (sepsis adults; COVID without the 5 adapter-training samples) |
| `prev4k` | previous K4000 Harmony-set cells (Poisson read-out, 5 draws pooled), same pipeline and sizes |

Every source and the reference have **20,000 cells** (equal number per sample), the same 5,961 (sepsis) / 5,745 (COVID) panel genes. **V0** = reference HVG + reference PCA (the report's main pipeline); **V1** = HVG chosen jointly (batch_key = source), PCA on the reference. Numbers are therefore like-for-like within this notebook, **not** comparable with the report's pooled values (1.15 / 1.06, 1.74 ceiling), which used other cell numbers.""")
code("""import json, glob, os, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)
W = '%s'; rows = []
for f in sorted(glob.glob(f'{W}/results/*_ilisi_*.json')):
    c = os.path.basename(f).split('_')[0]
    for v, d in json.load(open(f)).items():
        for s, x in d.items(): rows.append(dict(cohort=c, variant=v, source=s, **x))
R = pd.DataFrame(rows); ORD = ['realU', 'realT1', 'floor', 'prev4k']; R['source'] = pd.Categorical(R.source, ORD, ordered=True)
R = R.sort_values(['cohort', 'variant', 'source']).reset_index(drop=True); R.to_csv(f'{W}/results/ilisi_all.csv', index=False)
print(len(R), 'runs;', 'n per source:', sorted(set(R.n_source)), 'ref:', sorted(set(R.n_ref)))
main = R[['cohort', 'variant', 'source', 'iLISI_pre', 'iLISI_post', 'iLISI_post_sd_query', 'own_share_post']].copy()
main['share_of_real_ceiling_post'] = main.apply(lambda r: (r.iLISI_post - 1) / (R[(R.cohort == r.cohort) & (R.variant == r.variant) & (R.source == 'realU')].iLISI_post.iloc[0] - 1), axis=1)
main.round(3)""" % W)
md("""**Table 1.** iLISI before (`pre`, PCA projection) and after Harmony (`post`). `own_share_post` = share of the 90 neighbours that come from the cell's own source (0.5 = perfect mixing). `share_of_real_ceiling_post` = (iLISI − 1) / (real half B iLISI − 1).""")
code("""L = R[R.variant == 'V0'][['cohort', 'source', 'iLISI_post_T', 'iLISI_post_NK', 'iLISI_post_B', 'iLISI_post_Myeloid']].copy(); L.columns = ['cohort', 'source', 'T', 'NK', 'B', 'Myeloid']
L.to_csv(f'{W}/results/ilisi_post_by_lineage_V0.csv', index=False); L.round(3)""")
md("""**Table 2.** Post-Harmony iLISI within coarse lineage (V0). Lineage = frozen scANVI six-type for generated cells, author labels mapped to the same six types for real cells; Other (plasma, platelets etc.) not shown.""")
code("""do = pd.read_csv(f'{W}/results/covid_donor_overlap_ilisi.csv'); do.round(3)""")
md("""**Table 3 (COVID, donor-overlap check).** 9 of the 17 reference patients and 8 of the 18 control patients also appear as donors in the generated set (bulk-conditioned). Per-cell iLISI of generated cells (query = generated cells only, V0 post) by donor group: generated cells from donors that have real cells in the reference are not better mixed than those from other donors (floor 1.652 vs 1.628), so shared donors do not explain the result. Sepsis generated donors (121) have no overlap with the 40 real patients. Generated-side-only values are lower for `prev4k` (1.09) than the pooled-query Table 1 (1.37) because the previous cells form their own neighbourhoods.""")
code("""COL = {'realU': '#1f1f1e', 'realT1': '#8a8984', 'floor': '#1b8f6a', 'prev4k': '#eb6834'}; LB = {'realU': 'real half B (control)', 'realT1': 'real half B, depth-matched', 'floor': 'floor read-out (draw0)', 'prev4k': 'previous K4000'}
INK2, MUTED, GRID = '#5f5e5a', '#8a8984', '#e6e5e0'
mpl.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={'width_ratios': [1.1, 1.1, 1.5]})
for ax, c, t in zip(axs[:2], ['sepsis', 'covid'], ['Sepsis (Ye 2026)', 'COVID (Imperial ESKD)']):
    d = R[(R.cohort == c) & (R.variant == 'V0')].set_index('source'); x = np.arange(len(ORD))
    ax.bar(x - 0.2, d.loc[ORD, 'iLISI_pre'], 0.38, color=[COL[s] for s in ORD], alpha=0.35, label='before Harmony'); ax.bar(x + 0.2, d.loc[ORD, 'iLISI_post'], 0.38, color=[COL[s] for s in ORD], label='after Harmony')
    for i, s in enumerate(ORD): ax.text(i + 0.2, d.loc[s, 'iLISI_post'] + 0.02, f"{d.loc[s, 'iLISI_post']:.2f}", ha='center', fontsize=8, color='#1f1f1e')
    ax.set_xticks(x, ['real\\n(half B)', 'real, depth-\\nmatched', 'floor\\nread-out', 'previous\\nK4000'], fontsize=8); ax.set_ylim(1, 2); ax.axhline(2, color=MUTED, lw=0.6, ls=':'); ax.set_ylabel('iLISI (source label; 2 = perfect mixing)'); ax.set_title(t, loc='left', fontsize=10)
    ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
axs[0].legend(handles=[mpl.patches.Patch(color=MUTED, alpha=0.35, label='before Harmony'), mpl.patches.Patch(color=MUTED, label='after Harmony')], frameon=False, fontsize=8, loc='upper left')
ax = axs[2]; LIN = ['T', 'NK', 'B', 'Myeloid']; k = 0; xt = []; xl = []
for c, t in [('sepsis', 'Sepsis'), ('covid', 'COVID')]:
    for j, ln in enumerate(LIN):
        for i, s in enumerate(['realU', 'floor', 'prev4k']):
            v = L[(L.cohort == c) & (L.source == s)][ln].iloc[0]; ax.bar(k + i * 0.26, v, 0.24, color=COL[s], label=LB[s] if (k == 0) else None)
        xt.append(k + 0.26); xl.append(f'{ln}\\n{t}'); k += 1
    k += 0.5
ax.set_xticks(xt, xl, fontsize=8); ax.set_ylim(1, 2); ax.set_ylabel('iLISI after Harmony (V0), by lineage'); ax.set_title('By lineage', loc='left', fontsize=10); ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=8, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.0))
fig.tight_layout(); os.makedirs(f'{W}/figs', exist_ok=True)
for e in ('png', 'svg', 'pdf'): fig.savefig(f'{W}/figs/ilisi_floor.{e}', dpi=200 if e == 'png' else None)
plt.show()""")
import os
cf = os.path.join(W, 'code', 'conclusion.md')
if os.path.exists(cf): md(open(cf, encoding='utf-8').read())
nbf.write(nb, f'{W}/ilisi_floor.ipynb')
