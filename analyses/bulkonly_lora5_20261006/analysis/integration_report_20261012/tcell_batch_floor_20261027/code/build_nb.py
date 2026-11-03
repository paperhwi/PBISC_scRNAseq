"""Build tcell_batch.ipynb: per-T-subset expression differences between real A (reference), real B (real control) and floor-generated cells."""
import sys, nbformat as nbf
D, W = sys.argv[1], sys.argv[2]
nb = nbf.v4.new_notebook(); C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# T-cell subsets: how different is gene expression per batch? Real vs floor-generated (2026-10-27)

**Question.** In the integrated real + floor datasets, how much does gene expression of each T-cell subset differ by batch (data source)?

**Batches (sources).**
- **real A:** real held-out cells, patient half A, author fine labels. This is the reference that the floor cells were integrated with.
- **real B:** real held-out cells, the other patient half, labels transferred from real A by kNN (k = 30) after Harmony. This is the real-vs-real control: the size of a "batch" difference between two groups of real patients.
- **floor:** PBISC floor read-out (depth256-K4000-floor1e-4, draw0), same kNN label transfer from real A.

Inputs come from `ilisi_floor_20261021/data/big/` (~60,000 cells per source, raw counts on the panel) and the joint Harmony label transfer in `data/big_joint/`. Real B and floor use the same labelling route, so the real-B row is also the control for label-transfer effects.

**Method.**
- **Pseudobulk.** One pseudobulk per sample × subset. Each sums the raw counts of 20 randomly drawn cells, so every pseudobulk has the same cell number and noise level. Samples with fewer than 20 cells in the subset are dropped. A subset is kept if every source has at least 4 samples.
- **Genes.** Expressed genes only: mean CPM ≥ 20 in real A for that subset. Values are log2(CPM + 1).
- **Metrics per subset (source vs real A).**
  - (i) Pearson r of the mean profiles.
  - (ii) **Batch R²**: the share of between-sample variance explained by source. Floor samples are subsampled to match real A's sample count (20 repeats, mean), so R² is comparable across sources.
  - (iii) DEGs from a per-gene OLS model `log2CPM ~ source + disease group`, with |log2FC| ≥ 1 and BH q < 0.01.
  - (iv) **Model-specific DEGs**: floor-vs-A DEGs that are not also real-B-vs-A DEGs (q < 0.05 and |log2FC| ≥ 0.5).
- **Visuals.** Marker dot plots use cell-level log1p CP10k means and detected fractions.""")
code(f"""import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, matplotlib as mpl, matplotlib.pyplot as plt
from scipy import stats
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30); pd.set_option('display.max_rows', 80)
D = '{D}'; W = '{W}'
COH = ['sepsis', 'covid']; SRC = ['real A', 'real B', 'floor']; FILE = {{'real A': 'ref', 'real B': 'realU', 'floor': 'floor'}}
COL = {{'real A': '#3E4C59', 'real B': '#9AA5B1', 'floor': '#2B83BA'}}
TITLE = {{'sepsis': 'Sepsis (Ye 2026)', 'covid': 'COVID (Stephenson 2025)'}}
INK, INK2, MUT, GRID = '#1F2933', '#3E4C59', '#7B8794', '#E8ECF0'
mpl.rcParams.update({{'font.family': 'DejaVu Sans', 'font.size': 9.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.titlelocation': 'left', 'axes.titlesize': 10.5, 'svg.fonttype': 'none'}})
NCELL, MINS, MINCPM = 20, 4, 20
SEXG = {'XIST', 'TSIX', 'RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY', 'USP9Y', 'ZFY', 'TXLNGY', 'NLGN4Y', 'PRKY', 'RPS4Y2'}   # sex-linked: sample sex mix differs between sources, not a model effect
rng = np.random.default_rng(20261027)""")
code("""def load(c):
    ref_obs = ad.read_h5ad(f'{D}/big/{c}_ref.h5ad', backed='r').obs
    tsub = sorted(ref_obs.loc[ref_obs.lineage == 'T', 'fine'].astype(str).unique())
    out = {}
    for s in SRC:
        a = ad.read_h5ad(f'{D}/big/{c}_{FILE[s]}.h5ad')
        if s == 'real A':
            lab = a.obs['fine'].astype(str).to_numpy()
        else:
            J = ad.read_h5ad(f'{D}/big_joint/{c}_joint_{FILE[s]}.h5ad', backed='r').obs
            q = J[J.source == FILE[s]]
            assert len(q) == a.n_obs and (q['sample'].astype(str).to_numpy() == a.obs['sample'].astype(str).to_numpy()).all()
            lab = q['fine_tr'].astype(str).to_numpy()          # kNN transfer from real A (k = 30) on joint Harmony
        m = np.isin(lab, tsub)
        X = sp.csr_matrix(a.X)[np.where(m)[0]]
        obs = a.obs.iloc[np.where(m)[0]][['sample', 'patient', 'group']].astype(str).copy(); obs['subset'] = lab[m]; obs['source'] = s
        out[s] = (X, obs)
    return out, tsub, a.var_names.to_numpy()
DATA, TSUB, GENES = {}, {}, {}
for c in COH:
    DATA[c], TSUB[c], GENES[c] = load(c)
    print(c, {s: DATA[c][s][0].shape[0] for s in SRC}, 'T subsets in real A:', len(TSUB[c]))""")
md("""## 1. Pseudobulks: 20 cells per sample × subset""")
code("""PB = {}
for c in COH:
    rows, mats = [], []
    for s in SRC:
        X, obs = DATA[c][s]
        for (smp, sub), idx in obs.groupby(['sample', 'subset']).indices.items():
            if len(idx) < NCELL: continue
            take = rng.choice(idx, NCELL, replace=False)
            mats.append(np.asarray(X[take].sum(0)).ravel())
            rows.append(dict(source=s, sample=smp, subset=sub, group=obs['group'].iloc[idx[0]], n_cells_available=len(idx)))
    P = ad.AnnData(np.vstack(mats).astype(np.float32), obs=pd.DataFrame(rows), var=pd.DataFrame(index=GENES[c]))
    P.obs.index = [f'{c}|{r.source}|{r.sample}|{r.subset}' for r in P.obs.itertuples()]
    P.uns['description'] = f'Pseudobulk raw counts: sum of {NCELL} random cells per sample x T subset; seed 20261027; sources real A (author labels), real B and floor (kNN-transferred labels).'
    PB[c] = P
    P.write_h5ad(f'{W}/data/{c}_tcell_pseudobulk_20cells.h5ad', compression='gzip')
    cnt = P.obs.groupby(['subset', 'source']).size().unstack(fill_value=0)[SRC]
    print(c); print(cnt.to_string())""")
md("""## 2. Per-subset batch metrics""")
code("""def logcpm(M): return np.log2(M / np.maximum(M.sum(1, keepdims=True), 1) * 1e6 + 1)
def ols_source(Y, src_ind, group):
    # per-gene OLS: log2CPM ~ intercept + source + disease-group dummies -> source coefficient (log2FC), two-sided p
    G = pd.get_dummies(pd.Series(group), drop_first=True).to_numpy(float)
    G = G[:, G.std(0) > 0] if G.size else np.zeros((len(Y), 0))
    Xd = np.column_stack([np.ones(len(Y)), src_ind.astype(float), G])
    XtXi = np.linalg.pinv(Xd.T @ Xd); beta = XtXi @ Xd.T @ Y; res = Y - Xd @ beta
    dof = max(len(Y) - np.linalg.matrix_rank(Xd), 1); s2 = (res ** 2).sum(0) / dof
    t = beta[1] / np.sqrt(s2 * XtXi[1, 1] + 1e-12); p = 2 * stats.t.sf(np.abs(t), dof)
    return beta[1], p
def bh(p):
    o = np.argsort(p); q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty_like(q); out[o] = np.minimum(q, 1); return out
def r2_source(Ya, Yx, reps=20):
    vals = []
    for _ in range(reps if len(Yx) > len(Ya) else 1):
        Yb = Yx[rng.choice(len(Yx), len(Ya), replace=False)] if len(Yx) > len(Ya) else Yx
        Y = np.vstack([Ya, Yb]); mu = Y.mean(0); tot = ((Y - mu) ** 2).sum()
        btw = len(Ya) * ((Ya.mean(0) - mu) ** 2).sum() + len(Yb) * ((Yb.mean(0) - mu) ** 2).sum()
        vals.append(btw / tot)
    return float(np.mean(vals))
SUM, DEG = [], []
for c in COH:
    P = PB[c]; o = P.obs
    for sub in TSUB[c]:
        n = {s: int(((o.subset == sub) & (o.source == s)).sum()) for s in SRC}
        if min(n.values()) < MINS: continue
        L = {s: logcpm(P.X[((o.subset == sub) & (o.source == s)).to_numpy()]) for s in SRC}
        cpmA = 2 ** L['real A'] - 1; keep = (cpmA.mean(0) >= MINCPM) & ~np.isin(GENES[c], list(SEXG))
        L = {s: v[:, keep] for s, v in L.items()}; g = GENES[c][keep]
        grp = {s: o.loc[(o.subset == sub) & (o.source == s), 'group'].to_numpy() for s in SRC}
        row = dict(cohort=c, subset=sub, **{f'n_samples {s}': n[s] for s in SRC}, n_genes=int(keep.sum()))
        res = {}
        for s in ['real B', 'floor']:
            Y = np.vstack([L['real A'], L[s]]); ind = np.r_[np.zeros(len(L['real A'])), np.ones(len(L[s]))]
            lfc, p = ols_source(Y, ind, np.r_[grp['real A'], grp[s]]); q = bh(p); res[s] = (lfc, q)
            row[f'r mean vs A | {s}'] = float(np.corrcoef(L['real A'].mean(0), L[s].mean(0))[0, 1])
            row[f'batch R2 vs A | {s}'] = r2_source(L['real A'], L[s])
            row[f'DEG vs A | {s}'] = int(((np.abs(lfc) >= 1) & (q < 0.01)).sum())
        lB, qB = res['real B']; lF, qF = res['floor']
        bsig = (qB < 0.05) & (np.abs(lB) >= 0.5)
        ms = (np.abs(lF) >= 1) & (qF < 0.01) & ~bsig
        row['model-specific DEG'] = int(ms.sum())
        row['model-specific up'] = int((ms & (lF > 0)).sum()); row['model-specific down'] = int((ms & (lF < 0)).sum())
        # power-matched: floor subsampled to the real-B sample count (20 repeats, median), same model
        nB = len(L['real B']); pm, pms = [], []
        for _ in range(20 if len(L['floor']) > nB else 1):
            ii = rng.choice(len(L['floor']), nB, replace=False) if len(L['floor']) > nB else np.arange(len(L['floor']))
            Y = np.vstack([L['real A'], L['floor'][ii]]); ind = np.r_[np.zeros(len(L['real A'])), np.ones(len(ii))]
            l2, p2 = ols_source(Y, ind, np.r_[grp['real A'], grp['floor'][ii]]); sig = (np.abs(l2) >= 1) & (bh(p2) < 0.01)
            pm.append(sig.sum()); pms.append((sig & ~bsig).sum())
        row['DEG vs A | floor, n matched to real B'] = float(np.median(pm)); row['model-specific DEG, n matched'] = float(np.median(pms))
        d = pd.DataFrame({'cohort': c, 'subset': sub, 'gene': g, 'log2FC floor vs A': lF, 'q floor vs A': qF, 'log2FC realB vs A': lB, 'q realB vs A': qB,
                          'mean log2CPM A': L['real A'].mean(0), 'model_specific': ms})
        DEG.append(d)
        top = d[d.model_specific].sort_values('log2FC floor vs A')
        row['top up in floor'] = ', '.join(top[top['log2FC floor vs A'] > 0].tail(5).gene[::-1]); row['top down in floor'] = ', '.join(top[top['log2FC floor vs A'] < 0].head(5).gene)
        SUM.append(row)
S = pd.DataFrame(SUM); DG = pd.concat(DEG)
S.to_csv(f'{W}/results/tcell_batch_summary.csv', index=False); DG.to_csv(f'{W}/results/tcell_deg_floor_vs_realA_all_genes.csv.gz', index=False)
S[['cohort', 'subset', 'n_samples real A', 'n_samples real B', 'n_samples floor', 'n_genes', 'r mean vs A | real B', 'r mean vs A | floor',
   'batch R2 vs A | real B', 'batch R2 vs A | floor', 'DEG vs A | real B', 'DEG vs A | floor', 'DEG vs A | floor, n matched to real B', 'model-specific DEG', 'model-specific DEG, n matched']].round(3)""")
md("""**Table 1.** For each subset: mean-profile r and batch R² against real A for the real control (real B) and for floor, DEG counts (|log2FC| ≥ 1, q < 0.01, disease group as covariate), and model-specific DEGs (floor DEGs that are not real-B DEGs). Floor usually has many more samples than real B and so more power; the 'n matched' columns subsample floor to the real-B sample count (median of 20 repeats) and are the fair comparison. Sex-linked genes (XIST, chrY) are excluded because the sex mix differs between sample sets. Read real B as the size of an ordinary patient-group batch effect and floor as that plus the model effect.""")
code("""S[['cohort', 'subset', 'model-specific up', 'model-specific down', 'top up in floor', 'top down in floor']]""")

md("""## 3. Figures""")
code("""fig, axs = plt.subplots(2, 3, figsize=(15, 10.5), gridspec_kw={'width_ratios': [1.1, 1.1, 1.2]})
for i, c in enumerate(COH):
    s = S[S.cohort == c].sort_values('batch R2 vs A | floor'); y = np.arange(len(s))
    ax = axs[i, 0]
    ax.barh(y + 0.2, s['batch R2 vs A | real B'], 0.38, color=COL['real B'], label='real B vs real A (patient batch)')
    ax.barh(y - 0.2, s['batch R2 vs A | floor'], 0.38, color=COL['floor'], label='floor vs real A')
    ax.set_yticks(y, s.subset); ax.set_xlabel('batch R² (share of between-sample variance explained by source)'); ax.set_xlim(0, 1)
    ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True); ax.set_title(f'{TITLE[c]} · batch R²')
    for yy, v in zip(y - 0.2, s['batch R2 vs A | floor']): ax.text(v + 0.01, yy, f'{v:.2f}', va='center', fontsize=8, color=INK)
    if i == 0: ax.legend(frameon=False, fontsize=8.5, loc='lower right')
    ax = axs[i, 1]
    ax.barh(y + 0.2, s['DEG vs A | real B'], 0.38, color=COL['real B'], label='real B vs real A')
    ax.barh(y - 0.2, s['DEG vs A | floor, n matched to real B'], 0.38, color=COL['floor'], label='floor vs real A (samples matched to real B)')
    ax.set_yticks(y, [''] * len(y)); ax.set_xlabel('DEGs (|log2FC| ≥ 1, q < 0.01), floor samples matched to real B'); ax.xaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.set_title('DEG count');
    for yy, v, m in zip(y - 0.2, s['DEG vs A | floor, n matched to real B'], s['model-specific DEG, n matched']): ax.text(v + 2, yy, f'{v:.0f} ({m:.0f} model-specific)', va='center', fontsize=8, color=INK)
    if i == 0: ax.legend(frameon=False, fontsize=8.5, loc='lower right')
    ax = axs[i, 2]
    ax.scatter(s['r mean vs A | real B'], s['r mean vs A | floor'], s=36, color=COL['floor'], zorder=3)
    for _, r in s.iterrows(): ax.annotate(r.subset.replace('T_', '').replace('-GPR183', ''), (r['r mean vs A | real B'], r['r mean vs A | floor']), fontsize=7.5, color=INK2, xytext=(3, 2), textcoords='offset points')
    lo = min(s['r mean vs A | floor'].min(), s['r mean vs A | real B'].min()) - 0.01
    ax.plot([lo, 1], [lo, 1], color=MUT, lw=0.8, ls='--'); ax.set_xlim(lo, 1.0); ax.set_ylim(lo, 1.0)
    ax.set_xlabel('r of mean profile, real B vs real A'); ax.set_ylabel('r of mean profile, floor vs real A'); ax.set_title('mean-profile agreement (dashed = equal)')
fig.tight_layout(); [fig.savefig(f'{W}/figs/tcell_batch_metrics.{e}', dpi=200) for e in ('png', 'svg', 'pdf')]; plt.show()""")
md("""**Figure 1.** Per T subset: batch R² (left), DEG counts (middle; label = floor DEGs and how many are model-specific) and mean-profile r (right) against real A. Grey = real B, the effect of moving to another group of real patients; blue = floor.""")
code("""from sklearn.decomposition import PCA
for c in COH:
    s = S[S.cohort == c].sort_values('n_samples real A', ascending=False); subs = s.subset.tolist(); nc = 5; nr = int(np.ceil(len(subs) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.1 * nc, 2.9 * nr)); axs = np.atleast_2d(axs)
    P = PB[c]; o = P.obs
    for k, sub in enumerate(subs):
        ax = axs.flat[k]; m = (o.subset == sub).to_numpy(); L = logcpm(P.X[m]); src = o.source[m].to_numpy()
        keep = (2 ** L[src == 'real A'] - 1).mean(0) >= MINCPM; Z = L[:, keep]; Z = Z - Z[src == 'real A'].mean(0)
        pc = PCA(2, random_state=0).fit(Z); E = pc.transform(Z)
        for sname in ['floor', 'real B', 'real A']:
            mm = src == sname; ax.scatter(E[mm, 0], E[mm, 1], s=10 if sname == 'floor' else 16, color=COL[sname], alpha=0.55 if sname == 'floor' else 0.9, lw=0, label=sname)
        r2 = s.set_index('subset').loc[sub]
        ax.set_title(f"{sub}\\nR² B {r2['batch R2 vs A | real B']:.2f} · floor {r2['batch R2 vs A | floor']:.2f}", fontsize=8.5)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(f'PC1 {pc.explained_variance_ratio_[0]*100:.0f}%', fontsize=8); ax.set_ylabel(f'PC2 {pc.explained_variance_ratio_[1]*100:.0f}%', fontsize=8)
    for k in range(len(subs), nr * nc): axs.flat[k].axis('off')
    axs.flat[0].legend(frameon=False, fontsize=8, loc='best', markerscale=1.5)
    fig.suptitle(f'{TITLE[c]} · pseudobulk PCA per T subset (one point = one sample, 20 cells)', x=0.01, ha='left', fontsize=11.5, color=INK)
    fig.tight_layout(); [fig.savefig(f'{W}/figs/tcell_pseudobulk_pca_{c}.{e}', dpi=200) for e in ('png', 'svg', 'pdf')]; plt.show()""")
md("""**Figure 2.** Pseudobulk PCA within each subset (expressed genes, centred on real A). Where floor (blue) sits apart from real A and real B (dark and light grey) along PC1, source is the main axis of variation in that subset.""")
code("""for c in COH:
    d = DG[(DG.cohort == c) & DG.model_specific]
    subs = S[S.cohort == c].sort_values('n_samples real A', ascending=False).subset.tolist()
    top = (d.assign(a=d['log2FC floor vs A'].abs()).sort_values('a', ascending=False).groupby('subset').head(4).gene.unique().tolist())
    freq = d.groupby('gene').size().sort_values(ascending=False); top = list(dict.fromkeys(freq.index[:14].tolist() + top))[:30]
    F = DG[(DG.cohort == c) & DG.gene.isin(top)].pivot_table(index='gene', columns='subset', values='log2FC floor vs A').reindex(index=top, columns=subs)
    B = DG[(DG.cohort == c) & DG.gene.isin(top)].pivot_table(index='gene', columns='subset', values='log2FC realB vs A').reindex(index=top, columns=subs)
    fig, axs = plt.subplots(1, 2, figsize=(14, 0.32 * len(top) + 2.6), sharey=True)
    for ax, M, t in [(axs[0], F, 'floor vs real A'), (axs[1], B, 'real B vs real A (control)')]:
        im = ax.imshow(M.to_numpy(), cmap='RdBu_r', vmin=-3, vmax=3, aspect='auto')
        ax.set_xticks(range(len(subs)), [x.replace('T_', '') for x in subs], rotation=60, ha='right', fontsize=8); ax.set_yticks(range(len(top)), top, fontsize=8.5)
        ax.set_title(t); ax.spines[:].set_visible(False)
    fig.colorbar(im, ax=axs, shrink=0.6, label='log2FC vs real A (grey = not expressed)')
    fig.suptitle(f'{TITLE[c]} · most frequent and strongest model-specific genes across T subsets', x=0.01, ha='left', fontsize=11.5, color=INK)
    [fig.savefig(f'{W}/figs/tcell_model_specific_heatmap_{c}.{e}', dpi=200, bbox_inches='tight') for e in ('png', 'svg', 'pdf')]; plt.show()
    print(c, 'genes model-specific in >= half of the subsets:', ', '.join(freq[freq >= len(subs) / 2].index[:40]))""")
md("""**Figure 3.** log2FC against real A for model-specific genes: those shared across many subsets plus the strongest per subset. Left = floor, right = the real-B control on the same genes and colour scale. Genes red or blue in many floor columns but pale on the right are model-wide biases, not subset biology.""")
code("""MK = {'sepsis': ['CD3E', 'CD4', 'CD8A', 'CCR7', 'SELL', 'LEF1', 'TCF7', 'IL7R', 'GPR183', 'NR4A2', 'ANXA2', 'GZMK', 'GZMB', 'GNLY', 'PRF1', 'NKG7', 'CCL4', 'IFNG', 'FOXP3', 'IL2RA', 'KLRB1', 'SLC4A10', 'HLA-DRA', 'MKI67', 'CX3CR1', 'CXCR4'],
      'covid': ['CD3E', 'CD4', 'CD8A', 'CCR7', 'SELL', 'LEF1', 'TCF7', 'IL7R', 'GZMK', 'GZMB', 'GNLY', 'PRF1', 'NKG7', 'FOXP3', 'IL2RA', 'KLRB1', 'SLC4A10', 'TRDC', 'IFI6', 'ISG15', 'MX1', 'MKI67', 'CX3CR1', 'HLA-DRA']}
for c in COH:
    genes = [g for g in MK[c] if g in set(GENES[c])]; gi = pd.Index(GENES[c]).get_indexer(genes)
    subs = S[S.cohort == c].sort_values('n_samples real A', ascending=False).subset.tolist()
    rows = []
    for sub in subs:
        for s in SRC:
            X, obs = DATA[c][s]; m = np.where(obs.subset.to_numpy() == sub)[0]; Xs = X[m]
            tot = np.asarray(Xs.sum(1)).ravel(); Y = Xs[:, gi].toarray() / np.maximum(tot, 1)[:, None] * 1e4
            rows += [dict(subset=sub, source=s, gene=g, mean=np.log1p(Y[:, j]).mean(), det=(Y[:, j] > 0).mean()) for j, g in enumerate(genes)]
    T = pd.DataFrame(rows); T.to_csv(f'{W}/results/{c}_tcell_marker_dot.csv', index=False)
    ylab = [f'{sub.replace("T_", "")} · {s}' for sub in subs for s in SRC]
    fig, ax = plt.subplots(figsize=(0.42 * len(genes) + 3.2, 0.24 * len(ylab) + 1.8))
    vmax = T.groupby('gene')['mean'].max()
    for i, (sub, s) in enumerate([(a, b) for a in subs for b in SRC]):
        t = T[(T.subset == sub) & (T.source == s)].set_index('gene').loc[genes]
        ax.scatter(np.arange(len(genes)), np.full(len(genes), i), s=t.det.to_numpy() * 70 + 0.5, c=(t['mean'] / vmax.loc[genes]).to_numpy(), cmap='Blues', vmin=0, vmax=1, edgecolors=COL[s], linewidths=0.7)
        if s == 'floor' and i < len(ylab) - 1: ax.axhline(i + 0.5, color=GRID, lw=0.8)
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=8.5); ax.set_yticks(range(len(ylab)), ylab, fontsize=7.5); ax.invert_yaxis()
    for tl, s in zip(ax.get_yticklabels(), [b for a in subs for b in SRC]): tl.set_color(COL[s] if s != 'real B' else '#6B7785')
    ax.set_title(f'{TITLE[c]} · T-subset markers by source (size = % detected, shade = mean log1p CP10k scaled per gene)')
    fig.tight_layout(); [fig.savefig(f'{W}/figs/tcell_marker_dot_{c}.{e}', dpi=200) for e in ('png', 'svg', 'pdf')]; plt.show()""")
md("""**Figure 4.** Marker dot plot: for each subset, three rows in the order real A, real B, floor. Large departures of the floor row from the two real rows mark the subset-defining programmes that the model reproduces poorly.""")
md("""## 4. Summary""")
code("""Q = S.assign(R2_excess=S['batch R2 vs A | floor'] - S['batch R2 vs A | real B'])
print(Q.groupby('cohort')[['batch R2 vs A | real B', 'batch R2 vs A | floor', 'DEG vs A | real B', 'DEG vs A | floor, n matched to real B', 'model-specific DEG, n matched', 'r mean vs A | real B', 'r mean vs A | floor']].median().round(3).T)
Q.sort_values(['cohort', 'R2_excess'], ascending=[True, False])[['cohort', 'subset', 'batch R2 vs A | real B', 'batch R2 vs A | floor', 'R2_excess', 'DEG vs A | floor, n matched to real B', 'model-specific DEG', 'top up in floor', 'top down in floor']].round(3)""")

nbf.write(nb, f'{W}/tcell_batch.ipynb')
