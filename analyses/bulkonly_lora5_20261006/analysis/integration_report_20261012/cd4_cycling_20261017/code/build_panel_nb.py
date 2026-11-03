import nbformat as nbf
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# COVID T CD4 cycling — figures for the deck panel (2026-10-18)
Companion to `cd4_cycling.ipynb` (same folder, same `_extract/` inputs). Builds the panel figures that mirror the sepsis macrophage section:
F4 T/NK UMAPs with labels and cycle-gene counts; F5 where real T CD4 cycling cells go in label transfer; F6 how many generated cells get the label;
F7 dot plot by label group; F8 label-free sub-cluster heatmaps; F9 panel coverage of the real markers.
Display sets: real = held-out cells (all T CD4 cycling + random 12k CD4 T, 4k CD8 T, 4k NK); generated = all AIFI "Proliferating T/NK" + random 12k T, 4k NK.""")
code("""import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, json, os, warnings, collections
import matplotlib.pyplot as plt, matplotlib as mpl
warnings.filterwarnings('ignore'); os.makedirs('figs', exist_ok=True); os.makedirs('results', exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
GS = json.load(open('code/gene_set.json'))
C14 = ['MKI67','TYMS','CENPF','TOP2A','BIRC5','CENPM','TPX2','PCLAF','CDC20','MYBL2','CENPE','UBE2C','ASPM','DLGAP5']
def load_obs(s): o = np.load(f'_extract/{s}_obs.npz', allow_pickle=True); return {k: o[k] for k in o.files}
def rows(s):
    r = np.load(f'_extract/{s}_rows.npz', allow_pickle=True)
    X = sp.csr_matrix((r['data'], r['indices'], r['indptr']), shape=tuple(r['shape'])); var = pd.Index(r['var']).astype(str)
    keep = ~var.duplicated(); return X[:, np.where(keep)[0]].tocsr(), r['rows'], var[keep]
REAL_MAP = {'T CD4 Naive': 'CD4 T', 'T CD4 Naive IFN stim': 'CD4 T', 'T CD4 Helper': 'CD4 T', 'T CD4 Helper IFN stim': 'CD4 T', 'T CD4 CTL': 'CD4 T',
            'T reg': 'Treg', 'T CD8 CTL': 'CD8 T', 'T CD8 CM': 'CD8 T', 'T CD8 Naive': 'CD8 T', 'T CD8 CTL IFN stim': 'CD8 T',
            'NK': 'NK', 'NK CD56 bright': 'NK', 'NK IFN stim': 'NK', 'T CD4 cycling': 'T CD4 cycling'}
GEN_MAP = {'Naive CD4 T cell': 'CD4 T', 'Memory CD4 T cell': 'CD4 T', 'Treg': 'Treg', 'Naive CD8 T cell': 'CD8 T', 'Memory CD8 T cell': 'CD8 T', 'CD8aa': 'CD8 T',
           'CD56dim NK cell': 'NK', 'CD56bright NK cell': 'NK', 'Proliferating T cell': 'Proliferating T', 'Proliferating NK cell': 'Proliferating NK',
           'MAIT': 'other T', 'gdT': 'other T', 'DN T cell': 'other T'}
def make(src):
    X, rw, var = rows(src); o = load_obs(src)
    if src == 'real_heldout':
        fine = o['fine'][rw].astype(str); maj = o['major'][rw].astype(str); lab = pd.Series(maj).map(REAL_MAP).fillna('other').values; fine_lab = fine
    else:
        l2 = o['aifi_l2'][rw].astype(str); lab = pd.Series(l2).map(GEN_MAP).fillna('other').values; fine_lab = l2
    A = ad.AnnData(X.astype(np.float32), obs=pd.DataFrame({'label': lab, 'fine': fine_lab, 'group': o['group'][rw].astype(str), 'sample': o['sample_id'][rw].astype(str)}), var=pd.DataFrame(index=var))
    A = A[A.obs.label != 'other'].copy()
    A.obs['n_cycle'] = np.asarray((A[:, [g for g in C14 if g in A.var_names]].X > 0).sum(1)).ravel()
    A.obs['cycling'] = A.obs.n_cycle >= 2
    return A
D = {'Real (held-out)': make('real_heldout'), 'PBISC K = 1,000': make('gen_k1000'), 'PBISC K = 4,000': make('gen_k4000')}
for k, A in D.items(): print(k, A.shape, A.obs.label.value_counts().to_dict(), 'cycling', int(A.obs.cycling.sum()))""")
md("## F4. UMAP of T/NK cells per dataset: labels and cycle-gene counts")
code("""for k, A in D.items():
    A.layers['counts'] = A.X.copy()
    sc.pp.highly_variable_genes(A, flavor='seurat_v3', n_top_genes=2000, layer='counts')
    sc.pp.normalize_total(A, 1e4); sc.pp.log1p(A); A.raw = A
    B = A[:, A.var.highly_variable].copy(); sc.pp.scale(B, max_value=10); sc.tl.pca(B, n_comps=30, random_state=0)
    A.obsm['X_pca'] = B.obsm['X_pca']; sc.pp.neighbors(A, n_neighbors=15, random_state=0); sc.tl.umap(A, random_state=0, min_dist=0.3)
    sc.tl.leiden(A, resolution=1.0, random_state=0, key_added='leiden', flavor='igraph', n_iterations=2, directed=False)
PAL = {'CD4 T': '#2A78D6', 'Treg': '#7FB2EA', 'CD8 T': '#EB6834', 'NK': '#1BAF7A', 'T CD4 cycling': '#C4373A', 'Proliferating T': '#C4373A', 'Proliferating NK': '#7A5CC2', 'other T': '#B9BEC7'}
fig, ax = plt.subplots(2, 3, figsize=(17, 10.5))
for j, (k, A) in enumerate(D.items()):
    U = A.obsm['X_umap']; lab = A.obs.label.values
    order = [l for l in PAL if l in set(lab)]
    for l in order:
        m = lab == l; big = l in ('T CD4 cycling', 'Proliferating T', 'Proliferating NK')
        ax[0, j].scatter(U[m, 0], U[m, 1], s=6 if big else 2, c=PAL[l], label=f'{l} ({m.sum():,})', alpha=0.9 if big else 0.5, lw=0, zorder=3 if big else 1)
    ax[0, j].set_title(f'{k}: label', fontsize=12); ax[0, j].legend(markerscale=4, frameon=False, fontsize=9, loc='best')
    n = A.obs.n_cycle.clip(upper=3).values; cmap = {0: '#E3E6EB', 1: '#F2C1A8', 2: '#E06A3A', 3: '#8E1B10'}
    for v in [0, 1, 2, 3]:
        m = n == v; ax[1, j].scatter(U[m, 0], U[m, 1], s=2 if v < 2 else 7, c=cmap[v], lw=0, zorder=v, label=f"{'≥3' if v == 3 else v} cycle genes ({m.sum():,})")
    ax[1, j].set_title(f'{k}: number of 14 cycle genes detected', fontsize=12); ax[1, j].legend(markerscale=4, frameon=False, fontsize=9, loc='best')
    for a in ax[:, j]: a.set_xticks([]); a.set_yticks([]); a.set_xlabel('UMAP 1'); a.set_ylabel('UMAP 2')
plt.tight_layout(); plt.savefig('figs/f4_tnk_umap.png', dpi=170); plt.show()""")
md("## F5. Where real T CD4 cycling cells go in label transfer (PBISC read-out of real cells, patient-disjoint)")
code("""cc = pd.concat([pd.read_csv(f'../../paired_patient_validation_20261014/results/patient_disjoint_fold{f}_confusion_counts.csv.gz') for f in [0, 1]])
cc = cc[(cc.source == 'realT2') & (cc.level == 'fine')]
dest = {}; rec = []
for lab in ['T CD4 cycling', 'B_ASC_dividing']:
    for r in ['R1', 'R2', 'R3']:
        t = cc[(cc.true == lab) & (cc.route == r)]; tot = t.n.sum(); hit = t[t.pred == lab].n.sum(); called = cc[(cc.pred == lab) & (cc.route == r)].n.sum()
        rec.append({'label': lab, 'route': r, 'true cells': int(tot), 'recall': hit / tot if tot else np.nan, 'precision': hit / called if called else np.nan})
        if lab == 'T CD4 cycling': dest[r] = t.groupby('pred').n.sum() / tot
rec = pd.DataFrame(rec).round(3); display(rec); rec.to_csv('results/t12_mapping_recall_patient_disjoint.csv', index=False)
dd = pd.DataFrame(dest).fillna(0); dd = dd.loc[dd.max(1).sort_values(ascending=False).index]; display((dd * 100).round(1))
dd.to_csv('results/t13_cycling_destinations.csv')
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6), gridspec_kw={'width_ratios': [1.35, 1]})
cols = {'T CD4 cycling': '#C4373A', 'T CD4 Helper': '#2A78D6', 'T reg': '#7FB2EA', 'T CD4 Naive': '#A9C8EE', 'T CD8 CTL': '#EB6834', 'T CD8 Naive': '#F3A983', 'NKT': '#1BAF7A', 'NK': '#7BD3B0'}
left = np.zeros(3)
for p in dd.index:
    v = dd.loc[p, ['R1', 'R2', 'R3']].values * 100
    ax[0].barh(['R1', 'R2', 'R3'], v, left=left, color=cols.get(p, '#C9CED8'), label=p, edgecolor='white'); left += v
ax[0].invert_yaxis(); ax[0].set_xlim(0, 100); ax[0].set_xlabel('% of true T CD4 cycling cells (PBISC read-out of real, patient-disjoint)')
ax[0].legend(frameon=False, fontsize=9, ncol=4, bbox_to_anchor=(0, -0.18), loc='upper left'); ax[0].set_title('Real T CD4 cycling cells: assigned label', fontsize=12)
x = np.arange(3); w = 0.38
for i, (lab, c) in enumerate([('T CD4 cycling', '#C4373A'), ('B_ASC_dividing', '#4A4F5C')]):
    r = rec[rec.label == lab]; ax[1].bar(x + (i - 0.5) * w, r.recall.fillna(0), w, color=c, label=f'{lab} (n = {r["true cells"].max()})')
    for xi, v in zip(x + (i - 0.5) * w, r.recall.fillna(0)): ax[1].text(xi, v + 0.02, f'{v:.2f}', ha='center', fontsize=9)
ax[1].set_xticks(x); ax[1].set_xticklabels(['R1', 'R2', 'R3']); ax[1].set_ylim(0, 1.05); ax[1].set_ylabel('recall'); ax[1].legend(frameon=False, fontsize=9)
ax[1].set_title('Recall: cycling CD4 T vs dividing plasmablast', fontsize=12)
plt.tight_layout(); plt.savefig('figs/f5_mapping_destinations.png', dpi=180); plt.show()""")
md("## F6. How many generated cells get the label")
code("""def mapping_counts(path):
    f = h5py.File(path, 'r'); o = f['obs']
    def col(g):
        if isinstance(g, h5py.Group) and 'codes' in g:
            cats = np.array([x.decode() for x in g['categories'][:]] + ['NA'], dtype=object); return cats[g['codes'][:]]
        return g[:]
    out = {}
    for r in ['R1', 'R3']:
        lab = col(o[f'{r}_fine']); out[r] = int((lab == 'T CD4 cycling').sum())
    if 'true_fine' in o: out['truth'] = int((col(o['true_fine']) == 'T CD4 cycling').sum())
    out['n'] = len(lab); return out
import h5py
base = '../proportion_qualification/results/'
M = {'Real, PBISC read-out (truth)': mapping_counts(base + 'covid_realT2_mapping.h5ad'), 'PBISC K = 1,000': mapping_counts(base + 'covid_k1000_mapping.h5ad'),
     'PBISC K = 4,000': mapping_counts(base + 'covid_k4000_mapping.h5ad')}
Mt = pd.DataFrame(M).T; display(Mt); Mt.to_csv('results/t14_label_counts_mapping.csv')
t5 = pd.read_csv('results/t5_cycling_gate_by_lineage.csv'); t10 = pd.read_csv('results/t10_format_aware_classifier.csv')
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6))
names = list(M); x = np.arange(3); w = 0.26
for i, (key, c) in enumerate([('truth', '#4A4F5C'), ('R1', '#2A78D6'), ('R3', '#EB6834')]):
    v = [M[n].get(key, np.nan) / M[n]['n'] * 100 for n in names]; cnt = [M[n].get(key, np.nan) for n in names]
    ax[0].bar(x + (i - 1) * w, np.nan_to_num(v), w, color=c, label={'truth': 'true label', 'R1': 'R1 call', 'R3': 'R3 call'}[key])
    for xi, vi, ci in zip(x + (i - 1) * w, v, cnt):
        if not np.isnan(vi): ax[0].text(xi, vi + 0.0006, f'{int(ci)}', ha='center', fontsize=9)
ax[0].set_xticks(x); ax[0].set_xticklabels(['Real\\n(PBISC read-out)', 'PBISC\\nK = 1,000', 'PBISC\\nK = 4,000']); ax[0].set_ylabel('% of all mapped cells')
ax[0].set_title('Cells labelled "T CD4 cycling" (numbers = cells; 37k real, 82k each generated)', fontsize=11); ax[0].legend(frameon=False, fontsize=9)
gate = [t5[(t5.lineage == 'T')]['real depth-matched %'].mean()] + list(t5[(t5.lineage == 'T')]['generated %'])
clf = list(t10['called % @0.9'])
for i, (vals, c, lab) in enumerate([(gate, '#C4373A', 'label-free: ≥2 of 14 cycle genes'), (clf, '#7A5CC2', 'format-aware classifier (P ≥ 0.9)')]):
    ax[1].bar(x + (i - 0.5) * 0.38, vals, 0.38, color=c, label=lab)
    for xi, v in zip(x + (i - 0.5) * 0.38, vals): ax[1].text(xi, v + 0.006, f'{v:.2f}%', ha='center', fontsize=9)
ax[1].set_xticks(x); ax[1].set_xticklabels(['Real\\n(depth-matched)', 'PBISC\\nK = 1,000', 'PBISC\\nK = 4,000']); ax[1].set_ylabel('% of T cells')
ax[1].set_title('Without the mapping: cycling T cells', fontsize=11); ax[1].legend(frameon=False, fontsize=9)
plt.tight_layout(); plt.savefig('figs/f6_generated_counts.png', dpi=180); plt.show()""")
md("## F7. Dot plot by label group: cycle genes, CD4 activation / Treg-Tfh genes, interferon genes")
code("""GENES = [('cell cycle', ['MKI67', 'TOP2A', 'TYMS', 'BIRC5', 'CENPF', 'STMN1']), ('CD4 / Treg-Tfh', ['CD4', 'CD40LG', 'FOXP3', 'CTLA4', 'ICOS', 'CXCL13', 'DUSP4', 'CD27']),
         ('CD8 / NK', ['CD8A', 'CD8B', 'NKG7']), ('interferon / activation', ['ISG15', 'IFI27', 'HLA-DRA', 'LGALS1'])]
genes = [g for _, gl in GENES for g in gl]
GROUPS = [('Real (held-out)', 'fine', 'T CD4 cycling', 'T CD4 cycling'), ('Real (held-out)', 'fine', 'T CD4 Helper', 'CD4 Helper'), ('Real (held-out)', 'fine', 'T reg', 'Treg'), ('Real (held-out)', 'fine', 'T CD8 CTL', 'CD8 CTL'),
          ('PBISC K = 1,000', 'pt_cyc', True, 'Prolif T, cycling'), ('PBISC K = 1,000', 'pt_cyc', False, 'Prolif T, not cycling'), ('PBISC K = 1,000', 'fine', 'Memory CD4 T cell', 'Memory CD4'), ('PBISC K = 1,000', 'fine', 'Treg', 'Treg'), ('PBISC K = 1,000', 'fine', 'Memory CD8 T cell', 'Memory CD8'),
          ('PBISC K = 4,000', 'pt_cyc', True, 'Prolif T, cycling'), ('PBISC K = 4,000', 'pt_cyc', False, 'Prolif T, not cycling'), ('PBISC K = 4,000', 'fine', 'Memory CD4 T cell', 'Memory CD4'), ('PBISC K = 4,000', 'fine', 'Treg', 'Treg'), ('PBISC K = 4,000', 'fine', 'Memory CD8 T cell', 'Memory CD8')]
REF = {'Real (held-out)': ('fine', 'T CD4 Helper'), 'PBISC K = 1,000': ('fine', 'Memory CD4 T cell'), 'PBISC K = 4,000': ('fine', 'Memory CD4 T cell')}
for k, A in D.items(): A.obs['pt_cyc'] = np.where(A.obs.fine == 'Proliferating T cell', A.obs.cycling, np.nan).astype(object)
def stats(A, m):
    R = A.raw[m.values if hasattr(m, 'values') else m][:, genes].X
    det = np.asarray((R > 0).mean(0)).ravel(); mean = np.asarray(np.expm1(R).mean(0)).ravel(); return det, mean
rows_ = []
for ds, key, val, name in GROUPS:
    A = D[ds]; m = (A.obs[key] == val).values
    det, mean = stats(A, m); rk, rv = REF[ds]; _, ref = stats(A, (A.obs[rk] == rv).values)
    rows_.append({'dataset': ds, 'group': name, 'n': int(m.sum()), 'det': det, 'l2r': np.log2((mean + 0.05) / (ref + 0.05))})
fig, ax = plt.subplots(figsize=(15, 8.2)); norm = mpl.colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3); cm = plt.cm.RdBu_r
for i, r in enumerate(rows_):
    ax.scatter(np.arange(len(genes)), np.full(len(genes), i), s=r['det'] * 420 + 2, c=cm(norm(np.clip(r['l2r'], -3, 3))), edgecolor='#8A909C', lw=0.3)
ax.set_yticks(range(len(rows_))); ax.set_yticklabels([f"{r['dataset'].replace('Real (held-out)', 'Real')} · {r['group']} (n = {r['n']:,})" for r in rows_], fontsize=10)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=55, ha='right', fontsize=10); ax.invert_yaxis()
x0 = 0
for name, gl in GENES:
    ax.axvline(x0 - 0.5, c='#D5D9E0', lw=1); ax.text(x0 + len(gl) / 2 - 0.5, -1.2, name, ha='center', fontsize=10, fontweight='bold'); x0 += len(gl)
for y in [3.5, 8.5]: ax.axhline(y, c='#4A4F5C', lw=0.8)
sm = plt.cm.ScalarMappable(norm=norm, cmap=cm); cb = plt.colorbar(sm, ax=ax, fraction=0.025, pad=0.01); cb.set_label('log2 mean vs CD4 reference (real: CD4 Helper; PBISC: Memory CD4)')
for s in [0.1, 0.5, 1.0]: ax.scatter([], [], s=s * 420 + 2, c='#C9CED8', edgecolor='#8A909C', label=f'{int(s*100)}% detected')
ax.legend(frameon=False, loc='lower left', bbox_to_anchor=(1.07, 0.0), fontsize=9); ax.set_ylim(len(rows_) - 0.5, -1.8)
plt.tight_layout(); plt.savefig('figs/f7_dotplot_label_groups.png', dpi=170); plt.show()
pd.DataFrame([{'dataset': r['dataset'], 'group': r['group'], 'n': r['n'], **{f'det_{g}': d for g, d in zip(genes, r['det'])}, **{f'l2r_{g}': l for g, l in zip(genes, r['l2r'])}} for r in rows_]).round(3).to_csv('results/t15_dotplot_values.csv', index=False)""")
md("## F8. Label-free sub-clusters (Leiden r = 1.0 on T cells of each display set): heatmap")
code("""HG = ['MKI67', 'TOP2A', 'TYMS', 'BIRC5', 'STMN1', 'CD3E', 'CD4', 'CD40LG', 'IL7R', 'CCR7', 'FOXP3', 'CTLA4', 'ICOS', 'CXCL13', 'DUSP4', 'CD27', 'CD8A', 'CD8B', 'GZMK', 'GZMB', 'NKG7', 'GNLY', 'ISG15', 'IFI27', 'HLA-DRA', 'LGALS1']
fig, axes = plt.subplots(1, 3, figsize=(20, 9.5), gridspec_kw={'wspace': 0.42}); top = []
for a, (k, A) in zip(axes, D.items()):
    T = A[~A.obs.label.isin(['NK', 'Proliferating NK'])]
    cl = T.obs.leiden.astype(str); g_ = [g for g in HG if g in T.raw.var_names]
    E = pd.DataFrame(T.raw[:, g_].X.toarray(), columns=g_, index=T.obs_names); E['cl'] = cl.values
    mu = E.groupby('cl').mean(); n = E.groupby('cl').size(); cyc = T.obs.groupby(cl.values).cycling.mean() * 100
    keep = n[n >= 40].index; mu, n, cyc = mu.loc[keep], n.loc[keep], cyc.loc[keep]
    z = (mu - mu.mean()) / (mu.std() + 1e-9); order = cyc.sort_values(ascending=False).index; z = z.loc[order]
    im = a.imshow(z.values, cmap='RdBu_r', vmin=-2.5, vmax=2.5, aspect='auto')
    a.set_xticks(range(len(g_))); a.set_xticklabels(g_, rotation=70, ha='right', fontsize=9)
    a.set_yticks(range(len(order))); a.set_yticklabels([f'c{c} · n {n[c]:,} · cyc {cyc[c]:.0f}%' for c in order], fontsize=9)
    a.set_title(f'{k}\\nT-cell sub-clusters (sorted by % cycling)', fontsize=12)
    a.add_patch(plt.Rectangle((-0.5, -0.5), len(g_), 1, fill=False, ec='#C4373A', lw=2.5))
    c0 = order[0]; top.append({'dataset': k, 'cluster': c0, 'n': int(n[c0]), 'pct_cycling': round(cyc[c0], 1), **{f'mean_{g}': round(mu.loc[c0, g], 3) for g in g_}})
cb = fig.colorbar(im, ax=axes, fraction=0.012, pad=0.01); cb.set_label('cluster mean, z-score across clusters')
plt.savefig('figs/f8_subcluster_heatmap.png', dpi=160, bbox_inches='tight'); plt.show()
TOP = pd.DataFrame(top); display(TOP[['dataset', 'cluster', 'n', 'pct_cycling'] + [c for c in TOP.columns if c.startswith('mean_') and c[5:] in ['MKI67','CD4','CD40LG','FOXP3','CTLA4','CXCL13','CD8A','IFI27']]]); TOP.to_csv('results/t16_top_cycling_cluster.csv', index=False)""")
code("""# composition of the top cycling cluster in each dataset (labels, CD4/CD8 markers)
for k, A in D.items():
    T = A[~A.obs.label.isin(['NK', 'Proliferating NK'])]; c0 = TOP[TOP.dataset == k].cluster.iloc[0]; m = (T.obs.leiden.astype(str) == c0).values
    R = T.raw[m][:, ['CD4', 'CD40LG', 'FOXP3', 'CD8A', 'CD8B']].X.toarray() > 0
    cd4 = (R[:, :3].any(1)) & ~(R[:, 3:].any(1)); cd8 = (R[:, 3:].any(1)) & ~(R[:, :3].any(1))
    print(k, 'cluster', c0, 'n', m.sum(), '| labels', dict(collections.Counter(T.obs.fine[m]).most_common(4)), '| CD4-marked', round(cd4.mean() * 100, 1), '% CD8-marked', round(cd8.mean() * 100, 1), '%')""")
md("## F9. Panel coverage of the real T CD4 cycling markers")
code("""_, _, varg = rows('gen_k1000'); PANEL = set(varg)
de = pd.read_csv('results/de_recap_cd4cycling_vs_cd4.csv')
up = de[(de.logfoldchanges > 1) & (de.pvals_adj < 0.05) & (de.det_cyc >= 0.25) & (de.det_cyc - de.det_cd4 >= 0.2)].sort_values('scores', ascending=False).reset_index(drop=True)
up['in_panel'] = up.names.isin(PANEL)
cyc_set = set(GS['S']) | set(GS['G2M']) | set(GS['EXTRA']); upc = up[up.names.isin(cyc_set)]
cov = pd.DataFrame({'set': ['top 10', 'top 25', 'top 50', f'all ({len(up)})', f'cell-cycle markers ({len(upc)})', 'background'],
                    'in panel': [up.head(10).in_panel.sum(), up.head(25).in_panel.sum(), up.head(50).in_panel.sum(), up.in_panel.sum(), upc.in_panel.sum(), np.nan],
                    'fraction': [up.head(10).in_panel.mean(), up.head(25).in_panel.mean(), up.head(50).in_panel.mean(), up.in_panel.mean(), upc.in_panel.mean(), 0.293]})
display(cov.round(3)); cov.to_csv('results/t17_panel_coverage.csv', index=False)
fig, ax = plt.subplots(1, 2, figsize=(16, 6.4), gridspec_kw={'width_ratios': [1, 1.7]})
ax[0].barh(cov.set, cov.fraction * 100, color=['#2A78D6'] * 4 + ['#C4373A', '#9AA3B2']); ax[0].invert_yaxis()
for i, (f, c) in enumerate(zip(cov.fraction, cov['in panel'])): ax[0].text(f * 100 + 1, i, f'{f*100:.0f}%' + ('' if np.isnan(c) else f' ({int(c)})'), va='center', fontsize=10)
ax[0].set_xlim(0, 100); ax[0].set_xlabel('% of markers in the 6k panel'); ax[0].set_title('Real T CD4 cycling markers in the PBISC panel', fontsize=12)
t40 = up.head(40)
ax[1].bar(range(40), t40.logfoldchanges, color=np.where(t40.in_panel, '#2A78D6', '#C4373A'))
ax[1].set_xticks(range(40)); ax[1].set_xticklabels(t40.names, rotation=75, ha='right', fontsize=8.5); ax[1].set_ylabel('log2FC vs other CD4 T (real, full transcriptome)')
from matplotlib.patches import Patch
ax[1].legend(handles=[Patch(color='#2A78D6', label='in panel'), Patch(color='#C4373A', label='not in panel')], frameon=False); ax[1].set_title('Top 40 markers (recap label, 192 cells)', fontsize=12)
plt.tight_layout(); plt.savefig('figs/f9_panel_coverage.png', dpi=170); plt.show()
print(pd.read_csv('results/t3_panel_only_classifier.csv'))""")
nb = nbf.v4.new_notebook(); nb['cells'] = C; nb['metadata'] = {'kernelspec': {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}}
nbf.write(nb, '/home/claude/cd4cyc/cd4_cycling_panel_figs.ipynb'); print('ok')
