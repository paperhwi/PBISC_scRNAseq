import nbformat as nbf, sys
from nbclient import NotebookClient
OUT = "/mnt/user-data/outputs/integration_report_20261012/nk_ilc_markers"
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# NK / ILC marker genes in the PBISC 6k panel and in generated cells

Question: among the ~6,000 genes PBISC generates, are the genes known to be NK- or ILC-specific present, and do generated cells express them in the right cells?

1. Panel membership of a curated NK / ILC marker list (sepsis panel 5,975 genes, COVID panel 5,927 genes).
2. Detection rate (fraction of cells with count > 0) and mean log1p CP10k by cell group, real vs generated (K = 1,000 and K = 4,000).
   Real = recapitulated references (Ye 2026 sepsis in-scheme cells; Stephenson 2025 COVID, 187-sample roster), up to 3,000 cells per group, all ILC cells.
   Generated = the 100-cells-per-sample-per-draw joint subsets; groups from the format-aware classifier (route R3, fine level).
   Both normalised on the shared panel genes only.
3. Specificity call per marker and source.
4. Label-free checks: co-detection of NK markers with CD3 genes, and an ILC-like gate.""")

code("""import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, sys, os, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/work'); import store
U = '/mnt/user-data/uploads/PBISC/.claude_pylib/extract'; P = '/home/claude/work/prop'; R = '/home/claude/work/report'
OUT = '%s'; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(11)
pd.set_option('display.width', 220); pd.set_option('display.max_rows', 200)""" % OUT)

code("""# curated markers: gene, target, note (protein name / known role)
M = [
 # NK core identity (receptors largely absent from T cells)
 ('NCAM1','NK','CD56'), ('FCGR3A','NK','CD16; also CD16 monocytes'), ('KLRF1','NK','NKp80'), ('NCR1','NK','NKp46'),
 ('NCR3','NK','NKp30'), ('SH2D1B','NK','EAT-2'), ('KLRC1','NK','NKG2A'), ('KLRD1','NK','CD94; also CD8/NKT'),
 ('KLRK1','NK','NKG2D; also CD8'), ('KIR2DL1','NK','KIR'), ('KIR2DL3','NK','KIR'), ('KIR2DL4','NK','KIR'),
 ('KIR3DL1','NK','KIR'), ('KIR3DL2','NK','KIR'), ('KIR2DS4','NK','KIR'), ('CD160','NK','also CD8'),
 ('MYOM2','NK','PBMC scRNA NK marker'), ('AKR1C3','NK','PBMC scRNA NK marker'), ('SPON2','NK','PBMC scRNA NK marker'),
 ('FGFBP2','NK','also CTL'), ('B3GAT1','NK','CD57; also CTL'), ('TYROBP','NK','DAP12; also myeloid'),
 ('FCER1G','NK','FcRg; also myeloid'), ('IL2RB','NK','CD122'), ('TBX21','NK','T-bet'), ('EOMES','NK','also CD8'),
 # cytotoxic programme (NK + CTL)
 ('GNLY','NK','granulysin'), ('NKG7','NK','NKG7'), ('PRF1','NK','perforin'), ('GZMB','NK','granzyme B'),
 ('GZMH','NK','granzyme H'), ('GZMA','NK','granzyme A'),
 # CD56bright NK
 ('XCL1','NK','CD56bright'), ('XCL2','NK','CD56bright'), ('GZMK','NK','CD56bright; also CD8'), ('SELL','NK','CD62L, CD56bright'),
 # ILC (helper-like innate lymphoid)
 ('IL7R','ILC','CD127; also T'), ('KIT','ILC','CD117, ILCP / ILC3'), ('PTGDR2','ILC','CRTH2, ILC2'), ('HPGDS','ILC','ILC2'),
 ('IL1RL1','ILC','ST2, ILC2'), ('IL17RB','ILC','IL-25R, ILC2'), ('GATA3','ILC','ILC2; also Th2'), ('KLRG1','ILC','ILC2; also NK/CTL'),
 ('RORC','ILC','RORgt, ILC3'), ('AHR','ILC','ILC3'), ('IL23R','ILC','ILC3'), ('IL22','ILC','ILC3'), ('CD200R1','ILC','ILCP'),
 ('ID2','ILC','ILC/NK TF'), ('KLRB1','ILC','CD161; also MAIT/NK'), ('IL4','ILC','ILC2 cytokine'), ('IL5','ILC','ILC2 cytokine'),
 ('IL13','ILC','ILC2 cytokine'), ('AREG','ILC','ILC2'), ('ICOS','ILC','also T'), ('TNFRSF4','ILC','OX40; also T'),
 # lineage exclusion markers
 ('CD3E','exclude','T'), ('CD3D','exclude','T'), ('CD3G','exclude','T'), ('TRDC','exclude','gdT; NK transcripts'),
 ('CD8A','exclude','CD8 T'), ('CD4','exclude','CD4 T / mono'), ('MS4A1','exclude','B'), ('CD79A','exclude','B'),
 ('CD14','exclude','monocyte'), ('LYZ','exclude','myeloid')]
M = pd.DataFrame(M, columns=['gene', 'target', 'note'])
gen_var = {c: np.load(f'{U}/{c}_k1000_counts.npz', allow_pickle=True)['var'].astype(str) for c in ['sepsis', 'covid']}
M['in_sepsis_panel'] = M.gene.isin(gen_var['sepsis']); M['in_covid_panel'] = M.gene.isin(gen_var['covid'])
print({c: len(v) for c, v in gen_var.items()})
M""")

code("""s = M.groupby('target')[['in_sepsis_panel', 'in_covid_panel']].agg(['sum', 'count']); print(s)
print('NOT in sepsis panel:', M.loc[~M.in_sepsis_panel, 'gene'].tolist())
print('NOT in COVID panel :', M.loc[~M.in_covid_panel, 'gene'].tolist())
M.to_csv(f'{OUT}/marker_panel_membership.csv', index=False)""")

code("""GROUPS = ['NK', 'ILC', 'NKT', 'T', 'B/Plasma', 'Myeloid', 'Other']
def grp(x):
    x = str(x)
    if x == 'ILC': return 'ILC'
    if x.startswith('NKT'): return 'NKT'
    if x.startswith('NK'): return 'NK'
    if x in ('Unresolved (low-quality T/NK)', 'Platelet-T/NK complex', 'nan', ''): return None
    if x.startswith(('T_', 'T ', 'gdT', 'MAIT')): return 'T'
    if x.startswith('B_'): return 'B/Plasma'
    if x.startswith(('Mono', 'CD14mono', 'CD16mono', 'Int.mono', 'DC', 'cDC', 'pDC', 'ASDC', 'Macrophage')): return 'Myeloid'
    return 'Other'

def lognorm(X):
    X = X.tocsr().astype(np.float32); tot = np.asarray(X.sum(1)).ravel()
    Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y

def load_real(c, genes, cap=3000):
    if c == 'sepsis':
        st = store.Store('/home/claude/work/sepsis/store'); sym = st.var.gene_symbol.astype(str).to_numpy()
        ro = pd.concat([pd.read_csv(f'/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz', keep_default_na=False,
                        usecols=['celltype_ye2026', 'in_paper_scheme']) for i in range(3)], ignore_index=True)
        lab = ro.celltype_ye2026.astype(str).to_numpy(); ok = ro.in_paper_scheme.astype(str).eq('True').to_numpy()
    else:
        st = store.Store('/home/claude/work/covid/store'); sym = st.var.gene.astype(str).to_numpy()
        L = pd.read_csv('/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz', index_col=0, keep_default_na=False)
        lab = L.cell_type_detailed_recap.astype(str).to_numpy()
        roster = set(pd.read_csv(f'{R}/covid_roster187_local_ids.csv').sample_id.astype(str))
        ok = st.obs.sample_id.astype(str).isin(roster).to_numpy()
    g = np.array([grp(x) for x in lab], dtype=object)
    idx = []
    for G in GROUPS:
        w = np.where(ok & (g == G))[0]; idx.append(rng.permutation(w)[:cap])
    idx = np.sort(np.concatenate(idx))
    pos = pd.Series(np.arange(len(sym)), index=sym); pos = pos[~pos.index.duplicated()]
    keep = list(dict.fromkeys(x for x in genes if x in pos.index))
    X = st.rows(idx)[:, pos.loc[keep].to_numpy()]
    return X, keep, g[idx], lab[idx]

def load_gen(c, K):
    z = np.load(f'{U}/{c}_{K}_counts.npz', allow_pickle=True)
    X = sp.csr_matrix((z['data'].astype(np.float32), z['indices'].astype(np.int32), z['indptr']), shape=(len(z['rows']), len(z['var'])))
    lab = ad.read_h5ad(f'{P}/{c}_{K}_mapping.h5ad', backed='r').obs['R3_fine'].astype(str).to_numpy()
    return X, z['var'].astype(str), np.array([grp(x) for x in lab], dtype=object), lab""")

code("""res = []; data = {}
for c in ['sepsis', 'covid']:
    gens = {K: load_gen(c, K) for K in ['k1000', 'k4000']}
    gvar = gens['k1000'][1]
    Xr, shared, gr, lr = load_real(c, list(gvar))                      # real restricted to panel genes it shares
    Yr = lognorm(Xr); posr = pd.Series(np.arange(len(shared)), index=shared)
    srcs = {'real': (Xr, Yr, posr, gr)}
    for K, (Xg, gv, gg, lg) in gens.items():
        gpos = pd.Series(np.arange(len(gv)), index=gv); gpos = gpos[~gpos.index.duplicated()]   # 6 duplicated symbols in the panel: keep first
        cols = gpos.loc[shared].to_numpy(); Xs = Xg[:, cols].tocsr()
        srcs[K] = (Xs, lognorm(Xs), posr, gg)
    data[c] = srcs
    for src, (X, Y, pos, g) in srcs.items():
        for G in GROUPS:
            m = g == G; n = int(m.sum())
            for gene in M.gene:
                if gene not in pos.index: continue
                j = pos[gene]
                det = float((X[m][:, j].toarray().ravel() > 0).mean()) if n else np.nan
                mu = float(Y[m][:, j].toarray().ravel().mean()) if n else np.nan
                res.append(dict(cohort=c, source=src, group=G, n_cells=n, gene=gene, det=det, mean_lognorm=mu))
    print(c, {s: pd.Series(v[3]).value_counts().reindex(GROUPS).fillna(0).astype(int).to_dict() for s, v in srcs.items()})
res = pd.DataFrame(res); res.to_csv(f'{OUT}/marker_detection_by_group.csv.gz', index=False)""")

code("""# specificity call: target group detection >= 0.10 and >= 2x the highest detection among clearly different lineages
OTHERS = {'NK': ['T', 'B/Plasma', 'Myeloid', 'Other'], 'ILC': ['NK', 'T', 'B/Plasma', 'Myeloid', 'Other']}
rows = []
for (c, src, gene), d in res.groupby(['cohort', 'source', 'gene']):
    t = M.set_index('gene').loc[gene, 'target']
    if t not in OTHERS: continue
    d = d.set_index('group')
    if d.loc[t, 'n_cells'] < 20: tdet = np.nan
    else: tdet = d.loc[t, 'det']
    oth = d.loc[OTHERS[t], 'det']; best = oth.idxmax()
    rows.append(dict(cohort=c, source=src, gene=gene, target=t, det_target=tdet, best_other=best, det_best_other=oth.max(),
                     ratio=(tdet + 0.005) / (oth.max() + 0.005) if tdet == tdet else np.nan))
S = pd.DataFrame(rows); S['specific'] = (S.det_target >= 0.10) & (S.ratio >= 2)
W = S.pivot_table(index=['target', 'gene'], columns=['cohort', 'source'], values=['det_target', 'ratio', 'specific'], aggfunc='first')
S.to_csv(f'{OUT}/marker_specificity.csv', index=False)
pd.set_option('display.float_format', '{:.2f}'.format)
for c in ['sepsis', 'covid']:
    T = S[S.cohort == c].pivot_table(index=['target', 'gene'], columns='source', values=['det_target', 'ratio'], aggfunc='first')
    T.columns = [f'{a}_{b}' for a, b in T.columns]
    sp_ = S[S.cohort == c].pivot_table(index=['target', 'gene'], columns='source', values='specific', aggfunc='first')
    T['specific_real'] = sp_['real']; T['specific_k1000'] = sp_['k1000']; T['specific_k4000'] = sp_['k4000']
    print('=====', c); display(T[['det_target_real', 'det_target_k1000', 'det_target_k4000', 'ratio_real', 'ratio_k1000', 'ratio_k4000', 'specific_real', 'specific_k1000', 'specific_k4000']])""")

code("""# heatmap: detection by group, real vs K = 1,000 (K = 4,000 in the csv)
show = M[M.target.isin(['NK', 'ILC'])].gene.tolist()
fig, axes = plt.subplots(1, 4, figsize=(17, 15), sharey=True)
for k, (c, src) in enumerate([('sepsis', 'real'), ('sepsis', 'k1000'), ('covid', 'real'), ('covid', 'k1000')]):
    d = res[(res.cohort == c) & (res.source == src)].pivot(index='gene', columns='group', values='det').reindex(index=show, columns=GROUPS)
    ax = axes[k]; im = ax.imshow(d.to_numpy(), aspect='auto', cmap='Blues', vmin=0, vmax=1)
    ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels(GROUPS, rotation=60, ha='right', fontsize=10)
    ax.set_yticks(range(len(show))); ax.set_yticklabels(show, fontsize=9)
    n = res[(res.cohort == c) & (res.source == src)].drop_duplicates('group').set_index('group').n_cells.reindex(GROUPS)
    ax.set_title(f'{c} · {\"real\" if src == \"real\" else \"generated K=1,000\"}\\n' + ' '.join(f'{g}:{int(v)}' for g, v in n.items() if g in ('NK', 'ILC')), fontsize=11)
    for i in range(d.shape[0]):
        if d.iloc[i].isna().all(): ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(GROUPS), 1, color='#dddddd'))
    ax.axhline(M[M.target == 'NK'].shape[0] - 0.5, color='k', lw=1)
fig.colorbar(im, ax=axes, shrink=0.4, label='fraction of cells detected (count > 0)')
fig.savefig(f'{OUT}/nk_ilc_marker_detection_heatmap.png', dpi=130, bbox_inches='tight'); plt.show()""")

code("""# label-free check 1: cells detecting an NK-specific receptor gene, and what fraction also detect a CD3 gene
def frac(X, pos, genes):
    genes = [g for g in genes if g in pos.index]
    return (X[:, pos.loc[genes].to_numpy()].sum(1).A1 > 0) if genes else None
out = []
for c in data:
    for src, (X, Y, pos, g) in data[c].items():
        nkr = frac(X, pos, ['KLRF1', 'NCR1', 'SH2D1B']); cd3 = frac(X, pos, ['CD3E', 'CD3D', 'CD3G'])
        out.append(dict(cohort=c, source=src, cells=X.shape[0], nk_receptor_pos=nkr.mean(), cd3_pos_among_nkr=cd3[nkr].mean(),
                        nk_receptor_pos_in_T=nkr[g == 'T'].mean(), nk_receptor_pos_in_NK=nkr[g == 'NK'].mean()))
L1 = pd.DataFrame(out); L1""")

code("""# label-free check 2: ILC-like gate = IL7R+ and (KIT or PTGDR2 or IL1RL1 or RORC)+ ; CD3D/E/G-, MS4A1/CD79A-, CD14/LYZ-
out = []
for c in data:
    for src, (X, Y, pos, g) in data[c].items():
        il7r = frac(X, pos, ['IL7R']); ilc = frac(X, pos, ['KIT', 'PTGDR2', 'IL1RL1', 'RORC'])
        neg = ~frac(X, pos, ['CD3E', 'CD3D', 'CD3G']) & ~frac(X, pos, ['MS4A1', 'CD79A']) & ~frac(X, pos, ['CD14', 'LYZ'])
        gate = il7r & ilc & neg
        out.append(dict(cohort=c, source=src, cells=X.shape[0], ilc_like=int(gate.sum()),
                        ilc_like_in_labelled_ILC=float(gate[g == 'ILC'].mean()) if (g == 'ILC').sum() else np.nan,
                        labelled_ILC=int((g == 'ILC').sum()),
                        IL7R_pos=il7r.mean(), KIT_or_PTGDR2_or_IL1RL1_or_RORC_pos=ilc.mean()))
L2 = pd.DataFrame(out)
print('Note: real cells were sampled per group (ILC enriched), so ilc_like counts are not population frequencies; compare ilc_like_in_labelled_ILC and marker positivity.')
L2""")

code("""L1.to_csv(f'{OUT}/labelfree_nk_receptor_cd3.csv', index=False); L2.to_csv(f'{OUT}/labelfree_ilc_gate.csv', index=False)
# h5ad: per-cell marker counts for the NK/ILC genes, all sources, for inspection
mk = M.gene.tolist(); obj = []
for c in data:
    for src, (X, Y, pos, g) in data[c].items():
        gg = [x for x in mk if x in pos.index]
        a = ad.AnnData(sp.csr_matrix(X[:, pos.loc[gg].to_numpy()]), obs=pd.DataFrame({'cohort': c, 'source': src, 'group': pd.Series(g).astype(str).values}),
                       var=pd.DataFrame(index=gg))
        a.obs_names = [f'{c}_{src}_{i}' for i in range(a.n_obs)]; obj.append(a)
A = ad.concat(obj, join='outer', fill_value=0); A.obs = A.obs.astype('category'); A.write_h5ad(f'{OUT}/nk_ilc_marker_counts.h5ad')
A""")

nb = nbf.v4.new_notebook(); nb.cells = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
NotebookClient(nb, timeout=3600, kernel_name='python3').execute()
nbf.write(nb, f'{OUT}/nk_ilc_markers.ipynb'); print('done')
