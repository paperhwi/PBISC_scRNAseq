import nbformat as nbf
from nbclient import NotebookClient
OUT = "/mnt/user-data/outputs/integration_report_20261012/macrophage_markers"
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# Macrophage marker genes in the PBISC 6k panel and in generated cells

Question: among the ~6,000 genes PBISC generates, are genes that define a macrophage profile present, and do generated cells express them in a macrophage-like pattern?

1. Panel membership of a curated macrophage list (complement C1Q, scavenger receptors, lipid / APOE, tissue-resident, lysosomal and M2-like genes), with monocyte contrast genes.
2. Detection rate by cell group, real vs generated (K = 1,000 and K = 4,000), normalised on shared panel genes.
   - Real: recapitulated references. Ye 2026 sepsis has a 'Macrophage' label (3,974 cells); the Stephenson 2025 COVID scheme has none.
   - Generated: groups from the format-aware classifier (route R3, fine). Route-3 recovery of real macrophages in PBISC format: recall 0.51, precision 0.56 (sepsis major level).
3. Specificity call: macrophage group vs CD14 monocytes and vs the best other group.
4. Label-free: co-detection of C1Q genes and a macrophage gate (C1Q + APOE/CD163/MRC1/MARCO/SPP1, no lymphoid genes) across all cells.""")

code("""import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/claude/work'); import store
U = '/mnt/user-data/uploads/PBISC/.claude_pylib/extract'; P = '/home/claude/work/prop'; R = '/home/claude/work/report'
OUT = '%s'; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(13)
pd.set_option('display.width', 230); pd.set_option('display.max_rows', 200); pd.set_option('display.float_format', '{:.2f}'.format)""" % OUT)

code("""M = [
 # complement C1Q (tissue / blood macrophages; also some CD16 monocytes)
 ('C1QA', 'macrophage', 'complement C1q'), ('C1QB', 'macrophage', 'complement C1q'), ('C1QC', 'macrophage', 'complement C1q'),
 # scavenger / pattern receptors
 ('CD163', 'macrophage', 'haemoglobin scavenger; also some monocytes'), ('MRC1', 'macrophage', 'CD206, M2-like'), ('MSR1', 'macrophage', 'CD204'),
 ('MARCO', 'macrophage', 'scavenger receptor'), ('SIGLEC1', 'macrophage', 'CD169; IFN-induced in monocytes'), ('CD209', 'macrophage', 'DC-SIGN'),
 ('STAB1', 'macrophage', 'stabilin-1'), ('LYVE1', 'macrophage', 'tissue-resident'), ('FOLR2', 'macrophage', 'tissue-resident'),
 ('TIMD4', 'macrophage', 'tissue-resident'), ('VSIG4', 'macrophage', 'complement receptor'), ('CD5L', 'macrophage', 'splenic/BM macrophage'),
 ('MERTK', 'macrophage', 'efferocytosis'), ('AXL', 'macrophage', 'efferocytosis'), ('MS4A4A', 'macrophage', 'M2-like'),
 # lipid / lysosomal / tissue macrophage programme
 ('APOE', 'macrophage', 'lipid'), ('APOC1', 'macrophage', 'lipid'), ('TREM2', 'macrophage', 'lipid-associated'), ('SPP1', 'macrophage', 'osteopontin'),
 ('GPNMB', 'macrophage', 'lysosomal'), ('LGMN', 'macrophage', 'legumain'), ('CTSK', 'macrophage', 'cathepsin K'), ('ACP5', 'macrophage', 'TRAP'),
 ('LIPA', 'macrophage', 'lysosomal lipase'), ('FABP5', 'macrophage', 'lipid'), ('CHIT1', 'macrophage', 'chitotriosidase'), ('CHI3L1', 'macrophage', 'YKL-40; also neutrophil'),
 ('PLA2G7', 'macrophage', 'Lp-PLA2'), ('CCL18', 'macrophage', 'M2 chemokine'), ('F13A1', 'macrophage', 'M2-like'), ('SELENOP', 'macrophage', 'selenoprotein P'),
 ('SLC40A1', 'macrophage', 'ferroportin, iron recycling'), ('HMOX1', 'macrophage', 'haem oxygenase'), ('NUPR1', 'macrophage', 'stress / tissue'),
 # transcription factors
 ('MAF', 'macrophage', 'c-Maf'), ('MAFB', 'macrophage', 'MafB; also monocytes'),
 # shared myeloid, for context
 ('CD68', 'myeloid', 'also monocytes / DC'), ('CSF1R', 'myeloid', 'also monocytes'), ('FCGR1A', 'myeloid', 'CD64'), ('ITGAM', 'myeloid', 'CD11b'),
 # monocyte contrast genes (should be lower in macrophages)
 ('CD14', 'monocyte', 'classical monocyte'), ('LYZ', 'monocyte', ''), ('S100A8', 'monocyte', ''), ('S100A9', 'monocyte', ''), ('VCAN', 'monocyte', ''),
 ('FCN1', 'monocyte', ''), ('SELL', 'monocyte', ''), ('CCR2', 'monocyte', '')]
M = pd.DataFrame(M, columns=['gene', 'target', 'note'])
gen_var = {c: np.load(f'{U}/{c}_k1000_counts.npz', allow_pickle=True)['var'].astype(str) for c in ['sepsis', 'covid']}
M['in_sepsis_panel'] = M.gene.isin(gen_var['sepsis']); M['in_covid_panel'] = M.gene.isin(gen_var['covid'])
print(M.groupby('target')[['in_sepsis_panel', 'in_covid_panel']].agg(['sum', 'count']))
print('NOT in sepsis panel:', M.loc[~M.in_sepsis_panel, 'gene'].tolist()); print('NOT in COVID panel :', M.loc[~M.in_covid_panel, 'gene'].tolist())
M.to_csv(f'{OUT}/macrophage_marker_panel_membership.csv', index=False); M""")

code("""GROUPS = ['Macrophage', 'CD14 mono', 'CD16 mono', 'DC', 'T/NK', 'B/Plasma', 'Other']
def grp(x):
    x = str(x)
    if x in ('Unresolved (low-quality T/NK)', 'Platelet-T/NK complex', 'nan', '', 'NA'): return None
    if x == 'Macrophage': return 'Macrophage'
    if x.startswith(('Mono_c07', 'CD16mono')) or x.startswith('Mono_c08') or x == 'Int.mono': return 'CD16 mono'
    if x.startswith(('Mono_c0', 'CD14mono')): return 'CD14 mono'
    if x.startswith(('cDC', 'pDC', 'DC', 'ASDC')): return 'DC'
    if x.startswith(('T_', 'T ', 'gdT', 'MAIT', 'NK', 'ILC')): return 'T/NK'
    if x.startswith('B_'): return 'B/Plasma'
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
        roster = set(pd.read_csv(f'{R}/covid_roster187_local_ids.csv').sample_id.astype(str)); ok = st.obs.sample_id.astype(str).isin(roster).to_numpy()
    g = np.array([grp(x) for x in lab], dtype=object)
    idx = np.sort(np.concatenate([rng.permutation(np.where(ok & (g == G))[0])[:cap] for G in GROUPS]))
    pos = pd.Series(np.arange(len(sym)), index=sym); pos = pos[~pos.index.duplicated()]
    keep = list(dict.fromkeys(x for x in genes if x in pos.index))
    return st.rows(idx)[:, pos.loc[keep].to_numpy()], keep, g[idx]
def load_gen(c, K):
    z = np.load(f'{U}/{c}_{K}_counts.npz', allow_pickle=True)
    X = sp.csr_matrix((z['data'].astype(np.float32), z['indices'].astype(np.int32), z['indptr']), shape=(len(z['rows']), len(z['var'])))
    lab = ad.read_h5ad(f'{P}/{c}_{K}_mapping.h5ad', backed='r').obs['R3_fine'].astype(str).to_numpy()
    return X, z['var'].astype(str), np.array([grp(x) for x in lab], dtype=object)""")

code("""res = []; data = {}
for c in ['sepsis', 'covid']:
    gens = {K: load_gen(c, K) for K in ['k1000', 'k4000']}; gvar = gens['k1000'][1]
    Xr, shared, gr = load_real(c, list(gvar)); posr = pd.Series(np.arange(len(shared)), index=shared)
    srcs = {'real': (Xr, gr)}
    for K, (Xg, gv, gg) in gens.items():
        gpos = pd.Series(np.arange(len(gv)), index=gv); gpos = gpos[~gpos.index.duplicated()]   # 6 duplicated symbols in the panel
        srcs[K] = (Xg[:, gpos.loc[shared].to_numpy()].tocsr(), gg)
    data[c] = (srcs, posr)
    for src, (X, g) in srcs.items():
        Y = lognorm(X)
        for G in GROUPS:
            m = g == G; n = int(m.sum())
            for gene in M.gene:
                if gene not in posr.index: continue
                j = posr[gene]
                res.append(dict(cohort=c, source=src, group=G, n_cells=n, gene=gene,
                                det=float((X[m][:, j].toarray().ravel() > 0).mean()) if n else np.nan,
                                mean_lognorm=float(Y[m][:, j].toarray().ravel().mean()) if n else np.nan))
    print(c, {s: pd.Series(v[1]).value_counts().reindex(GROUPS).fillna(0).astype(int).to_dict() for s, v in srcs.items()})
res = pd.DataFrame(res); res.to_csv(f'{OUT}/macrophage_marker_detection_by_group.csv.gz', index=False)""")

md("""## Specificity (sepsis: real macrophage label vs generated cells labelled macrophage)

`ratio_vs_mono` = detection in macrophages / detection in CD14 monocytes; `ratio_vs_best` = vs the highest of all other groups (CD16 monocytes included). Called specific when macrophage detection ≥ 0.10 and ratio_vs_best ≥ 2.""")

code("""rows = []
for (c, src, gene), d in res.groupby(['cohort', 'source', 'gene']):
    d = d.set_index('group')
    if d.loc['Macrophage', 'n_cells'] < 20: continue
    t = d.loc['Macrophage', 'det']; oth = d.drop('Macrophage').det
    rows.append(dict(cohort=c, source=src, gene=gene, target=M.set_index('gene').loc[gene, 'target'], det_macro=t, det_CD14mono=d.loc['CD14 mono', 'det'],
                     det_CD16mono=d.loc['CD16 mono', 'det'], best_other=oth.idxmax(), det_best_other=oth.max(),
                     ratio_vs_mono=(t + .005) / (d.loc['CD14 mono', 'det'] + .005), ratio_vs_best=(t + .005) / (oth.max() + .005)))
S = pd.DataFrame(rows); S['specific'] = (S.det_macro >= .10) & (S.ratio_vs_best >= 2)
S.to_csv(f'{OUT}/macrophage_marker_specificity.csv', index=False)
T = S[S.cohort == 'sepsis'].pivot_table(index=['target', 'gene'], columns='source', values=['det_macro', 'det_CD14mono', 'ratio_vs_best', 'specific'], aggfunc='first')
T.columns = [f'{a}|{b}' for a, b in T.columns]
T[['det_macro|real', 'det_macro|k1000', 'det_macro|k4000', 'det_CD14mono|real', 'det_CD14mono|k1000', 'ratio_vs_best|real', 'ratio_vs_best|k1000', 'ratio_vs_best|k4000', 'specific|real', 'specific|k1000', 'specific|k4000']]""")

code("""# where do macrophage genes appear in COVID generated cells (no macrophage label in the COVID scheme): detection by group
x = res[(res.cohort == 'covid') & res.gene.isin(M[M.target == 'macrophage'].gene)].pivot_table(index='gene', columns=['source', 'group'], values='det')
x[[('real', 'CD14 mono'), ('real', 'CD16 mono'), ('real', 'DC'), ('k1000', 'CD14 mono'), ('k1000', 'CD16 mono'), ('k1000', 'DC')]]""")

md("""## Label-free checks

- C1Q co-detection: cells with ≥ 2 of C1QA/C1QB/C1QC (in real blood this is a macrophage / C1Q+ CD16-monocyte trait).
- Macrophage gate: ≥ 2 C1Q genes **and** ≥ 1 of APOE, CD163, MRC1, MARCO, SPP1, SELENOP, **and** no CD3E/CD3D/MS4A1/CD79A counts.
Real cells were sampled per group (macrophages enriched), so compare rates within groups, not totals.""")

code("""def anyg(X, pos, genes, k=1):
    gi = [pos[g] for g in genes if g in pos.index]
    return ((X[:, gi] > 0).sum(1).A1 >= k) if gi else np.zeros(X.shape[0], bool)
rows = []
for c, (srcs, pos) in data.items():
    for src, (X, g) in srcs.items():
        c1q = anyg(X, pos, ['C1QA', 'C1QB', 'C1QC'], 2)
        gate = c1q & anyg(X, pos, ['APOE', 'CD163', 'MRC1', 'MARCO', 'SPP1', 'SELENOP']) & ~anyg(X, pos, ['CD3E', 'CD3D', 'MS4A1', 'CD79A'])
        for G in GROUPS:
            m = g == G
            if m.sum() < 20: continue
            rows.append(dict(cohort=c, source=src, group=G, cells=int(m.sum()), C1Q_2plus=c1q[m].mean(), macro_gate=gate[m].mean()))
        rows.append(dict(cohort=c, source=src, group='ALL cells', cells=X.shape[0], C1Q_2plus=c1q.mean(), macro_gate=gate.mean()))
LF = pd.DataFrame(rows); LF.to_csv(f'{OUT}/macrophage_labelfree.csv', index=False)
LF.pivot_table(index=['cohort', 'group'], columns='source', values=['C1Q_2plus', 'macro_gate'])""")

code("""show = M[M.target.isin(['macrophage', 'myeloid', 'monocyte'])].gene.tolist()
fig, axes = plt.subplots(1, 4, figsize=(17, 14), sharey=True)
for k, (c, src) in enumerate([('sepsis', 'real'), ('sepsis', 'k1000'), ('covid', 'real'), ('covid', 'k1000')]):
    d = res[(res.cohort == c) & (res.source == src)].pivot(index='gene', columns='group', values='det').reindex(index=show, columns=GROUPS)
    ax = axes[k]; im = ax.imshow(d.to_numpy(), aspect='auto', cmap='Blues', vmin=0, vmax=1)
    ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels(GROUPS, rotation=60, ha='right', fontsize=10)
    ax.set_yticks(range(len(show))); ax.set_yticklabels(show, fontsize=9)
    n = res[(res.cohort == c) & (res.source == src)].drop_duplicates('group').set_index('group').n_cells.reindex(GROUPS)
    ax.set_title(f'{c} · {\"real\" if src == \"real\" else \"generated K=1,000\"}\\nMacrophage cells: {int(n[\"Macrophage\"])}', fontsize=11)
    for i in range(d.shape[0]):
        if d.iloc[i].isna().all(): ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(GROUPS), 1, color='#dddddd'))
    for j in range(len(GROUPS)):
        if d.iloc[:, j].isna().all(): ax.add_patch(plt.Rectangle((j - 0.5, -0.5), 1, len(show), color='#f0f0f0'))
    ax.axhline((M.target == 'macrophage').sum() - 0.5, color='k', lw=1); ax.axhline((M.target.isin(['macrophage', 'myeloid'])).sum() - 0.5, color='k', lw=1)
fig.colorbar(im, ax=axes, shrink=0.4, label='fraction of cells detected (count > 0)')
fig.savefig(f'{OUT}/macrophage_marker_detection_heatmap.png', dpi=130, bbox_inches='tight'); plt.show()""")

code("""obj = []
for c, (srcs, pos) in data.items():
    gg = [x for x in M.gene if x in pos.index]
    for src, (X, g) in srcs.items():
        a = ad.AnnData(sp.csr_matrix(X[:, pos.loc[gg].to_numpy()]), obs=pd.DataFrame({'cohort': c, 'source': src, 'group': pd.Series(g).astype(str).values}), var=pd.DataFrame(index=gg))
        a.obs_names = [f'{c}_{src}_{i}' for i in range(a.n_obs)]; obj.append(a)
A = ad.concat(obj, join='outer', fill_value=0); A.obs = A.obs.astype('category'); A.write_h5ad(f'{OUT}/macrophage_marker_counts.h5ad', compression='gzip'); A""")

nb = nbf.v4.new_notebook(); nb.cells = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
NotebookClient(nb, timeout=3600, kernel_name='python3').execute()
nbf.write(nb, f'{OUT}/macrophage_markers.ipynb'); print('done')
