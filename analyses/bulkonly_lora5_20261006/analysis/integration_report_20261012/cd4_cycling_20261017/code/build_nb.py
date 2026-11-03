import nbformat as nbf
C = []
def md(s): C.append(nbf.v4.new_markdown_cell(s))
def code(s): C.append(nbf.v4.new_code_cell(s))

md("""# COVID (Stephenson 2025): why PBISC output has no *T CD4 cycling* cells — 2026-10-17

Question: in the COVID (Imperial ESKD) analysis, real cells contain a **T CD4 cycling** population, but no PBISC-generated cell carries that label.
Is it (1) the 6k gene panel, (2) the labeller / label-transfer route, or (3) the generated cells themselves?

Inputs (all on the PC, extracted by `code/extract_cd4cyc.py`, h5py + numpy only):
- `recap`: Stephenson 2025 recap object (578k–588k cells, full transcriptome) — labels, 219-gene matrix for all cells, full rows for all T CD4 cycling cells + 6,000 other CD4 T + 3,000 others.
- `real_heldout`: real held-out panel counts (35 patients / 97 samples / 333,503 cells; 0 overlap with LoRA training patients) — labels, 219-gene matrix, full panel rows for selected cells.
- `gen_k1000`, `gen_k4000`: PBISC bulk-only COVID sets (410,000 cells each, 164 samples × 5 draws × 500 cells) — frozen AIFI labels, 219-gene matrix, full panel rows for selected cells.
- Existing results reused: R1/R3 mapping (`proportion_qualification`), patient-disjoint mapping check (`paired_patient_validation_20261014`), paired-test calibration (`features/imperial_eskd/calibration.csv`, K1k & K4k).
""")

code("""import numpy as np, pandas as pd, json, scipy.sparse as sp, collections, os, warnings
import scanpy as sc, anndata as ad, h5py
import matplotlib; import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
warnings.filterwarnings('ignore')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':110})
os.makedirs('results', exist_ok=True); os.makedirs('figs', exist_ok=True); os.makedirs('data', exist_ok=True)
GS = json.load(open('code/gene_set.json'))
def load(s):
    o = np.load(f'_extract/{s}_obs.npz', allow_pickle=True); G = np.load(f'_extract/{s}_G.npz', allow_pickle=True)
    return {k: o[k] for k in o.files}, G['G'], list(G['genes']), G['present']
def rows(s):
    r = np.load(f'_extract/{s}_rows.npz', allow_pickle=True)
    X = sp.csr_matrix((r['data'], r['indices'], r['indptr']), shape=tuple(r['shape']))
    return X, r['rows'], pd.Index(r['var']).astype(str)
oc, Gc, genes, pc = load('recap'); orr, Gr, _, pr = load('real_heldout')
og1, Gg1, _, pg = load('gen_k1000'); og4, Gg4, _, _ = load('gen_k4000')
gi = {g: i for i, g in enumerate(genes)}
_, _, var_gen = rows('gen_k1000'); PANEL = set(var_gen)
print('cells: recap', len(oc['lib']), '| real held-out', len(orr['lib']), '| K1k', len(og1['lib']), '| K4k', len(og4['lib']))
print('generated panel genes (unique symbols):', len(PANEL))""")

md("## 1. What is *T CD4 cycling* in the real data?")
code("""dep = oc['cell_type_detailed']; rec = oc['cell_type_detailed_recap']; cc = oc['case_control']
n_tot = len(dep); roster = None
t = pd.DataFrame({
  'label set': ['deposited (authors)', 'recap (our re-annotation)', 'real held-out (recap labels)'],
  'T CD4 cycling cells': [(dep=='T CD4 cycling').sum(), (rec=='T CD4 cycling').sum(), (orr['major']=='T CD4 cycling').sum()],
  'all cells': [n_tot, n_tot, len(orr['major'])]})
t['% of all cells'] = (t['T CD4 cycling cells']/t['all cells']*100).round(3)
cd4set = ['T CD4 Helper','T CD4 Naive','T CD4 CTL','T reg','T CD4 cycling']
t['% of CD4 T'] = [round((dep=='T CD4 cycling').sum()/np.isin(dep,cd4set).sum()*100,3), round((rec=='T CD4 cycling').sum()/np.isin(rec,cd4set).sum()*100,3),
                   round((orr['major']=='T CD4 cycling').sum()/np.isin(orr['major'],cd4set).sum()*100,3)]
display(t)
both = ((dep=='T CD4 cycling')&(rec=='T CD4 cycling')).sum()
print('cells labelled cycling in BOTH label sets:', both)
print('deposited cycling -> recap label:', collections.Counter(rec[dep=='T CD4 cycling']).most_common(5))
print('by case_control (deposited):', dict(collections.Counter(cc[dep=='T CD4 cycling'])), '| (recap):', dict(collections.Counter(cc[rec=='T CD4 cycling'])))
mk = ['MKI67','TOP2A','TYMS','BIRC5','CENPF','STMN1','CD4','CD40LG','FOXP3','CD8A','CD8B','NKG7']
rows_ = []
for name, m in [('deposited cycling', dep=='T CD4 cycling'), ('recap cycling', rec=='T CD4 cycling'),
                ('other CD4 T (recap)', np.isin(rec, cd4set[:-1]) & (dep!='T CD4 cycling'))]:
    r = {'group': name, 'n': int(m.sum()), 'median UMI': int(np.median(oc['lib'][m]))}
    r.update({g: round((Gc[m, gi[g]]>0).mean(), 2) for g in mk}); rows_.append(r)
t1 = pd.DataFrame(rows_); display(t1); t.to_csv('results/t1_real_label_counts.csv', index=False); t1.to_csv('results/t1b_real_label_detection.csv', index=False)""")
md("""The class is tiny (0.03–0.09% of all cells) and the two label sets agree on only a minority of cells. The recap label is MKI67-positive by construction; the deposited label is broader (MKI67+ 39%, 3.7× the library size of other CD4 T).""")

md("## 2. Is it the 6k panel?")
code("""S = [g for g in GS['S'] if pc[gi[g]]]; G2M = [g for g in GS['G2M'] if pc[gi[g]]]
cov = pd.DataFrame({'set': ['S phase (Tirosh/Seurat)', 'G2/M (Tirosh/Seurat)'],
  'genes in recap': [len(S), len(G2M)], 'in PBISC panel': [sum(g in PANEL for g in S), sum(g in PANEL for g in G2M)]})
cov['fraction'] = (cov['in PBISC panel']/cov['genes in recap']).round(2); display(cov)
print('in panel  :', [g for g in S+G2M if g in PANEL])
print('not in panel (examples):', [g for g in ['CDK1','CCNB1','CCNB2','HMGB2','TUBB','PCNA','MCM2','TK1','AURKA','AURKB','NUSAP1','CKS1B'] if g not in PANEL])

# full-transcriptome DE: recap T CD4 cycling vs other CD4 T (recap rows)
X, rws, var = rows('recap'); r_ = oc['cell_type_detailed_recap'][rws]; d_ = oc['cell_type_detailed'][rws]; pid = oc['patient_id'][rws]
de_res = {}; clf_res = []
for lab, arr in [('recap', r_), ('deposited', d_)]:
    y = arr=='T CD4 cycling'; keep = y | (np.isin(r_, cd4set[:-1]) & (d_!='T CD4 cycling') & (r_!='T CD4 cycling'))
    A = ad.AnnData(X[keep].astype(np.float32), obs=pd.DataFrame({'grp': np.where(y[keep],'cyc','cd4')}), var=pd.DataFrame(index=var)); A.var_names_make_unique()
    det = lambda g: np.asarray((A.X[(A.obs.grp==g).values]>0).mean(0)).ravel()
    dc, dr = det('cyc'), det('cd4')
    sc.pp.normalize_total(A, 1e4); sc.pp.log1p(A)
    sc.tl.rank_genes_groups(A, 'grp', groups=['cyc'], reference='cd4', method='wilcoxon')
    df = sc.get.rank_genes_groups_df(A, 'cyc').set_index('names')
    df['det_cyc'] = pd.Series(dc, index=A.var_names); df['det_cd4'] = pd.Series(dr, index=A.var_names); df = df.reset_index()
    df['in_panel'] = df.names.isin(PANEL)
    up = df[(df.logfoldchanges>1)&(df.pvals_adj<0.05)&(df.det_cyc>=0.25)&(df.det_cyc-df.det_cd4>=0.2)].sort_values('scores', ascending=False)
    bg = pd.Index(A.var_names)[dr>=0.01]
    de_res[lab] = {'n_cycling': int(y.sum()), 'markers': len(up), 'markers in panel': round(up.in_panel.mean(),3),
                   'top20 in panel': round(up.head(20).in_panel.mean(),2), 'top50 in panel': round(up.head(50).in_panel.mean(),2),
                   'background (genes det. >=1% CD4 T) in panel': round(np.mean(bg.isin(list(PANEL))),3)}
    df.to_csv(f'results/de_{lab}_cd4cycling_vs_cd4.csv', index=False)
    if lab == 'recap':
        print('top 30 recap markers ((x) = not in panel):', ', '.join(f"{a}{'' if b else '(x)'}" for a, b in zip(up.head(30).names, up.head(30).in_panel)))
    # panel-only classifier, grouped CV by patient
    Xk = X[keep]; yk = y[keep]; g = pid[keep]; lib = np.asarray(Xk.sum(1)).ravel(); N = sp.diags(1e4/lib) @ Xk; N.data = np.log1p(N.data); N = N.tocsr()
    for nm, cols in [('all genes', np.ones(len(var), bool)), ('panel genes', var.isin(list(PANEL))), ('panel genes without MKI67', var.isin(list(PANEL-{'MKI67'})))]:
        Nc = N[:, np.nonzero(cols)[0]]; p = np.zeros(len(yk))
        for tr, te in GroupKFold(5).split(Nc, yk, g):
            p[te] = LogisticRegression(C=0.1, max_iter=2000, class_weight='balanced').fit(Nc[tr], yk[tr]).predict_proba(Nc[te])[:, 1]
        clf_res.append({'label set': lab, 'features': nm, 'n genes': int(cols.sum()), 'AUROC': round(roc_auc_score(yk, p),3), 'AP': round(average_precision_score(yk, p),3)})
t2 = pd.DataFrame(de_res).T; display(t2); t3 = pd.DataFrame(clf_res); display(t3)
t2.to_csv('results/t2_marker_panel_coverage.csv'); t3.to_csv('results/t3_panel_only_classifier.csv', index=False)""")
md("""Panel verdict: cell-cycle genes are covered at the background rate (~30%), with the most specific ones (MKI67, TOP2A, TYMS, BIRC5, CENPF, PCLAF) present; several highly expressed cycle genes (CDK1, CCNB1, HMGB2, TUBB, PCNA) are absent. The DE marker list is *depleted* in the panel, but the missing genes are mostly housekeeping / metabolic genes that rise with activation. A panel-only classifier separates the class as well as the full transcriptome (AUROC 0.99–1.00) → **the panel is not the cause**.""")

md("## 3. Is it the labeller or the label-transfer route?")
code("""l2 = og1['aifi_l2']
print('PBISC native labeller (frozen AIFI L2) classes:', sorted(set(l2)))
lab_tab = pd.DataFrame({k: {'Proliferating T cell': (o['aifi_l2']=='Proliferating T cell').sum(), 'Proliferating NK cell': (o['aifi_l2']=='Proliferating NK cell').sum(), 'all cells': len(o['aifi_l2'])} for k, o in [('K1k', og1), ('K4k', og4)]}).T
lab_tab['% Prolif T'] = (lab_tab['Proliferating T cell']/lab_tab['all cells']*100).round(2); display(lab_tab)
rc = pd.read_csv('../../paired_patient_validation_20261014/results/covid_mapping_recoverability.csv')
r1 = rc[(rc.label.isin(['T CD4 cycling','B_ASC_dividing']))&(rc.source=='realT2')&(rc.level=='fine')][['protocol','route','label','recall','precision','n_true_cells']]
display(r1)
rp = pd.read_csv('../results/rare_parent_denominators.csv'); rp = rp[(rp.cohort=='covid')&(rp.subset.isin(['T CD4 cycling','B_ASC_dividing']))]
rp['expected K1k cells at real share'] = (rp.real_pooled_pct/100*rp.k1000_n_parent).round(1)
display(rp[['subset','group','real_n_target','real_n_parent','real_pooled_pct','k1000_n_target','k1000_n_parent','expected K1k cells at real share','k4000_n_target']])
r1.to_csv('results/t4_mapping_recovery.csv', index=False); rp.to_csv('results/t4b_r1_counts.csv', index=False)""")
md("""- The PBISC native labeller has no CD4-specific cycling class, only **Proliferating T cell** (CD4 and CD8 together), and it *does* assign that label to generated cells (~0.8%).
- The label-transfer routes cannot recover *T CD4 cycling* even from **real** cells in PBISC read-out: R1/R2 recall 0, R3 recall 0.04–0.09 (patient-disjoint, 9–19 true cells), while the dividing plasmablast class next to it is recovered (recall 0.79–0.87).
- So the zero in the integration report (R1) is expected **whatever PBISC generates**: at the real share, ~22 K1k cells would be expected, and R1 would place them in T CD4 Helper.""")

md("## 4. Do generated cells carry the cell-cycle programme? (label-free gate)")
code("""C14 = ['MKI67','TYMS','CENPF','TOP2A','BIRC5','CENPM','TPX2','PCLAF','CDC20','MYBL2','CENPE','UBE2C','ASPM','DLGAP5']
de_rec = pd.read_csv('results/de_recap_cd4cycling_vs_cd4.csv').set_index('names')
display(de_rec.loc[C14, ['det_cyc','det_cd4','logfoldchanges']].round(3).rename(columns={'det_cyc':'detected in recap cycling','det_cd4':'detected in other CD4 T'}))
idx = [gi[g] for g in C14]
Tlab = ['T CD4 Helper','T CD4 Naive','T CD4 CTL','T reg','T CD8 CTL','T CD8 CM','T CD8 Naive','T CD4 cycling','NKT','MAIT','T g/d','T CD8 EM','T CD8 EMRA']
NKlab = ['NK','NK CD56 bright']; ASC = ['B_ASC_IgA','B_ASC_IgG','B_ASC_IgM','B_ASC_dividing','B_ASC_other']
maj = orr['major']
rng = np.random.default_rng(20261017)
def gate(M): return (M>0).sum(1) >= 2
def lineage_masks_real(): return {'T': np.isin(maj, Tlab), 'NK': np.isin(maj, NKlab), 'plasmablast': np.isin(maj, ASC)}
def lineage_masks_gen(o): return {'T': o['celltype_l1_pred']=='T cell', 'NK': o['celltype_l1_pred']=='NK cell', 'plasmablast': o['aifi_l2']=='Plasma cell'}
Lr = lineage_masks_real(); res = []; gated = {}
gated['real'] = gate(Gr[:, idx])
for name, o, G in [('K1k', og1, Gg1), ('K4k', og4, Gg4)]:
    Lg = lineage_masks_gen(o); gated[name] = gate(G[:, idx])
    for lin in ['T','NK','plasmablast']:
        p = min(1.0, np.median(o['lib'][Lg[lin]]) / np.median(orr['lib'][Lr[lin]]))
        thin = rng.binomial(Gr[Lr[lin]][:, idx].astype(np.int64), p)
        res.append({'lineage': lin, 'comparison': name, 'real %': gated['real'][Lr[lin]].mean()*100,
                    f'real depth-matched %': gate(thin).mean()*100, 'generated %': gated[name][Lg[lin]].mean()*100,
                    'depth factor': round(p,3), 'n real': int(Lr[lin].sum()), 'n generated': int(Lg[lin].sum())})
t5 = pd.DataFrame(res); t5['generated / depth-matched real'] = (t5['generated %']/t5['real depth-matched %']).round(2)
display(t5.round(3)); t5.to_csv('results/t5_cycling_gate_by_lineage.csv', index=False)
print('real gated cells by label:', collections.Counter(maj[gated['real']]).most_common(10))""")
code("""# by disease group (per-sample mean of % cycling among T cells) and by draw
def per_sample(o, T, g):
    d = pd.DataFrame({'s': o['sample_id'], 'grp': o['group'], 'T': T, 'cyc': T & g})
    p = d.groupby(['s','grp'])[['T','cyc']].sum().reset_index(); p['pct'] = p.cyc/p['T']*100; return p
ps = {'real': per_sample(orr, Lr['T'], gated['real'])}
for name, o in [('K1k', og1), ('K4k', og4)]: ps[name] = per_sample(o, lineage_masks_gen(o)['T'], gated[name])
t6 = pd.concat({k: v.groupby('grp').pct.agg(['mean','median','count']) for k, v in ps.items()}, axis=1).round(3); display(t6)
for name, o in [('K1k', og1), ('K4k', og4)]:
    T = lineage_masks_gen(o)['T']; d = pd.DataFrame({'draw': o['draw'], 'T': T, 'c': T & gated[name]}).groupby('draw').sum(); d['%'] = (d.c/d['T']*100).round(3)
    print(name, 'cycling % of T by draw:', d['%'].to_dict())
t6.to_csv('results/t6_cycling_by_group.csv')""")

md("## 5. Which T cells cycle? CD4 vs CD8 identity of the cycling cells")
code("""def lin_call(G):
    cd8 = (G[:, gi['CD8A']]>0) | (G[:, gi['CD8B']]>0); cd4 = (G[:, gi['CD4']]>0) | (G[:, gi['CD40LG']]>0) | (G[:, gi['FOXP3']]>0)
    return np.where(cd8 & ~cd4, 'CD8', np.where(cd4 & ~cd8, 'CD4', np.where(cd4 & cd8, 'both', 'neither')))
rows_l = []
for name, o, G, T in [('real', orr, Gr, Lr['T']), ('K1k', og1, Gg1, lineage_masks_gen(og1)['T']), ('K4k', og4, Gg4, lineage_masks_gen(og4)['T'])]:
    L = lin_call(G); g = gated[name]
    for subset, m in [('all T', T), ('cycling T', T & g)]:
        c = collections.Counter(L[m]); n = m.sum()
        rows_l.append({'set': name, 'cells': subset, 'n': int(n), **{k: round(c.get(k,0)/n*100,1) for k in ['CD4','CD8','both','neither']}})
t7 = pd.DataFrame(rows_l); display(t7); t7.to_csv('results/t7_cycling_T_lineage.csv', index=False)
for name in ['real','K1k','K4k']:
    a = t7[(t7.set==name)&(t7.cells=='all T')].iloc[0]; b = t7[(t7.set==name)&(t7.cells=='cycling T')].iloc[0]
    print(f"{name}: CD4-marked share all T {a.CD4}% -> cycling T {b.CD4}% (enrichment {b.CD4/a.CD4:.2f}x); CD8 {a.CD8}% -> {b.CD8}%")
# CD4-marked cycling T as % of all T
for name, o, G, T in [('real', orr, Gr, Lr['T']), ('K1k', og1, Gg1, lineage_masks_gen(og1)['T']), ('K4k', og4, Gg4, lineage_masks_gen(og4)['T'])]:
    L = lin_call(G); print(name, 'CD4-marked cycling cells per 100 T:', round(((L=='CD4') & gated[name] & T).sum()/T.sum()*100, 3))""")

md("## 6. Programme of the cycling cells: cell cycle kept, CD4 activation identity lost")
code("""Xr, rr, varr = rows('real_heldout'); Xg, rg, varg = rows('gen_k1000')
shared = [g for g in varr.unique() if g in set(varg)]
def sub(X, var, gl):
    f = pd.Series(np.arange(len(var)), index=var); f = f[~f.index.duplicated()]; return X[:, f[gl].values].tocsr()
Xr_s = sub(Xr, varr, shared); Xg_s = sub(Xg, varg, shared); si = {g: i for i, g in enumerate(shared)}
gate_rows = lambda X: (X[:, [si[g] for g in C14]]>0).sum(1).A1 >= 2
majr = orr['major'][rr]; l2r = og1['aifi_l2'][rg]; l1r = og1['celltype_l1_pred'][rg]
Tr = np.isin(majr, [x for x in Tlab if x not in ['NKT','MAIT','T g/d','T CD8 EM','T CD8 EMRA']])
gr_, gg_ = gate_rows(Xr_s), gate_rows(Xg_s)
grp = {'real cycling T': Tr & (gr_ | (majr=='T CD4 cycling')), 'real other T': Tr & ~gr_ & (majr!='T CD4 cycling'),
       'gen cycling T (gated)': ((l1r=='T cell') | (l2r=='Proliferating T cell')) & gg_,
       'gen "Proliferating T" not gated': (l2r=='Proliferating T cell') & ~gg_,
       'gen other T': (l1r=='T cell') & (l2r!='Proliferating T cell') & ~gg_}
print({k: int(v.sum()) for k, v in grp.items()})
def lmean(X, m):
    lib = np.asarray(X[m].sum(1)).ravel(); N = sp.diags(1e4/lib) @ X[m]; return np.log1p(np.asarray(N.mean(0)).ravel())
mu = {k: lmean(Xr_s if k.startswith('real') else Xg_s, m) for k, m in grp.items()}
lib_med = {k: int(np.median(np.asarray((Xr_s if k.startswith('real') else Xg_s)[m].sum(1)).ravel())) for k, m in grp.items()}
print('median panel UMI:', lib_med)
sig = pd.DataFrame({'real: cycling vs other T': mu['real cycling T']-mu['real other T'],
                    'gen: cycling (gated) vs other T': mu['gen cycling T (gated)']-mu['gen other T'],
                    'gen: "Proliferating T" not gated vs other T': mu['gen "Proliferating T" not gated']-mu['gen other T']}, index=shared)
CYC = set(GS['S'])|set(GS['G2M'])|set(GS['EXTRA'])
ISG = ['ISG15','IFI6','MX1','IFI44L','XAF1','OAS2','EPSTI1','IFI27','IFIT3','IFI44','OAS1','RSAD2','IFIT1','STAT1','IFI16','IFITM3','LY6E','MX2','HERC5','IRF7']
CD4ACT = ['CTLA4','DUSP4','TRIB2','CD27','CXCR3','CXCL13','ITM2A','SIRPG','FOXP3','ICOS','IL2RA','TNFRSF4','TNFRSF18','MAF','TIGIT','BATF','CD40LG']
sig['category'] = ['cell cycle' if g in CYC else 'CD4 activation / Treg-Tfh' if g in CD4ACT else 'interferon' if g in ISG else 'other' for g in shared]
sig.to_csv('results/t8_signature_shift_log1p.csv')
summ = sig.groupby('category').mean(numeric_only=True).round(2); summ['n genes'] = sig.category.value_counts(); display(summ)
key = ['MKI67','TOP2A','TYMS','CENPF','BIRC5','STMN1','CTLA4','DUSP4','TRIB2','CD27','CXCR3','CXCL13','SIRPG','FOXP3','ICOS','ITM2A','IFI27','ISG15','HLA-DRA','LGALS1']
display(sig.loc[[k for k in key if k in sig.index]].round(2))
summ.to_csv('results/t8b_signature_by_category.csv')""")

md("## 7. What does the labeller call *Proliferating T* in generated cells?")
code("""rows_p = []
for name, o, G in [('K1k', og1, Gg1), ('K4k', og4, Gg4)]:
    for lab in ['Proliferating T cell','Proliferating NK cell','Plasma cell']:
        m = o['aifi_l2']==lab
        rows_p.append({'set': name, 'AIFI L2 label': lab, 'n': int(m.sum()), 'median UMI': int(np.median(o['lib'][m])),
                       'median confidence': round(float(np.median(o['aifi_l2_confidence'][m].astype(float))),2),
                       'MKI67+': round((G[m, gi['MKI67']]>0).mean(),3), 'cycle gate +': round(gated[name][m].mean(),3)})
    T = lineage_masks_gen(o)['T']; rows_p.append({'set': name, 'AIFI L2 label': 'all generated T', 'n': int(T.sum()), 'median UMI': int(np.median(o['lib'][T])),
        'median confidence': np.nan, 'MKI67+': round((G[T, gi['MKI67']]>0).mean(),4), 'cycle gate +': round(gated[name][T].mean(),4)})
rows_p.append({'set': 'real', 'AIFI L2 label': 'T CD4 cycling (recap)', 'n': int((maj=='T CD4 cycling').sum()), 'median UMI': int(np.median(orr['lib'][maj=='T CD4 cycling'])),
               'median confidence': np.nan, 'MKI67+': round((Gr[maj=='T CD4 cycling', gi['MKI67']]>0).mean(),3), 'cycle gate +': round(gated['real'][maj=='T CD4 cycling'].mean(),3)})
t9 = pd.DataFrame(rows_p); display(t9); t9.to_csv('results/t9_generated_prolif_labels.csv', index=False)
top = sig.sort_values('gen: "Proliferating T" not gated vs other T', ascending=False).head(15)
display(top.round(2))""")

md("## 8. Format-aware classifier (R3-like) for *T CD4 cycling*, applied to generated T cells")
code("""cols = np.nonzero(pr)[0]; T = Lr['T']; p = min(1.0, np.median(og1['lib'][lineage_masks_gen(og1)['T']])/np.median(orr['lib'][T]))
M = rng.binomial(Gr[T][:, cols].astype(np.int64), p).astype(np.float32); lib = orr['lib'][T]*p
F = np.log1p(M/lib[:, None]*1e4); y = maj[T]=='T CD4 cycling'; g = orr['patient_id'][T]
pp = np.zeros(len(y))
for tr, te in GroupKFold(5).split(F, y, g):
    pp[te] = LogisticRegression(C=0.1, max_iter=3000, class_weight='balanced').fit(F[tr], y[tr]).predict_proba(F[te])[:, 1]
clf = LogisticRegression(C=0.1, max_iter=3000, class_weight='balanced').fit(F, y)
out = [{'set': 'real held-out (CV, depth-matched)', 'called % @0.9': round(np.mean(pp>=0.9)*100,3), 'true label %': round(y.mean()*100,3),
        'recall': round(np.mean(pp[y]>=0.9),2), 'precision': round(np.mean(y[pp>=0.9]),2), 'AUROC': round(roc_auc_score(y, pp),3)}]
for name, o, G in [('K1k', og1, Gg1), ('K4k', og4, Gg4)]:
    Tg = lineage_masks_gen(o)['T']; q = clf.predict_proba(np.log1p(G[Tg][:, cols]/o['lib'][Tg][:, None]*1e4))[:, 1]
    out.append({'set': name, 'called % @0.9': round(np.mean(q>=0.9)*100,3)})
t10 = pd.DataFrame(out); display(t10); t10.to_csv('results/t10_format_aware_classifier.csv', index=False)
coef = pd.Series(clf.coef_[0], index=np.array(genes)[cols]).sort_values(); print('top +:', coef.tail(8).round(2).to_dict()); print('top -:', coef.head(5).round(2).to_dict())""")

md("## 9. Paired-test calibration (same samples, generated vs observed)")
code("""cal = []
for k in ['k1000','k4000']:
    c = pd.read_csv(f'inputs/calibration_imperial_{k}.csv'); c = c[(c.draw=='pooled')&(c.feature.isin(['l2 | Proliferating T cell','Proliferating T (label)','Cycling T/NK','Cycling activated T (CD38/HLA-DR)','l2 | Proliferating NK cell']))]
    c = c[['feature','gen_mean','obs_mean','ratio','spearman','COVID+ vs COVID-: gen diff','COVID+ vs COVID-: obs diff']].copy(); c.insert(0,'K',k); cal.append(c)
t11 = pd.concat(cal); display(t11.round(4)); t11.to_csv('results/t11_paired_calibration.csv', index=False)""")
md("""On the paired test samples, the labeller's *Proliferating T* fraction is at or above the observed level (ratio 1.3), but the marker-defined cycling T/NK state is at 0.41–0.44× and **cycling activated (CD38/HLA-DR) T at 0.14–0.16×** — the same picture as above: labels over-count, the programme is under-produced, and activated cycling T cells are the most depleted.""")

md("## 10. Figures")
code("""BL, OR, GR, DK = '#2A78D6', '#EB6834', '#9AA3B2', '#4A4F5C'
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for a, lin in zip(ax, ['T','NK','plasmablast']):
    r = t5[t5.lineage==lin]
    vals = [r['real %'].iloc[0], r['real depth-matched %'].mean(), r[r.comparison=='K1k']['generated %'].iloc[0], r[r.comparison=='K4k']['generated %'].iloc[0]]
    a.bar(['real', 'real\\n(depth-\\nmatched)', 'PBISC\\nK1k', 'PBISC\\nK4k'], vals, color=[DK, GR, BL, OR])
    for i, v in enumerate(vals): a.text(i, v, f'{v:.2f}%', ha='center', va='bottom', fontsize=9)
    a.set_title(f'{lin}: cells with \\u22652 of 14 cycle genes'); a.set_ylabel('% of lineage')
plt.tight_layout(); plt.savefig('figs/f1_cycling_gate_by_lineage.png', dpi=200); plt.savefig('figs/f1_cycling_gate_by_lineage.pdf'); plt.show()

fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
sets = ['real','K1k','K4k']; cols_ = {'CD4': BL, 'CD8': OR, 'both': DK, 'neither': GR}
for a, cells in zip(ax, ['all T', 'cycling T']):
    left = np.zeros(3)
    for k in ['CD4','CD8','both','neither']:
        v = [t7[(t7.set==s)&(t7.cells==cells)][k].iloc[0] for s in sets]
        a.barh(['real','PBISC K = 1,000','PBISC K = 4,000'], v, left=left, color=cols_[k], label=k); left += v
    a.set_title(f'{cells}: marker-called lineage (%)'); a.set_xlim(0, 100); a.invert_yaxis()
ax[1].legend(frameon=False, bbox_to_anchor=(1.0, 1.0), loc='upper left')
plt.tight_layout(); plt.savefig('figs/f2_cycling_T_lineage.png', dpi=200); plt.savefig('figs/f2_cycling_T_lineage.pdf'); plt.show()

fig, ax = plt.subplots(1, 2, figsize=(12, 5), sharex=True, sharey=True)
catc = {'cell cycle': OR, 'CD4 activation / Treg-Tfh': BL, 'interferon': '#7A5CC2', 'other': '#D5D9E0'}
for a, col in zip(ax, ['gen: cycling (gated) vs other T', 'gen: "Proliferating T" not gated vs other T']):
    for cat in ['other','interferon','cell cycle','CD4 activation / Treg-Tfh']:
        s_ = sig[sig.category==cat]; a.scatter(s_['real: cycling vs other T'], s_[col], s=10 if cat=='other' else 26, c=catc[cat], label=cat, alpha=0.8 if cat!='other' else 0.5, lw=0)
    for gname in (['MKI67','TOP2A','TYMS','CTLA4','FOXP3','CD27','IFI27','STMN1'] if col.startswith('gen: cycling') else ['MKI67','IFI27','STMN1','CTLA4']):
        if gname in sig.index: a.annotate(gname, (sig.loc[gname,'real: cycling vs other T'], sig.loc[gname, col]), fontsize=8, xytext=(3, 2), textcoords='offset points')
    a.axhline(0, c='k', lw=0.5); a.axvline(0, c='k', lw=0.5); a.plot([-1, 3], [-1, 3], ls='--', c=GR, lw=0.8)
    a.set_xlabel('real: cycling T vs other T (\\u0394 log1p CP10k)'); a.set_title(col.replace('gen: ', 'PBISC K = 1,000: '), fontsize=10)
ax[0].set_ylabel('generated (\\u0394 log1p CP10k)'); ax[1].legend(frameon=False, fontsize=8, loc='lower right')
plt.tight_layout(); plt.savefig('figs/f3_signature_shift.png', dpi=200); plt.savefig('figs/f3_signature_shift.pdf'); plt.show()""")

md("## 11. Save cells (h5ad) for inspection")
code("""obs_r = pd.DataFrame({'source': 'real_heldout', 'label': orr['major'][rr], 'group': orr['group'][rr], 'sample_id': orr['sample_id'][rr],
                      'cycle_gate': gr_, 'n_cycle_genes': (Xr_s[:, [si[g] for g in C14]]>0).sum(1).A1})
Xg4, rg4, varg4 = rows('gen_k4000'); Xg4_s = sub(Xg4, varg4, shared)
parts = [ad.AnnData(Xr_s.astype(np.float32), obs=obs_r)]
for name, o, Xs, rws_ in [('gen_k1000', og1, Xg_s, rg), ('gen_k4000', og4, Xg4_s, rg4)]:
    gg = (Xs[:, [si[g] for g in C14]]>0).sum(1).A1
    parts.append(ad.AnnData(Xs.astype(np.float32), obs=pd.DataFrame({'source': name, 'label': o['aifi_l2'][rws_], 'group': o['group'][rws_], 'sample_id': o['sample_id'][rws_],
                                                                      'cycle_gate': gg >= 2, 'n_cycle_genes': gg})))
for p_ in parts: p_.var_names = shared
A = ad.concat(parts, index_unique='-', merge='same')  # X = raw panel counts
A.uns['note'] = 'Selected cells only (all labelled cycling/prolif + random T/NK). X = raw panel counts on shared genes. cycle_gate = >=2 of 14 cycle genes.'
A.uns['cycle_genes'] = C14
A.write_h5ad('data/cd4_cycling_cells.h5ad', compression='gzip'); print(A)""")

md("## 12. Summary")
code("""summary = pd.DataFrame([
 ['Panel', 'Specific cycle genes present (MKI67, TOP2A, TYMS, BIRC5, CENPF); S/G2M coverage 30% = background; panel-only classifier AUROC ' + str(t3[(t3['label set']=='recap')&(t3.features=='panel genes')].AUROC.iloc[0]), 'not the cause'],
 ['Labeller', 'Native AIFI L2 has no CD4-cycling class (only Proliferating T); it assigns that label to ~0.8% of generated cells, mostly non-cycling ISG/HLA-DR cells', 'contributes'],
 ['Label transfer', 'Real T CD4 cycling recall in PBISC read-out: R1/R2 0, R3 0.04–0.09 (B_ASC_dividing 0.79–0.87); expected ~22 cells at real share', 'explains the zero'],
 ['Generation: amount', 'Cycling T (>=2 of 14 genes): generated ' + ', '.join(f"{x:.2f}%" for x in t5[t5.lineage=='T']['generated %']) + ' vs real depth-matched ' + ', '.join(f"{x:.2f}%" for x in t5[t5.lineage=='T']['real depth-matched %']), 'about one third'],
 ['Generation: identity', 'Real cycling T are CD4-dominant with CTLA4/FOXP3/CXCL13/DUSP4/CD27 programme; generated cycling T are CD8/unassigned and lack that programme though cycle genes are strong', 'main biological shortfall'],
], columns=['level', 'evidence', 'verdict']); display(summary); summary.to_csv('results/summary.csv', index=False)""")

nb = nbf.v4.new_notebook(); nb['cells'] = C
nb['metadata'] = {'kernelspec': {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}}
nbf.write(nb, '/home/claude/cd4cyc/cd4_cycling.ipynb')
print('ok', len(C))
