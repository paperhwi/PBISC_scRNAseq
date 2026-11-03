import nbformat as nbf
from nbclient import NotebookClient
OUT = "/mnt/user-data/outputs/integration_report_20261012/claim_classifiers"
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# New classifiers for Ye sepsis claims SEP-4A and SEP-3

- **SEP-4A**: Mono_c01-CD14-HLA-DRA as a share of CD14 monocytes (real, all ages, in-scheme: healthy 17.4% → sepsis 9.8%). The panel-matched mean-marker score reached AUC 0.49; the route-3 fine classifier recovered it with recall 0.26 and precision 0.28.
- **SEP-3**: NKT_c05-IFNG-CCL4hi as a share of NKT cells (healthy 3.3% → sepsis 22%). Generated cells had no usable NKT parent: route 3 recovers the NKT parent from PBISC-format real cells with recall 0.36 and precision 0.48.

**Approach**: format-aware, task-specific classifiers.
1. Train on real cells of half A of the samples passed through the PBISC-like read-out (realT2 A), true labels.
2. Choose the model by 5-fold cross-validation **grouped by sample** within A. Options: logistic regression over a C grid, or gradient boosting.
3. Test once on half B (realT2 B, different donors).
4. Estimate proportions by classify-and-count (CC) and adjusted classify-and-count, ACC = (CC − FPR) / (TPR − FPR), using the cross-validated TPR and FPR from A.
5. Apply the identical pipeline to generated cells. Also apply it to realT2 B, so generated and real are read out by the same instrument.

Features: log1p CP10k over the shared panel; 2,000 HVGs chosen within the parent on realT2 A; scaled with A mean and SD. Adults only for group contrasts (the generated set has one paediatric sample).""")

code("""import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, os, warnings, json
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
M = '/home/claude/work/mono'; P = '/home/claude/work/prop'; OUT = '%s'; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(21)
pd.set_option('display.width', 230); pd.set_option('display.float_format', '{:.3f}'.format)
r = ad.read_h5ad(f'{M}/sepsis_real_AB.h5ad'); g = ad.read_h5ad(f'{M}/sepsis_gen.h5ad')
r.obs['age'] = pd.read_csv(f'{M}/sepsis_real_AB_obs_age.csv.gz').age.astype(str).values
mapB = ad.read_h5ad(f'{P}/sepsis_realT2_mapping.h5ad', backed='r').obs
r.obs['R3_fine'] = 'NA'; r.obs['R3_lin6'] = 'NA'
b = (r.obs.half == 'B').to_numpy(); assert (r.obs.fine.values[b] == mapB.true_fine.astype(str).values).all()
r.obs.loc[b, 'R3_fine'] = mapB.R3_fine.astype(str).values; r.obs.loc[b, 'R3_lin6'] = mapB.R3_lin6.astype(str).values
grp_r = r.obs.group.map({'healthy control': 'HC', 'ICU control': 'ICU control', 'sepsis': 'sepsis'}).astype(str).to_numpy()
grp_g = g.obs.group.map({'adult_HC': 'HC', 'adult_ICU_ctrl': 'ICU control', 'adult_sepsis': 'sepsis', 'ped_ICU_ctrl': 'ped'}).astype(str).to_numpy()
adult_r = (r.obs.age == 'adult').to_numpy()
def lognorm(X):
    X = sp.csr_matrix(X, dtype=np.float32); Y = (sp.diags(1e4 / np.maximum(X.sum(1).A1, 1)) @ X).tocsr(); Y.data = np.log1p(Y.data); return Y
print(r, g)""" % OUT)

code("""class Task:
    '''binary classifier trained on realT2 A cells in `parent_A`, positives = `pos_labels`.'''
    def __init__(self, name, parent_A, y_A, groups_A):
        self.name = name
        XA = r.layers['t2'][parent_A]
        a = ad.AnnData(sp.csr_matrix(XA)); sc.pp.highly_variable_genes(a, flavor='seurat_v3', n_top_genes=2000)
        self.hv = np.where(a.var.highly_variable.to_numpy())[0]
        Y = lognorm(XA)[:, self.hv].toarray(); self.mu = Y.mean(0); self.sd = Y.std(0) + 1e-6
        self.XA = self.feat(XA); self.yA = y_A; self.gA = groups_A
    def feat(self, X): return np.clip((lognorm(X)[:, self.hv].toarray() - self.mu) / self.sd, -10, 10).astype(np.float32)
    def select(self):
        cv = GroupKFold(5); res = []
        cands = [('LR', C) for C in [0.003, 0.01, 0.03, 0.1, 0.3]] + [('HGB', None)]
        for kind, C in cands:
            p = np.zeros(len(self.yA))
            for tr, te in cv.split(self.XA, self.yA, self.gA):
                m = self.make(kind, C).fit(self.XA[tr], self.yA[tr]); p[te] = m.predict_proba(self.XA[te])[:, 1]
            res.append(dict(model=kind, C=C, cv_auc=roc_auc_score(self.yA, p), p=p))
        best = max(res, key=lambda d: d['cv_auc']); self.kind, self.C, self.pcv = best['model'], best['C'], best['p']
        fpr, tpr, thr = roc_curve(self.yA, self.pcv); j = np.argmax(tpr - fpr); self.thr = thr[j]
        pred = self.pcv >= self.thr; self.TPR = pred[self.yA == 1].mean(); self.FPR = pred[self.yA == 0].mean()
        self.model = self.make(self.kind, self.C).fit(self.XA, self.yA)
        self.cvtab = pd.DataFrame([{k: v for k, v in d.items() if k != 'p'} for d in res])
        return self.cvtab
    @staticmethod
    def make(kind, C):
        return LogisticRegression(C=C, max_iter=1000, class_weight='balanced') if kind == 'LR' else HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, class_weight='balanced', random_state=0)
    def predict(self, X): return self.model.predict_proba(self.feat(X))[:, 1] if X.shape[0] else np.zeros(0)
    def crossfit(self, parent_B, y_B, groups_B):
        XB = self.feat(r.layers['t2'][parent_B]); p = np.zeros(len(y_B))
        for tr, te in GroupKFold(5).split(XB, y_B, groups_B):
            p[te] = self.make(self.kind, self.C).fit(XB[tr], y_B[tr]).predict_proba(XB[te])[:, 1]
        fpr, tpr, thr = roc_curve(y_B, p); self.thrB = thr[np.argmax(tpr - fpr)]; self.cv_auc_B = roc_auc_score(y_B, p)
        self.modelB = self.make(self.kind, self.C).fit(XB, y_B)
    def call_real(self, idx):
        idx = np.asarray(idx); out = np.zeros(len(idx), bool); a = isA[idx]
        if (~a).any(): out[~a] = self.model.predict_proba(self.feat(r.layers['t2'][idx[~a]]))[:, 1] >= self.thr
        if a.any(): out[a] = self.modelB.predict_proba(self.feat(r.layers['t2'][idx[a]]))[:, 1] >= self.thrB
        return out
    def call_gen(self, X):
        if X.shape[0] == 0: return np.zeros(0)
        F = self.feat(X)
        return ((self.model.predict_proba(F)[:, 1] >= self.thr).astype(float) + (self.modelB.predict_proba(F)[:, 1] >= self.thrB).astype(float)) / 2
    def acc(self, cc): return float(np.clip((cc - self.FPR) / max(self.TPR - self.FPR, 1e-6), 0, 1))
    def test(self, XB, yB):
        p = self.predict(XB); pred = p >= self.thr
        return dict(task=self.name, model=f'{self.kind} C={self.C}', cv_auc_A=self.cvtab.cv_auc.max(), test_auc_B=roc_auc_score(yB, p),
                    recall_B=pred[yB == 1].mean(), precision_B=pred[yB == 1].sum() / max(pred.sum(), 1), threshold=self.thr, TPR_A=self.TPR, FPR_A=self.FPR)

def boot_frac(pos, smp, grp, groups=('HC', 'sepsis'), nb=1000, acc=None):
    '''pooled fraction per group with 95% CI by resampling samples; acc = function applied to the pooled CC fraction'''
    f = (lambda x: x) if acc is None else acc
    out = {}; per = pd.DataFrame({'pos': pos, 'smp': smp, 'grp': grp}).groupby(['grp', 'smp']).pos.agg(['sum', 'count']).reset_index()
    bs = {}
    for G in groups:
        d = per[per.grp == G]
        if len(d) == 0: out[G] = (np.nan, np.nan, np.nan, 0, 0); continue
        est = f(d['sum'].sum() / d['count'].sum())
        bb = np.array([f(d.iloc[ix]['sum'].sum() / d.iloc[ix]['count'].sum()) for ix in rng.integers(0, len(d), (nb, len(d)))]); bs[G] = bb
        out[G] = (est, np.percentile(bb, 2.5), np.percentile(bb, 97.5), len(d), int(d['count'].sum()))
    diff = bs[groups[1]] - bs[groups[0]] if all(G in bs for G in groups) else np.array([np.nan])
    return out, (out[groups[1]][0] - out[groups[0]][0], np.nanpercentile(diff, 2.5), np.nanpercentile(diff, 97.5))""")

md("""## SEP-4A: HLA-DR-high CD14 monocytes (Mono_c01) within CD14 monocytes""")

code("""fine = r.obs.fine.astype(str); isA = (r.obs.half == 'A').to_numpy(); isB = ~isA
cd14_true = fine.str.match(r'Mono_c0[1-6]').to_numpy()
t4 = Task('SEP-4A Mono_c01 vs other CD14 mono', np.where(isA & cd14_true)[0], (fine[isA & cd14_true] == 'Mono_c01-CD14-HLA-DRA').astype(int).to_numpy(),
          r.obs['sample'].to_numpy()[isA & cd14_true])
print(t4.select().to_string(index=False))
ixB = np.where(isB & cd14_true)[0]
T4 = t4.test(r.layers['t2'][ixB], (fine.to_numpy()[ixB] == 'Mono_c01-CD14-HLA-DRA').astype(int))
t4.crossfit(ixB, (fine.to_numpy()[ixB] == 'Mono_c01-CD14-HLA-DRA').astype(int), r.obs['sample'].to_numpy()[ixB]); T4['cv_auc_B_model'] = t4.cv_auc_B; pd.Series(T4)""")

code("""# Same instrument for real and generated: real cells are called out-of-sample (half A by the model trained on B, half B by the model trained on A);
# generated cells are called by both models and the two calls averaged. Real parent = true CD14 monocytes (route-3 CD14 recall/precision 0.98).
smp_r = r.obs['sample'].astype(str).to_numpy(); smp_g = g.obs['sample'].astype(str).to_numpy()
def cd14_true_like(lab): return pd.Series(lab).astype(str).str.match(r'Mono_c0[1-6]').to_numpy()
def collect(rows, contr, src, pos, smp, grp, **extra):
    o, d = boot_frac(pos, smp, grp)
    rows += [dict(source=src, group=G, est=v[0], lo=v[1], hi=v[2], samples=v[3], cells=v[4], **extra) for G, v in o.items()]
    contr.append(dict(source=src, diff=d[0], lo=d[1], hi=d[2], **extra))
rows, contr = [], []
m = cd14_true & adult_r
collect(rows, contr, 'real, true labels', fine.to_numpy()[m] == 'Mono_c01-CD14-HLA-DRA', smp_r[m], grp_r[m])
collect(rows, contr, 'realT2, classifier', t4.call_real(np.where(m)[0]), smp_r[m], grp_r[m])
for K in ['k1000', 'k4000']:
    m = (g.obs.K == K).to_numpy() & cd14_true_like(g.obs.R3_fine)
    collect(rows, contr, f'gen {K}, classifier', t4.call_gen(g.X[np.where(m)[0]]), smp_g[m], grp_g[m])
R4, C4 = pd.DataFrame(rows), pd.DataFrame(contr)
R4.to_csv(f'{OUT}/sep4a_fractions.csv', index=False); C4.to_csv(f'{OUT}/sep4a_contrast.csv', index=False)
print(R4.to_string(index=False)); print(C4.to_string(index=False))""")

code("""# threshold-free companion: MHC class II score per sample in CD14 monocytes (adult HC vs sepsis), effect in SD units
MHC = ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DMA', 'HLA-DMB', 'CD74']
gi_r = r.var_names.get_indexer([x for x in MHC if x in r.var_names])
def mhc_effect(X, smp, grp):
    Y = lognorm(X); s = np.asarray(Y[:, gi_r].mean(1)).ravel()
    d = pd.DataFrame({'s': s, 'smp': smp, 'grp': grp}).groupby(['grp', 'smp']).s.mean().reset_index()
    a, b = d[d.grp == 'HC'].s, d[d.grp == 'sepsis'].s
    return dict(HC=a.mean(), sepsis=b.mean(), d=(b.mean() - a.mean()) / np.sqrt((a.var() + b.var()) / 2), n_HC=len(a), n_sepsis=len(b))
rows = []
m = cd14_true & adult_r
rows.append(dict(source='real A+B (raw)', **mhc_effect(r.X[np.where(m)[0]], r.obs['sample'].to_numpy()[m], grp_r[m])))
rows.append(dict(source='realT2 A+B', **mhc_effect(r.layers['t2'][np.where(m)[0]], r.obs['sample'].to_numpy()[m], grp_r[m])))
for K in ['k1000', 'k4000']:
    m = (g.obs.K == K).to_numpy() & cd14_true_like(g.obs.R3_fine)
    rows.append(dict(source=f'gen {K}', **mhc_effect(g.X[np.where(m)[0]], g.obs['sample'].to_numpy()[m], grp_g[m])))
MH = pd.DataFrame(rows); MH.to_csv(f'{OUT}/sep4a_mhc2_score.csv', index=False); MH""")

md("""## SEP-3: IFNG/CCL4-high NKT (NKT_c05) within NKT cells

Two steps: (1) NKT parent among T/NK lymphocytes; (2) NKT_c05 within NKT. Real lymphocytes exclude 'Unresolved (low-quality T/NK)' and 'Platelet-T/NK complex'. Alternative target NKT_c04 + c05 (both IFNG) reported as well.""")

code("""lymph = (fine.str.match(r'^(T_|NK_|NKT_|gdT)')).to_numpy()
isnkt = fine.str.startswith('NKT_').to_numpy()
t3p = Task('SEP-3 step 1: NKT vs other T/NK', np.where(isA & lymph)[0], isnkt[isA & lymph].astype(int), r.obs['sample'].to_numpy()[isA & lymph])
print(t3p.select().to_string(index=False))
ixB = np.where(isB & lymph)[0]; T3p = t3p.test(r.layers['t2'][ixB], isnkt[ixB].astype(int)); t3p.crossfit(ixB, isnkt[ixB].astype(int), r.obs['sample'].to_numpy()[ixB]); T3p['cv_auc_B_model'] = t3p.cv_auc_B
r3_major_B = pd.Series(mapB.R3_major.astype(str).values)
y = isnkt[isB]; p = (r3_major_B == 'NKT cells').to_numpy()
print('R3 route, NKT major on realT2 B (reference): recall %.2f precision %.2f' % (p[y].mean(), y[p].mean()))
t3c = Task('SEP-3 step 2: NKT_c05 vs other NKT', np.where(isA & isnkt)[0], (fine[isA & isnkt] == 'NKT_c05-IFNG-CCL4hi').astype(int).to_numpy(), r.obs['sample'].to_numpy()[isA & isnkt])
print(t3c.select().to_string(index=False))
ixB = np.where(isB & isnkt)[0]; T3c = t3c.test(r.layers['t2'][ixB], (fine.to_numpy()[ixB] == 'NKT_c05-IFNG-CCL4hi').astype(int)); t3c.crossfit(ixB, (fine.to_numpy()[ixB] == 'NKT_c05-IFNG-CCL4hi').astype(int), r.obs['sample'].to_numpy()[ixB]); T3c['cv_auc_B_model'] = t3c.cv_auc_B
t3c2 = Task('SEP-3 alt: NKT_c04+c05 (IFNG) vs other NKT', np.where(isA & isnkt)[0], fine[isA & isnkt].isin(['NKT_c04-IFNG', 'NKT_c05-IFNG-CCL4hi']).astype(int).to_numpy(), r.obs['sample'].to_numpy()[isA & isnkt])
t3c2.select(); T3c2 = t3c2.test(r.layers['t2'][ixB], np.isin(fine.to_numpy()[ixB], ['NKT_c04-IFNG', 'NKT_c05-IFNG-CCL4hi']).astype(int)); t3c2.crossfit(ixB, np.isin(fine.to_numpy()[ixB], ['NKT_c04-IFNG', 'NKT_c05-IFNG-CCL4hi']).astype(int), r.obs['sample'].to_numpy()[ixB]); T3c2['cv_auc_B_model'] = t3c2.cv_auc_B
TT = pd.DataFrame([T4, T3p, T3c, T3c2]); TT.to_csv(f'{OUT}/classifier_performance.csv', index=False); TT""")

code("""lymph_like = lambda lab: pd.Series(lab).astype(str).isin(['T', 'NK']).to_numpy()
def readout_sep3(tt, labels, name):
    rows, contr = [], []
    m = isnkt & adult_r
    collect(rows, contr, 'real, true labels', np.isin(fine.to_numpy()[m], labels), smp_r[m], grp_r[m], target=name)
    m = lymph & adult_r; idx = np.where(m)[0]; par = t3p.call_real(idx); idx2 = idx[par]
    sh = pd.Series(par).groupby(grp_r[idx]).mean().to_dict()
    collect(rows, contr, 'realT2, classifier', tt.call_real(idx2), smp_r[idx2], grp_r[idx2], target=name, parent_share_HC=sh.get('HC'), parent_share_sepsis=sh.get('sepsis'))
    for K in ['k1000', 'k4000']:
        idx = np.where((g.obs.K == K).to_numpy() & lymph_like(g.obs.R3_lin6))[0]; par = t3p.call_gen(g.X[idx]) >= 0.5; idx2 = idx[par]
        sh = pd.Series(par).groupby(grp_g[idx]).mean().to_dict()
        collect(rows, contr, f'gen {K}, classifier', tt.call_gen(g.X[idx2]), smp_g[idx2], grp_g[idx2], target=name, parent_share_HC=sh.get('HC'), parent_share_sepsis=sh.get('sepsis'))
    return pd.DataFrame(rows), pd.DataFrame(contr)
R3a, C3a = readout_sep3(t3c, ['NKT_c05-IFNG-CCL4hi'], 'NKT_c05')
R3b, C3b = readout_sep3(t3c2, ['NKT_c04-IFNG', 'NKT_c05-IFNG-CCL4hi'], 'NKT_c04+c05')
R3 = pd.concat([R3a, R3b]); C3 = pd.concat([C3a, C3b]); R3.to_csv(f'{OUT}/sep3_fractions.csv', index=False); C3.to_csv(f'{OUT}/sep3_contrast.csv', index=False)
print(R3.to_string(index=False)); print(C3.to_string(index=False))
TT = pd.DataFrame([T4, T3p, T3c, T3c2]); TT.to_csv(f'{OUT}/classifier_performance.csv', index=False); print(TT.to_string(index=False))""")

code("""fig, axes = plt.subplots(1, 3, figsize=(20, 6.5))
def panel(ax, R, title, srcs):
    xs = np.arange(len(srcs)); w = 0.36
    for j, (G, col) in enumerate([('HC', '#7f7f7f'), ('sepsis', '#D7191C')]):
        d = R[R.group == G].set_index('source').reindex(srcs)
        ax.bar(xs + (j - .5) * w, d.est * 100, w, color=col, label=G, yerr=[(d.est - d.lo) * 100, (d.hi - d.est) * 100], capsize=3)
    ax.set_xticks(xs); ax.set_xticklabels([s.replace(', ', '\\n') for s in srcs], fontsize=9); ax.set_ylabel('% of parent'); ax.set_title(title); ax.legend()
srcs = ['real, true labels', 'realT2, classifier', 'gen k1000, classifier', 'gen k4000, classifier']
panel(axes[0], R4, f'SEP-4A Mono_c01 within CD14 monocytes\\n(classifier AUC on held-out donors {T4[\"test_auc_B\"]:.2f})', srcs)
panel(axes[1], R3[R3.target == 'NKT_c05'], f'SEP-3 NKT_c05 within NKT\\n(parent AUC {T3p[\"test_auc_B\"]:.2f}, child AUC {T3c[\"test_auc_B\"]:.2f})', srcs)
ax = axes[2]; d = MH.set_index('source'); ax.bar(range(len(d)), d.d, color=['#7f7f7f', '#E3A21A', '#2B83BA', '#D7191C'])
ax.set_xticks(range(len(d))); ax.set_xticklabels(d.index, fontsize=9); ax.axhline(0, c='k', lw=.5); ax.set_ylabel("Cohen's d, sepsis vs HC (per-sample)")
ax.set_title('MHC class II score in CD14 monocytes (threshold-free)')
fig.tight_layout(); fig.savefig(f'{OUT}/claim_classifiers.png', dpi=130); plt.show()""")

code("""# per-cell predictions as h5ad for inspection (generated + realT2 B), SEP-4A and SEP-3
rows = []
for src, X, obs, grp in [('realT2 B', r.layers['t2'], r.obs, grp_r), ('gen', g.X, g.obs, grp_g)]:
    keep = (obs.half == 'B').to_numpy() if src == 'realT2 B' else np.ones(len(obs), bool)
    ix = np.where(keep)[0]; o = obs.iloc[ix]
    lab = o.R3_fine.astype(str); lin = o.R3_lin6.astype(str)
    p4 = np.full(len(ix), np.nan); m = cd14_true_like(lab); p4[m] = t4.predict(X[ix[m]])
    pp = np.full(len(ix), np.nan); m2 = lin.isin(['T', 'NK']).to_numpy(); pp[m2] = t3p.predict(X[ix[m2]])
    pc = np.full(len(ix), np.nan); m3 = m2.copy(); m3[m2] = pp[m2] >= t3p.thr; pc[m3] = t3c.predict(X[ix[m3]])
    a = ad.AnnData(sp.csr_matrix(X[ix]), obs=pd.DataFrame({'source': src if src != 'gen' else o.K.astype(str).values, 'sample': o['sample'].astype(str).values,
          'group': grp[ix], 'R3_fine': lab.values, 'p_sep4a_monoc01': p4, 'p_nkt_parent': pp, 'p_nktc05_given_nkt': pc,
          'true_fine': o.fine.astype(str).values if src == 'realT2 B' else 'NA'}), var=pd.DataFrame(index=r.var_names))
    a.obs_names = [f'{src}|{i}' for i in range(a.n_obs)]; rows.append(a)
A = ad.concat(rows); A.uns['thresholds'] = json.dumps({'sep4a': float(t4.thr), 'nkt_parent': float(t3p.thr), 'nktc05': float(t3c.thr)})
A.obs['source'] = A.obs.source.astype('category'); A.obs['group'] = A.obs.group.astype('category')
A.write_h5ad(f'{OUT}/claim_classifier_predictions.h5ad', compression='gzip'); A""")

nb = nbf.v4.new_notebook(); nb.cells = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
NotebookClient(nb, timeout=10800, kernel_name='python3').execute()
nbf.write(nb, f'{OUT}/claim_classifiers.ipynb'); print('done')
