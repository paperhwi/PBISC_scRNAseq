import nbformat as nbf
C = []
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Sepsis real patient vs real patient (1:1) iLISI — baseline for integration report p.18 — 2026-10-18

Report p.18 shows median iLISI of **real vs generated** cells (sepsis K1k 1.22, K4k 1.22; COVID 1.16 / 1.17 after Harmony).
Here the same pipeline is applied to **one real patient vs another real patient**, so the generated gap can be read against
the patient-to-patient difference that exists in real data.

Pipeline (same as `code/c_joint.py` of the report, per pair):
1. Ye 2026 sepsis recap, cells in the paper scheme; one scRNA-seq sample per patient (133 samples = 133 patients: adult sepsis 93,
   HC 7, ICU control 8; paediatric sepsis 17, PHC 4, PIC 4). 800 random cells per patient extracted on the PC
   (`code/extract_sepsis_patients.py`), on the 5,853 generated-panel genes present in the real data.
2. Each repeat: draw 400 of the 800 cells per patient; random perfect matching of patients into pairs (66 pairs; one patient left out).
3. Per pair (patient p = reference, q = query): log1p CP10k on panel totals; 2,000 HVG (seurat_v3) on p; scale with p's mean/SD; PCA 30 fitted on p,
   q projected (**PCA projection**); Harmony (harmonypy, batch = patient, max_iter 20, seed 0) (**after Harmony**).
4. iLISI (perplexity 30, k = 90) on all 800 cells; median per pair; **mean over pairs = the repeat value**. 5 repeats.

Controls: same-patient split (400 vs 400 cells of one patient; ceiling for this cell number) and adult-only pairing.""")
code("""import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad, harmonypy, json, os, time, warnings
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore'); os.makedirs('results', exist_ok=True); os.makedirs('figs', exist_ok=True); os.makedirs('data', exist_ok=True)
z = np.load('_extract/sepsis_patients_cells.npz', allow_pickle=True)
X = sp.csr_matrix((z['data'], z['indices'], z['indptr']), shape=tuple(z['shape']))
obs = pd.DataFrame({k: z[k] for k in ['sample','group','age_group','celltype_major_ye2026']})
genes = z['genes']; print(X.shape, obs['sample'].nunique(), 'patients')
pat = obs.groupby('sample').agg(group=('group','first'), age=('age_group','first'), n=('group','size'))
display(pat.groupby(['age','group']).size().rename('patients'))""")
code("""# metric: identical to code/metrics.py of the integration report (LISI, Korsunsky 2019; perplexity 30, k = 90)
def _lisi_from_knn(d, idx, labels, perplexity=30, tol=1e-5, iters=50):
    lab = pd.Categorical(labels).codes; L = lab.max() + 1
    n = len(d); out = np.zeros(n); logU = np.log(perplexity)
    for i in range(n):
        di = d[i] ** 2; beta = 1.0; lo, hi = -np.inf, np.inf
        for _ in range(iters):
            P = np.exp(-di * beta); s = P.sum()
            if s == 0: H = 0.0; P = np.ones_like(di) / len(di)
            else: H = np.log(s) + beta * (di * P).sum() / s; P = P / s
            diff = H - logU
            if abs(diff) < tol: break
            if diff > 0: lo = beta; beta = beta * 2 if hi == np.inf else (beta + hi) / 2
            else: hi = beta; beta = beta / 2 if lo == -np.inf else (beta + lo) / 2
        w = np.bincount(lab[idx[i]], weights=P, minlength=L); out[i] = 1.0 / (w ** 2).sum()
    return out
def ilisi(E, b, k=90):
    d, idx = NearestNeighbors(n_neighbors=k).fit(E).kneighbors(E)
    return float(np.median(_lisi_from_knn(d, idx, b)))
def lognorm(M):
    tot = np.asarray(M.sum(1)).ravel(); Y = (sp.diags(1e4 / np.maximum(tot, 1)) @ M).tocsr(); Y.data = np.log1p(Y.data); return Y
def pair_ilisi(ia, ib, seed=0, XB=None):
    Xa = X[ia]; Xb = X[ib] if XB is None else XB[ib]
    a = ad.AnnData(Xa.copy())
    try: sc.pp.highly_variable_genes(a, flavor='seurat_v3', n_top_genes=2000)
    except Exception: sc.pp.highly_variable_genes(a, flavor='seurat_v3', n_top_genes=2000, span=0.5)
    hv = np.where(a.var.highly_variable.to_numpy())[0]
    Ya, Yb = lognorm(Xa), lognorm(Xb)
    mu = np.asarray(Ya[:, hv].mean(0)).ravel(); sd = np.sqrt(np.maximum(np.asarray(Ya[:, hv].multiply(Ya[:, hv]).mean(0)).ravel() - mu ** 2, 0)) + 1e-6
    scale = lambda Y: np.clip((Y[:, hv].toarray() - mu) / sd, -10, 10)
    pca = PCA(30, random_state=0).fit(scale(Ya))
    pre = np.vstack([pca.transform(scale(Ya)), pca.transform(scale(Yb))])
    b = np.array(['p'] * len(ia) + ['q'] * len(ib))
    ho = harmonypy.run_harmony(pre.astype(np.float64), pd.DataFrame({'b': b}), ['b'], random_state=seed, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); post = Z.T if Z.shape[0] == pre.shape[1] else Z
    return ilisi(pre, b), ilisi(post, b)
cells = {s: np.where(obs['sample'].to_numpy() == s)[0] for s in pat.index}
t = time.time(); print('test pair', pair_ilisi(cells[pat.index[0]][:400], cells[pat.index[1]][:400]), round(time.time() - t, 1), 's')""")
code("""def run_repeat(r, patients, tag):
    rng = np.random.default_rng(1000 + r)
    sub = {s: rng.choice(cells[s], 400, replace=False) for s in patients}
    perm = rng.permutation(patients); pairs = [(perm[i], perm[i + 1]) for i in range(0, len(perm) - 1, 2)]
    res = Parallel(n_jobs=2)(delayed(pair_ilisi)(sub[p], sub[q]) for p, q in pairs)
    return pd.DataFrame([{'set': tag, 'repeat': r + 1, 'patient_p': p, 'patient_q': q,
                          'group_p': pat.loc[p, 'group'], 'group_q': pat.loc[q, 'group'], 'age_p': pat.loc[p, 'age'], 'age_q': pat.loc[q, 'age'],
                          'iLISI_pre': a, 'iLISI_post': b} for (p, q), (a, b) in zip(pairs, res)])
t0 = time.time(); allp = list(pat.index); adults = list(pat.index[pat.age == 'adult'])
R = []
for r in range(5):
    R.append(run_repeat(r, allp, 'all patients')); print('all, repeat', r + 1, round(time.time() - t0), 's', flush=True)
for r in range(5):
    R.append(run_repeat(r, adults, 'adults only')); print('adults, repeat', r + 1, round(time.time() - t0), 's', flush=True)
# same-patient split control (one per patient, repeat 1 sampling)
def split_one(s, seed):
    rng = np.random.default_rng(seed); ix = rng.permutation(cells[s]); return pair_ilisi(ix[:400], ix[400:800])
sp_res = Parallel(n_jobs=2)(delayed(split_one)(s, 5000 + i) for i, s in enumerate(allp))
split = pd.DataFrame({'set': 'same patient split', 'repeat': 1, 'patient_p': allp, 'patient_q': allp, 'group_p': pat.group.values, 'group_q': pat.group.values,
                      'age_p': pat.age.values, 'age_q': pat.age.values, 'iLISI_pre': [a for a, _ in sp_res], 'iLISI_post': [b for _, b in sp_res]})
P = pd.concat(R + [split], ignore_index=True); P.to_csv('results/pair_ilisi_all_pairs.csv', index=False); print('total', round(time.time() - t0), 's')""")
md("## Result: mean iLISI over patient pairs, per repeat")
code("""per_rep = P.groupby(['set', 'repeat'])[['iLISI_pre', 'iLISI_post']].mean().round(3)
per_rep['n pairs'] = P.groupby(['set', 'repeat']).size()
display(per_rep); per_rep.to_csv('results/pair_ilisi_per_repeat.csv')
summ = per_rep.groupby('set')[['iLISI_pre', 'iLISI_post']].agg(['mean', 'std', 'min', 'max']).round(3); display(summ); summ.to_csv('results/pair_ilisi_summary.csv')
med = P.groupby('set')[['iLISI_pre', 'iLISI_post']].median().round(3); print('median over all pairs:'); display(med)""")
code("""# context: by pair type (all-patients pairing), and the report's pooled values
def kind(r):
    g = sorted([r.group_p, r.group_q]); return ' vs '.join(g)
A_ = P[P.set == 'all patients'].copy(); A_['pair type'] = A_.apply(kind, axis=1)
display(A_.groupby('pair type')[['iLISI_pre', 'iLISI_post']].agg(['mean', 'count']).round(3))
J = json.load(open('../results/sepsis_joint.json'))['mix']; Jc = json.load(open('../results/covid_joint.json'))['mix']
cal = pd.read_csv('../results/check/calibration_summary.csv')
ref = pd.DataFrame({'pre': [J['iLISI_real_vs_k1000_pre'], J['iLISI_real_vs_k4000_pre'], Jc['iLISI_real_vs_k1000_pre'], Jc['iLISI_real_vs_k4000_pre']],
                    'post': [J['iLISI_real_vs_k1000_post'], J['iLISI_real_vs_k4000_post'], Jc['iLISI_real_vs_k1000_post'], Jc['iLISI_real_vs_k4000_post']]},
                   index=['Sepsis real vs K=1,000', 'Sepsis real vs K=4,000', 'COVID real vs K=1,000', 'COVID real vs K=4,000']).round(3)
display(ref); display(cal[['cohort', 'source', 'iLISI_pre', 'iLISI_post']].dropna().round(3))""")
md("""## Same pipeline, real patient vs one generated sample (1:1) — like-for-like reference
The report's p.18 values pool all real and all generated cells; the patient-pair values use 400 + 400 cells. To compare at the same scale,
each real patient (reference) is paired with one PBISC bulk-only sepsis sample (400 cells, the sample's five draws pooled), K = 1,000 and K = 4,000, 5 repeats.""")
code("""GEN = {}
for K in ['k1000', 'k4000']:
    g = np.load(f'_extract/sepsis_{K}_samples_cells.npz', allow_pickle=True); assert (g['genes'] == genes).all()
    GEN[K] = (sp.csr_matrix((g['data'], g['indices'], g['indptr']), shape=tuple(g['shape'])), pd.DataFrame({'sample': g['sample'], 'draw': g['draw'], 'group': g['group']}))
def run_gen(r, K):
    XG, og = GEN[K]; rng = np.random.default_rng(2000 + r)
    gs = rng.permutation(og['sample'].unique()); rp = rng.choice(allp, len(gs), replace=False)
    gcells = {s: np.where(og['sample'].to_numpy() == s)[0] for s in gs}
    jobs = [(rng.choice(cells[p], 400, replace=False), rng.choice(gcells[q], 400, replace=False), p, q) for p, q in zip(rp, gs)]
    res = Parallel(n_jobs=2)(delayed(pair_ilisi)(a, b, 0, XG) for a, b, _, _ in jobs)
    return pd.DataFrame([{'set': f'real patient vs {K} sample', 'repeat': r + 1, 'patient_p': p, 'patient_q': q, 'group_p': pat.loc[p, 'group'],
                          'group_q': og.loc[og['sample'] == q, 'group'].iloc[0], 'iLISI_pre': a, 'iLISI_post': b} for (_, _, p, q), (a, b) in zip(jobs, res)])
t0 = time.time(); G_ = []
for K in ['k1000', 'k4000']:
    for r in range(5): G_.append(run_gen(r, K)); print(K, 'repeat', r + 1, round(time.time() - t0), 's', flush=True)
PG = pd.concat(G_, ignore_index=True); PG.to_csv('results/pair_ilisi_real_vs_generated_pairs.csv', index=False)
per_rep_g = PG.groupby(['set', 'repeat'])[['iLISI_pre', 'iLISI_post']].mean().round(3); per_rep_g['n pairs'] = PG.groupby(['set', 'repeat']).size()
display(per_rep_g); per_rep_g.to_csv('results/pair_ilisi_real_vs_generated_per_repeat.csv')
both = pd.concat([per_rep, per_rep_g]); summ2 = both.groupby('set')[['iLISI_pre', 'iLISI_post']].agg(['mean', 'std']).round(3); display(summ2); summ2.to_csv('results/pair_ilisi_summary_all.csv')""")
md("## Figure for report p.18 (existing bars + real patient vs patient)")
code("""allr = per_rep.loc['all patients']; split_r = per_rep.loc['same patient split']
labels = ['Sepsis\\nreal vs K = 1,000', 'Sepsis\\nreal vs K = 4,000', 'COVID\\nreal vs K = 1,000', 'COVID\\nreal vs K = 4,000', 'Sepsis real patient\\nvs patient (1:1,\\n400+400 cells, 5 repeats)']
pre = list(ref.pre) + [allr.iLISI_pre.mean()]; post = list(ref.post) + [allr.iLISI_post.mean()]
err_pre = [None] * 4 + [allr.iLISI_pre.std()]; err_post = [None] * 4 + [allr.iLISI_post.std()]
cpost = ['#2A7AB8', '#D7191C', '#2A7AB8', '#D7191C', '#3C9D5D']
fig, ax = plt.subplots(figsize=(10.5, 5.2)); x = np.arange(5); w = 0.36
ax.bar(x - w/2, pre, w, color='#C8C8C8', label='PCA projection'); ax.bar(x + w/2, post, w, color=cpost, label='after Harmony')
ax.errorbar([4 - w/2, 4 + w/2], [pre[4], post[4]], yerr=[err_pre[4], err_post[4]], fmt='none', ecolor='k', capsize=4)
for i, v in enumerate(post): ax.text(i + w/2, v + (err_post[i] or 0) + 0.03, f'{v:.2f}', ha='center', fontsize=10)
ax.axhline(2, ls='--', c='k', lw=1); ax.text(4.45, 2.03, 'perfect mixing = 2', ha='right', fontsize=10)
ax.axhline(1, ls=':', c='k', lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10); ax.set_ylim(0, 2.6); ax.set_ylabel('median iLISI (two sources)')
ax.spines[['top', 'right']].set_visible(False)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#C8C8C8', label='PCA projection'), Patch(color='#555555', label='after Harmony')], frameon=False, loc='upper left', ncol=2)
ax.set_title('Report p.18 (pooled cells) + sepsis real patient vs patient (1:1 pairs, 400 + 400 cells; different scale)', fontsize=10)
plt.tight_layout(); plt.savefig('figs/p18_with_patient_pairs.png', dpi=200); plt.savefig('figs/p18_with_patient_pairs.pdf'); plt.show()
print('same-patient split ceiling (400 vs 400 cells): pre', round(split_r.iLISI_pre.mean(), 3), 'post', round(split_r.iLISI_post.mean(), 3))""")
code("""# Figure B: all 1:1 pairs at the same scale (400 + 400 cells), mean of 5 repeats +/- SD
sets = [('all patients', 'real patient vs\\nreal patient', '#3C9D5D'), ('real patient vs k1000 sample', 'real patient vs\\nK = 1,000 sample', '#2A7AB8'),
        ('real patient vs k4000 sample', 'real patient vs\\nK = 4,000 sample', '#D7191C'), ('same patient split', 'same patient\\n(400 vs 400 split)', '#8C8C8C')]
fig, ax = plt.subplots(figsize=(8.8, 5.0)); x = np.arange(len(sets)); w = 0.36
for i, (k, lab, c) in enumerate(sets):
    r = both.loc[k]; a, b = r.iLISI_pre.mean(), r.iLISI_post.mean(); ea, eb = (r.iLISI_pre.std(), r.iLISI_post.std()) if len(r) > 1 else (0, 0)
    ax.bar(i - w/2, a, w, color='#C8C8C8'); ax.bar(i + w/2, b, w, color=c)
    if len(r) > 1: ax.errorbar([i - w/2, i + w/2], [a, b], yerr=[ea, eb], fmt='none', ecolor='k', capsize=4)
    ax.text(i + w/2, b + 0.05, f'{b:.2f}', ha='center', fontsize=10)
ax.axhline(2, ls='--', c='k', lw=1); ax.axhline(1, ls=':', c='k', lw=0.8); ax.set_ylim(0, 2.3)
ax.set_xticks(x); ax.set_xticklabels([l for _, l, _ in sets], fontsize=10); ax.set_ylabel('median iLISI per pair (mean of pairs)')
ax.spines[['top', 'right']].set_visible(False); from matplotlib.patches import Patch
ax.legend(handles=[Patch(color='#C8C8C8', label='PCA projection'), Patch(color='#555555', label='after Harmony')], frameon=False, loc='upper left', ncol=2)
ax.set_title('Sepsis, 1:1 pairs at the same scale (5 repeats)', fontsize=11)
plt.tight_layout(); plt.savefig('figs/pairs_1to1_same_scale.png', dpi=200); plt.savefig('figs/pairs_1to1_same_scale.pdf'); plt.show()""")
code("""A = ad.AnnData(X, obs=obs.astype(str)); A.var_names = pd.Index(genes).astype(str); A.var_names_make_unique()
A.uns['note'] = 'Ye 2026 sepsis recap, 800 random in-paper-scheme cells per sample (133 samples/patients), raw counts on 5,853 PBISC panel genes.'
A.write_h5ad('data/sepsis_patient_cells_800.h5ad', compression='gzip'); print(A)""")
nb = nbf.v4.new_notebook(); nb['cells'] = C; nb['metadata'] = {'kernelspec': {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}}
nbf.write(nb, '/home/claude/ilisi/patient_pair_ilisi.ipynb'); print('ok')
