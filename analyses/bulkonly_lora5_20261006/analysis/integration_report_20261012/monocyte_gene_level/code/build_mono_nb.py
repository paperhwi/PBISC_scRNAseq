import nbformat as nbf
from nbclient import NotebookClient
OUT = "/mnt/user-data/outputs/integration_report_20261012/monocyte_gene_level"
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s))
code = lambda s: cells.append(nbf.v4.new_code_cell(s))

md("""# Monocyte residual at gene level: generated vs real CD14 monocytes

Goal: which genes and programmes differ between generated and real CD14 monocytes, beyond the read-out format and donor variation, and do generated monocytes carry the disease programme?

**Sources** (same cells and read-out draws as `prop/c_mapping.py`, seed 7):
- `gen K1000`, `gen K4000`: generated cells (100 per sample per draw); CD14 monocytes by the format-aware classifier (route R3, fine).
- `real B`: held-out real cells (half B of the samples), true labels, untouched.
- `realT2 B`: the same cells through the PBISC-like read-out (15-NN smoothed profile on the cell's own genes, generated-lineage depth, one Poisson draw).
- `realT2 A` / `real A`: the other half of the real samples. **Control** = A vs B: same study and format, different donors.

Adults only (the generated set has a single paediatric sample). Case group: adult sepsis / COVID positive; control group: adult healthy / COVID negative.

**Why pseudobulk first.** The read-out emulation is not a perfect format match: realT2 cells detect fewer genes than generated cells (sepsis monocytes: 244 vs 316 genes). Per-cell statistics (AUC, variance, classifiers) are therefore format-sensitive. The per-sample mean CP10k (pseudobulk) is insensitive to the count draw. So gene-level calls use per-sample pseudobulk, tested across samples, and must hold against both `real B` and `realT2 B`.""")

code("""import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, os, warnings
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
M = '/home/claude/work/mono'; OUT = '%s'; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(3)
pd.set_option('display.width', 230); pd.set_option('display.float_format', '{:.3f}'.format)""" % OUT)

code("""def cd14(lab):
    lab = pd.Series(lab).astype(str)
    return (lab.str.match(r'Mono_c0[1-6]') | lab.isin(['CD14mono', 'CD14mono_IFN', 'CD14mono_dex'])).to_numpy()

def cp10k(X):
    X = sp.csr_matrix(X, dtype=np.float32); return (sp.diags(1e4 / np.maximum(X.sum(1).A1, 1)) @ X).tocsr()

COH = {'sepsis': dict(case_r='sepsis', ctrl_r='healthy control', case_g='adult_sepsis', ctrl_g='adult_HC'),
       'covid': dict(case_r='covid_positive', ctrl_r='covid_negative', case_g='covid_positive', ctrl_g='covid_negative')}
D = {}
for c, k in COH.items():
    r = ad.read_h5ad(f'{M}/{c}_real_AB.h5ad'); g = ad.read_h5ad(f'{M}/{c}_gen.h5ad')
    age = pd.read_csv(f'{M}/{c}_real_AB_obs_age.csv.gz').age.astype(str).to_numpy() if c == 'sepsis' else np.full(r.n_obs, 'adult')
    mr = cd14(r.obs.fine) & (age == 'adult'); mg = cd14(g.obs.R3_fine)
    cellsets = []
    for half in ['A', 'B']:
        ix = np.where(mr & (r.obs.half == half).to_numpy())[0]
        for lay, name in [(None, f'real {half}'), ('t2', f'realT2 {half}')]:
            X = r.X[ix] if lay is None else r.layers[lay][ix]
            cellsets.append((name, X, r.obs['sample'].to_numpy()[ix], r.obs.group.to_numpy()[ix]))
    for K in ['k1000', 'k4000']:
        ix = np.where(mg & (g.obs.K == K).to_numpy())[0]
        grp = g.obs.group.to_numpy()[ix]
        grp = np.where(grp == k['case_g'], k['case_r'], np.where(grp == k['ctrl_g'], k['ctrl_r'], grp))
        cellsets.append(({'k1000': 'gen K1,000', 'k4000': 'gen K4,000'}[K], g.X[ix], g.obs['sample'].to_numpy()[ix], grp))
    D[c] = dict(genes=np.array(r.var_names), sets={n: dict(X=sp.csr_matrix(X, dtype=np.float32), sample=s.astype(str), group=gr.astype(str)) for n, X, s, gr in cellsets})
names = list(D['sepsis']['sets']); print(names)
rows = []
for c in D:
    for n, s in D[c]['sets'].items():
        for gname in [COH[c]['case_r'], COH[c]['ctrl_r']]:
            m = s['group'] == gname; X = s['X'][m]
            rows.append(dict(cohort=c, source=n, group=gname, cells=int(m.sum()), samples=len(set(s['sample'][m])),
                             UMI_median=np.median(X.sum(1).A1) if m.any() else np.nan, genes_detected_median=np.median(np.diff(X.indptr)) if m.any() else np.nan))
pd.DataFrame(rows)""")

md("""## 1. Gene-level differences in the case group (per-sample pseudobulk)

Per sample (≥ 20 CD14 monocytes): mean CP10k per gene. Test generated samples vs real samples (Mann–Whitney on log2(CP10k + 0.1), BH q). Control: real A vs real B (donor variation within one study).

**Model-specific gene** (expressed: mean CP10k ≥ 1 in either source): |log2FC| ≥ 1 and q < 0.01 against **both** real B and realT2 B, same sign, and |log2FC| < 0.5 in the A vs B control.""")

code("""def pseudobulk(s, group, minc=20):
    m = s['group'] == group; X = cp10k(s['X'][m]); smp = s['sample'][m]
    out = {}
    for u in np.unique(smp):
        w = smp == u
        if w.sum() >= minc: out[u] = np.asarray(X[w].mean(0)).ravel()
    return pd.DataFrame(out, index=None).T

def de(PA, PB, genes):
    la, lb = np.log2(PA.to_numpy() + 0.1), np.log2(PB.to_numpy() + 0.1)
    p = np.array([mannwhitneyu(la[:, j], lb[:, j]).pvalue if (la[:, j].std() + lb[:, j].std()) > 0 else 1.0 for j in range(la.shape[1])])
    p = np.nan_to_num(p, nan=1.0); q = multipletests(p, method='fdr_bh')[1]
    return pd.DataFrame({'gene': genes, 'log2FC': np.log2((PA.mean(0).to_numpy() + 0.1) / (PB.mean(0).to_numpy() + 0.1)), 'q': q,
                         'cp10k_a': PA.mean(0).to_numpy(), 'cp10k_b': PB.mean(0).to_numpy()})

PB = {}; GT = {}
for c in D:
    case = COH[c]['case_r']; genes = D[c]['genes']
    PB[c] = {n: pseudobulk(s, case) for n, s in D[c]['sets'].items()}
    print(c, 'samples with >=20 CD14 monocytes:', {n: len(v) for n, v in PB[c].items()})
    T = pd.DataFrame({'gene': genes})
    for a, b, tag in [('gen K1,000', 'real B', 'gen_vs_realB'), ('gen K1,000', 'realT2 B', 'gen_vs_realT2B'), ('gen K4,000', 'real B', 'gen4k_vs_realB'),
                      ('real A', 'real B', 'ctrl_realA_vs_realB'), ('realT2 B', 'real B', 'format_realT2B_vs_realB')]:
        d = de(PB[c][a], PB[c][b], genes); T[f'log2FC|{tag}'] = d.log2FC.values; T[f'q|{tag}'] = d.q.values
    T['cp10k|gen K1,000'] = PB[c]['gen K1,000'].mean(0).values; T['cp10k|real B'] = PB[c]['real B'].mean(0).values; T['cp10k|realT2 B'] = PB[c]['realT2 B'].mean(0).values
    T['expressed'] = T[['cp10k|gen K1,000', 'cp10k|real B']].max(1) >= 1
    x1, x2, xc = T['log2FC|gen_vs_realB'], T['log2FC|gen_vs_realT2B'], T['log2FC|ctrl_realA_vs_realB']
    T['model_specific'] = T.expressed & (x1.abs() >= 1) & (x2.abs() >= 1) & (np.sign(x1) == np.sign(x2)) & (T['q|gen_vs_realB'] < .01) & (T['q|gen_vs_realT2B'] < .01) & (xc.abs() < .5)
    T['direction'] = np.where(x1 > 0, 'higher in generated', 'lower in generated')
    GT[c] = T; T.to_csv(f'{OUT}/{c}_cd14mono_pseudobulk_de.csv.gz', index=False)
    e = T[T.expressed]
    print(f'   expressed {len(e)} | |log2FC|>=1 & q<.01 vs real B: {int(((e["log2FC|gen_vs_realB"].abs()>=1)&(e["q|gen_vs_realB"]<.01)).sum())}; '
          f'control A vs B: {int(((e["log2FC|ctrl_realA_vs_realB"].abs()>=1)&(e["q|ctrl_realA_vs_realB"]<.01)).sum())}; '
          f'format realT2 vs real: {int(((e["log2FC|format_realT2B_vs_realB"].abs()>=1)&(e["q|format_realT2B_vs_realB"]<.01)).sum())}; model-specific: {int(e.model_specific.sum())}')
    print('   K4000 vs K1000 agreement (Spearman of log2FC vs real B):', round(spearmanr(e['log2FC|gen_vs_realB'], e['log2FC|gen4k_vs_realB']).correlation, 3))
    cols = ['gene', 'cp10k|real B', 'cp10k|realT2 B', 'cp10k|gen K1,000', 'log2FC|gen_vs_realB', 'log2FC|gen_vs_realT2B', 'log2FC|ctrl_realA_vs_realB']
    for dname, asc in [('lower in generated', True), ('higher in generated', False)]:
        x = e[e.model_specific & (e.direction == dname)].sort_values('log2FC|gen_vs_realB', ascending=asc)
        print(f'--- {dname}: {len(x)} genes (top 30)'); print(x[cols].head(30).to_string(index=False))""")

md("""## 2. Gene programmes (per-sample pseudobulk)

Programme score per sample = mean log2(CP10k + 0.1) over programme genes in the panel. Effect = generated minus real B, in units of the real between-sample SD. Control = real A minus real B.""")

code("""P = {
 'MHC class II': ['HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DMA', 'HLA-DMB', 'CD74'],
 'Alarmin S100 / emergency myelopoiesis': ['S100A8', 'S100A9', 'S100A12', 'LCN2', 'MMP9', 'RETN', 'CYP1B1', 'MCEMP1'],
 'Inflammatory / NF-kB': ['IL1B', 'CXCL8', 'TNF', 'NFKBIA', 'CCL3', 'CCL4', 'CCL3L1', 'IER3', 'PTGS2', 'G0S2', 'SOD2'],
 'Interferon': ['ISG15', 'IFI6', 'IFI44L', 'IFIT1', 'IFIT2', 'IFIT3', 'MX1', 'RSAD2', 'OAS1', 'IFITM3', 'LY6E', 'XAF1', 'IFI27'],
 'Immediate early': ['FOS', 'FOSB', 'JUN', 'JUNB', 'EGR1', 'DUSP1', 'IER2', 'KLF6', 'ZFP36'],
 'Ribosomal': 'RP', 'Mitochondrial': 'MT',
 'Classical identity': ['CD14', 'LYZ', 'VCAN', 'FCN1', 'CSF1R', 'CST3', 'MNDA', 'S100A6', 'CSF3R', 'SPI1'],
 'Lysosome / antigen processing': ['CTSS', 'CTSB', 'CTSD', 'LGALS3', 'PSAP', 'GRN', 'CD68'],
 'Glucocorticoid response': ['TSC22D3', 'FKBP5', 'ZBTB16', 'CD163', 'VSIG4', 'ADAMTS2', 'MERTK'],
 'Histone (replication-dependent)': 'HIST',
 'Lymphoid genes (leakage)': ['CD3E', 'CD3D', 'CD2', 'IL7R', 'LTB', 'CD247', 'GNLY', 'NKG7', 'MS4A1', 'CD79A', 'IGHM']}
def pidx(genes, gl):
    g = pd.Index(genes)
    if gl == 'RP': return np.where(g.str.match(r'^RP[LS]\\d'))[0]
    if gl == 'MT': return np.where(g.str.startswith('MT-'))[0]
    if gl == 'HIST': return np.where(g.str.startswith('HIST'))[0]
    return g.get_indexer([x for x in gl if x in g])
rows = []
for c in D:
    genes = D[c]['genes']
    for p, gl in P.items():
        gi = pidx(genes, gl)
        sc_ = {n: np.log2(v.to_numpy()[:, gi] + 0.1).mean(1) for n, v in PB[c].items()}
        sd = sc_['real B'].std(ddof=1)
        for n in ['gen K1,000', 'gen K4,000', 'realT2 B', 'real A']:
            rows.append(dict(cohort=c, programme=p, n_genes=len(gi), source=n, effect_vs_realB_in_SD=(sc_[n].mean() - sc_['real B'].mean()) / sd,
                             p=mannwhitneyu(sc_[n], sc_['real B']).pvalue))
PR = pd.DataFrame(rows); PR.to_csv(f'{OUT}/programme_scores_pseudobulk.csv', index=False)
PR.pivot_table(index=['cohort', 'programme'], columns='source', values='effect_vs_realB_in_SD')[['gen K1,000', 'gen K4,000', 'realT2 B', 'real A']]""")

md("""## 3. Lineage leakage at the cell level

Fraction of CD14 monocytes with any count for T / NK / B genes. Generated cells detect somewhat more genes than realT2 cells, so the fair reference is between realT2 B (sparser) and real B (denser).""")

code("""LEAK = {'T (CD3E/CD3D/CD2/IL7R)': ['CD3E', 'CD3D', 'CD2', 'IL7R'], 'NK/CTL (GNLY/NKG7/PRF1/GZMB)': ['GNLY', 'NKG7', 'PRF1', 'GZMB'],
        'B (MS4A1/CD79A/CD79B)': ['MS4A1', 'CD79A', 'CD79B'], 'Platelet (PPBP/PF4)': ['PPBP', 'PF4']}
rows = []
for c in D:
    g = pd.Index(D[c]['genes'])
    for n, s in D[c]['sets'].items():
        m = s['group'] == COH[c]['case_r']; X = s['X'][m]
        for L, gl in LEAK.items():
            gi = g.get_indexer([x for x in gl if x in g]); rows.append(dict(cohort=c, source=n, lineage=L, frac_cells_detected=float((X[:, gi].sum(1).A1 > 0).mean())))
LK = pd.DataFrame(rows); LK.to_csv(f'{OUT}/lineage_leakage.csv', index=False)
LK.pivot_table(index=['cohort', 'lineage'], columns='source', values='frac_cells_detected')[['real B', 'realT2 B', 'gen K1,000', 'gen K4,000']]""")

md("""## 4. Do generated monocytes carry the disease programme?

Within CD14 monocytes, per-gene disease log2FC (case vs control, per-sample pseudobulk) in real (A + B together, untouched) and in generated cells. Agreement over expressed genes and over the real disease genes (|log2FC| ≥ 1, q < 0.05).""")

code("""DIS = {}
for c in D:
    case, ctrl = COH[c]['case_r'], COH[c]['ctrl_r']; genes = D[c]['genes']
    def pb_both(names, grp):
        return pd.concat([pseudobulk(D[c]['sets'][n], grp) for n in names])
    out = {}
    for tag, names in [('real', ['real A', 'real B']), ('gen K1,000', ['gen K1,000']), ('gen K4,000', ['gen K4,000'])]:
        out[tag] = de(pb_both(names, case), pb_both(names, ctrl), genes)
    T = pd.DataFrame({'gene': genes, **{f'log2FC|{t}': v.log2FC.values for t, v in out.items()}, **{f'q|{t}': v.q.values for t, v in out.items()}})
    T['expressed'] = GT[c].expressed.values; DIS[c] = T; T.to_csv(f'{OUT}/{c}_cd14mono_disease_logFC.csv.gz', index=False)
    e = T[T.expressed]; sig = e[(e['log2FC|real'].abs() >= 1) & (e['q|real'] < .05)]
    for t in ['gen K1,000', 'gen K4,000']:
        sl = np.polyfit(sig['log2FC|real'], sig[f'log2FC|{t}'], 1)[0] if len(sig) > 2 else np.nan
        print(f'{c} {t}: Spearman over {len(e)} expressed genes = {spearmanr(e["log2FC|real"], e[f"log2FC|{t}"]).correlation:.2f}; '
              f'real disease genes = {len(sig)}, sign agreement = {(np.sign(sig["log2FC|real"]) == np.sign(sig[f"log2FC|{t}"])).mean():.2f}, slope = {sl:.2f}')
    show = ['S100A8', 'S100A9', 'S100A12', 'HLA-DRA', 'HLA-DPB1', 'CD74', 'IL1B', 'CXCL8', 'TNF', 'ISG15', 'IFI27', 'IFITM3', 'MCEMP1', 'CD163', 'RETN', 'FKBP5', 'IL1R2', 'VCAN']
    print(e[e.gene.isin(show)][['gene', 'log2FC|real', 'q|real', 'log2FC|gen K1,000', 'log2FC|gen K4,000']].to_string(index=False))""")

md("""## 5. Cell-level separability and co-expression (format-sensitive; read with the controls)

5-fold logistic regression (L2, C = 0.1) on scaled log1p CP10k, 4,000 cells per class, case group. Co-expression: correlation among the 300 most expressed genes, similarity to realT2 B.""")

code("""def ln(X): Y = cp10k(X); Y.data = np.log1p(Y.data); return Y
def clf_auc(Xa, Xb, n=4000, C=0.1):
    ia = rng.choice(Xa.shape[0], min(n, Xa.shape[0]), replace=False); ib = rng.choice(Xb.shape[0], min(n, Xb.shape[0]), replace=False)
    X = sp.vstack([ln(Xa[ia]), ln(Xb[ib])]).toarray(); y = np.r_[np.ones(len(ia)), np.zeros(len(ib))]
    keep = X.std(0) > 0; X = (X[:, keep] - X[:, keep].mean(0)) / X[:, keep].std(0)
    p = cross_val_predict(LogisticRegression(C=C, max_iter=500), X, y, cv=StratifiedKFold(5, shuffle=True, random_state=0), method='predict_proba')[:, 1]
    return roc_auc_score(y, p)
CL = []; CO = []
for c in D:
    case = COH[c]['case_r']; ex = GT[c].expressed.to_numpy()
    Xs = {n: s['X'][s['group'] == case][:, ex] for n, s in D[c]['sets'].items()}
    for a, b in [('gen K1,000', 'realT2 B'), ('gen K1,000', 'real B'), ('real B', 'realT2 B'), ('realT2 A', 'realT2 B'), ('real A', 'real B')]:
        CL.append(dict(cohort=c, task=f'{a} vs {b}', cv_auc=clf_auc(Xs[a], Xs[b])))
    Y = {n: ln(v) for n, v in Xs.items()}; mu = np.asarray(Y['realT2 B'].mean(0)).ravel(); top = np.argsort(-mu)[:300]; iu = np.triu_indices(300, 1)
    C = {n: np.nan_to_num(np.corrcoef(v[rng.choice(v.shape[0], min(4000, v.shape[0]), replace=False)][:, top].toarray().T)) for n, v in Y.items()}
    for n in ['realT2 A', 'real B', 'gen K1,000', 'gen K4,000']:
        CO.append(dict(cohort=c, source=n, coexpr_similarity_to_realT2B=np.corrcoef(C[n][iu], C['realT2 B'][iu])[0, 1]))
CL = pd.DataFrame(CL); CO = pd.DataFrame(CO); CL.to_csv(f'{OUT}/classifier_auc.csv', index=False); CO.to_csv(f'{OUT}/coexpression_similarity.csv', index=False)
print(CL.to_string(index=False)); print(CO.to_string(index=False))""")

code("""fig, axes = plt.subplots(2, 3, figsize=(20, 12))
for i, c in enumerate(D):
    T = GT[c]; e = T[T.expressed]
    ax = axes[i, 0]
    ax.scatter(e['log2FC|ctrl_realA_vs_realB'], -np.log10(e['q|ctrl_realA_vs_realB'] + 1e-12), s=6, c='#bbbbbb', label='control: real A vs real B')
    ax.scatter(e['log2FC|gen_vs_realB'], -np.log10(e['q|gen_vs_realB'] + 1e-12), s=7, c='#2B83BA', alpha=.6, label='generated K1,000 vs real B')
    ms = e[e.model_specific]; ax.scatter(ms['log2FC|gen_vs_realB'], -np.log10(ms['q|gen_vs_realB'] + 1e-12), s=14, c='#D7191C', label=f'model-specific ({len(ms)})')
    lab = pd.concat([ms.nsmallest(9, 'log2FC|gen_vs_realB'), ms.nlargest(9, 'log2FC|gen_vs_realB')])
    for _, r in lab.iterrows(): ax.annotate(r.gene, (r['log2FC|gen_vs_realB'], -np.log10(r['q|gen_vs_realB'] + 1e-12)), fontsize=8)
    ax.set_xlabel('log2 fold change, per-sample pseudobulk (generated / real B)'); ax.set_ylabel('-log10 q'); ax.legend(fontsize=8)
    ax.set_title(f'{c}: CD14 monocytes, case group, expressed genes')
    ax = axes[i, 1]; d = PR[PR.cohort == c].pivot_table(index='programme', columns='source', values='effect_vs_realB_in_SD')[['gen K1,000', 'gen K4,000', 'realT2 B', 'real A']]
    d = d.loc[d['gen K1,000'].sort_values().index]; yy = np.arange(len(d))
    for j, (k, col) in enumerate([('gen K1,000', '#2B83BA'), ('gen K4,000', '#D7191C'), ('realT2 B', '#E3A21A'), ('real A', '#7f7f7f')]):
        ax.barh(yy + (j - 1.5) * 0.2, d[k], height=0.2, color=col, label=k)
    ax.set_yticks(yy); ax.set_yticklabels(d.index, fontsize=9); ax.axvline(0, c='k', lw=.5); ax.legend(fontsize=8)
    ax.set_xlabel('programme score vs real B (real between-sample SD units)'); ax.set_title('programmes (grey = donor control, orange = format control)')
    ax = axes[i, 2]; T2 = DIS[c]; e2 = T2[T2.expressed]
    ax.scatter(e2['log2FC|real'], e2['log2FC|gen K1,000'], s=6, c='#2B83BA', alpha=.5)
    sig = e2[(e2['log2FC|real'].abs() >= 1) & (e2['q|real'] < .05)]; ax.scatter(sig['log2FC|real'], sig['log2FC|gen K1,000'], s=12, c='#D7191C')
    for _, r in sig.reindex(sig['log2FC|real'].abs().sort_values().index[-14:]).iterrows(): ax.annotate(r.gene, (r['log2FC|real'], r['log2FC|gen K1,000']), fontsize=8)
    lim = np.nanmax(np.abs(e2[['log2FC|real', 'log2FC|gen K1,000']].to_numpy())) * 1.05; ax.plot([-lim, lim], [-lim, lim], 'k:', lw=.8); ax.axhline(0, c='k', lw=.4); ax.axvline(0, c='k', lw=.4)
    ax.set_xlabel(f'real disease log2FC ({COH[c]["case_r"]} / {COH[c]["ctrl_r"]})'); ax.set_ylabel('generated K1,000 disease log2FC')
    ax.set_title(f'disease programme within CD14 monocytes (red: real |log2FC|>=1, q<0.05, n={len(sig)})')
fig.tight_layout(); fig.savefig(f'{OUT}/monocyte_gene_level.png', dpi=130); plt.show()""")

code("""for c in D:
    objs = []
    for n, s in D[c]['sets'].items():
        objs.append(ad.AnnData(s['X'], obs=pd.DataFrame({'source': n, 'sample': s['sample'], 'group': s['group']}, index=[f'{n}|{i}' for i in range(s['X'].shape[0])]),
                               var=pd.DataFrame(index=D[c]['genes'])))
    A = ad.concat(objs); A.obs = A.obs.astype('category'); A.write_h5ad(f'{OUT}/{c}_cd14mono_cells.h5ad', compression='gzip'); print(c, A)""")

nb = nbf.v4.new_notebook(); nb.cells = cells
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
NotebookClient(nb, timeout=7200, kernel_name='python3').execute()
nbf.write(nb, f'{OUT}/monocyte_gene_level.ipynb'); print('done')
