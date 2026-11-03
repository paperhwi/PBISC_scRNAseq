"""Build lowexpr_floor.ipynb. Run: python build_nb.py <out_dir>; execute with nbclient. Conclusion text is read from code/conclusion.md if present."""
import sys, os, nbformat as nbf
W = sys.argv[1]; nb = nbf.v4.new_notebook(); C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Floor read-out vs previous K4000: low-expression cliff gate M1–M4 (2026-10-21)

**Question.** The new inference style (`depth256-K4000-floor1e-4`, draw0: bulk ×256 into the conditioner, no Poisson draw, integer rule `x<1e-4→0, 1e-4≤x<1→1, x≥1→floor(x)`) was built to cure the low-expression cliff (2026-10-20: detection floor below ~5% within a type, 1.7–2.1 UMI per detection, over-shoot at ranks 500–1,500). Does it pass gates M1–M4 of the model-fix plan?

| gate | definition (gene × six-type, depth-matched real = `realT1`) | pass |
|---|---|---|
| M1 | detection ratio (source / realT1) where realT1 detection is 2–5% | ≥ 0.7 (previous 0.18) |
| M2 | same, realT1 detection 0.5–2% | ≥ 0.5 |
| M3 | UMI per detection, genes at 5–30% realT1 detection | 1.0–1.4 (real 1.0–1.2; previous 1.7–2.1) |
| M4 | per-gene CP10k ratio to real: ranks 500–1,500 within 0.85–1.15; every bin beyond rank 2,000 ≥ 0.5 | both |

M5 (rare-state markers) is only approximated here by marker-gene detection in all cells; M6 (composition, logFC, DEG, iLISI, RF realism) needs the cell-level analyses and is **not tested**.

**Compared.** `prev K4000` = Harmony-set K=4,000 files from 2026-10-06 (Poisson read-out, 5 draws, 500 cells/sample/draw, adapters retrained at K=4,000); `floor` = the new file (single draw0, 4,000 cells/sample, deterministic floor read-out). They differ in several things at once (input depth ×256, read-out, draw design), so a difference cannot be attributed to the floor rule alone. Each is compared with real cells binomially thinned to its own depth. Same real sample (seed 20261020, reproduced to 3e-8), same gene panel (sepsis 5,961, COVID 5,745), 30,000 cells per source; sepsis adults only; COVID without the 5 adapter-training samples in the floor file.""")
code("""import numpy as np, pandas as pd, anndata as ad, matplotlib as mpl, matplotlib.pyplot as plt, os
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 40)
W = '%s'; PREV = W + '/../lowexpr_cliff_20261020'
COH = ['sepsis', 'covid']; TITLE = {'sepsis': 'Sepsis (Ye 2026)', 'covid': 'COVID (Imperial ESKD)'}
INK, INK2, MUTED, GRID = '#1f1f1e', '#5f5e5a', '#8a8984', '#e6e5e0'
COL = {'prev': '#eb6834', 'floor': '#1b8f6a', 'realT1': '#8a8984'}; LAB = {'prev': 'prev K4000 (Poisson read-out)', 'floor': 'floor read-out (draw0)', 'realT1': 'real, depth-matched'}
mpl.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
D = {}
for c in COH:
    A0 = ad.read_h5ad(f'{PREV}/data/{c}_gene_stats.h5ad'); A1 = ad.read_h5ad(f'{W}/data/{c}_floor_gene_stats.h5ad'); A0 = A0[:, A1.var_names].copy()
    D[c] = (A0, A1); print(c, 'genes', A1.n_vars, '| floor obs', A1.obs.level.value_counts().to_dict())
dep = pd.concat({c: pd.read_csv(f'{W}/results/{c}_floor_depth.csv').set_index('source') for c in COH}); dep.round(3)""" % W)
md("""**Depth and per-cell shape.** The floor cells are at real depth, but they detect more genes per cell than depth-matched real cells and carry more single-UMI entries; see the table (previous K4000: sepsis 274 genes/cell and 29% of non-zero entries equal 1; real 404 / 53%).""")
code("""def lv(A, src, layer):
    m = ((A.obs.level == 'all') & (A.obs.source == src)).values; return (A.X if layer == 'X' else A.layers[layer])[np.where(m)[0][0]]
def six(A, srcs, layer, types):
    out = {}
    for s in srcs:
        out[s] = np.concatenate([(A.X if layer == 'X' else A.layers[layer])[np.where(((A.obs.level == 'six_type') & (A.obs.source == s) & (A.obs.label == t)).values)[0][0]] for t in types])
    return pd.DataFrame(out)
def types_of(A0, A1):
    f = lambda A, s: set(A.obs.label[(A.obs.level == 'six_type') & (A.obs.source == s)])
    return sorted((f(A0, 'realT1') & f(A0, 'k4000') & f(A1, 'realT1') & f(A1, 'floor')) - {'Other'})
PAIRS = {}   # per cohort: gene x type table with det / umi per source
for c in COH:
    A0, A1 = D[c]; T = types_of(A0, A1)
    d0 = six(A0, ['realT1', 'k4000'], 'det', T); u0 = six(A0, ['realT1', 'k4000'], 'mean_umi', T)
    d1 = six(A1, ['realT1', 'floor'], 'det', T); u1 = six(A1, ['realT1', 'floor'], 'mean_umi', T)
    PAIRS[c] = pd.DataFrame({'det_T1_prev': d0.realT1, 'det_prev': d0.k4000, 'umi_T1_prev': u0.realT1, 'umi_prev': u0.k4000,
                             'det_T1_new': d1.realT1, 'det_floor': d1.floor, 'umi_T1_new': u1.realT1, 'umi_floor': u1.floor})
    print(c, 'types:', T, '| gene x type pairs:', len(PAIRS[c]))""")
md("""## M1 / M2 — detection inside a type, by real (depth-matched) detection

Ratio = median source detection / median realT1 detection in each bin (the statistic used on 2026-10-20), plus the ratio of summed detections (less sensitive to a few tiny values). Bins are defined on each comparison's own realT1 detection.""")
code("""EDG = [0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0]
rows = []
for c in COH:
    P = PAIRS[c]
    for name, (r, s) in {'prev K4000': ('det_T1_prev', 'det_prev'), 'floor': ('det_T1_new', 'det_floor')}.items():
        b = pd.cut(P[r], EDG, include_lowest=True); g = P.groupby(b, observed=True)
        t = pd.DataFrame({'n_pairs': g.size(), 'real_T1_det_median': g[r].median(), 'src_det_median': g[s].median(),
                          'ratio_of_medians': g[s].median() / g[r].median(), 'ratio_of_sums': g[s].sum() / g[r].sum(),
                          'share_pairs_src_zero': g[s].apply(lambda x: (x == 0).mean()), 'share_pairs_ratio_gt_1.3': g.apply(lambda x: (x[s] > 1.3 * x[r]).mean())})
        t.insert(0, 'cohort', c); t.insert(1, 'source', name); t.index = t.index.astype(str); rows.append(t)
T_det = pd.concat(rows); T_det.index.name = 'real_T1_detection_bin'; T_det.to_csv(f'{W}/results/b1_within_type_detection_ratio.csv')
T_det.round(3)""")
md("""**Table B1.** Inside each six-type, detection of PBISC relative to depth-matched real. Bins below 5% are the former floor zone.""")
code("""def gate_det(P, r, s, lo, hi):
    m = (P[r] > lo) & (P[r] <= hi); return dict(med=P.loc[m, s].median() / P.loc[m, r].median(), sums=P.loc[m, s].sum() / P.loc[m, r].sum(), n=int(m.sum()))
def umi_det(P, r, s, lo, hi):
    m = (P[r] > lo) & (P[r] <= hi); ur = P['umi_' + r.split('det_')[1]]; us = P['umi_' + s.split('det_')[1]] if False else None
GATE = []
for c in COH:
    P = PAIRS[c]
    for nm, (r, s, ur, us) in {'prev K4000': ('det_T1_prev', 'det_prev', 'umi_T1_prev', 'umi_prev'), 'floor': ('det_T1_new', 'det_floor', 'umi_T1_new', 'umi_floor')}.items():
        g1 = gate_det(P, r, s, 0.02, 0.05); g2 = gate_det(P, r, s, 0.005, 0.02)
        m3 = (P[r] > 0.05) & (P[r] <= 0.30) & (P[s] > 0)
        GATE.append(dict(cohort=c, source=nm, M1_ratio_2to5pct=g1['med'], M1_sums=g1['sums'], M1_n=g1['n'], M2_ratio_0p5to2pct=g2['med'], M2_sums=g2['sums'], M2_n=g2['n'],
                         M3_umi_per_det_src=(P.loc[m3, us] / P.loc[m3, s]).median(), M3_umi_per_det_real=(P.loc[(P[r] > 0.05) & (P[r] <= 0.30), ur] / P.loc[(P[r] > 0.05) & (P[r] <= 0.30), r]).median()))
G1 = pd.DataFrame(GATE).set_index(['cohort', 'source']); G1.round(3)""")
md("""## M3 — counts carried by each detection (inside a type)""")
code("""rows = []
for c in COH:
    P = PAIRS[c]
    for nm, (r, s, ur, us) in {'prev K4000': ('det_T1_prev', 'det_prev', 'umi_T1_prev', 'umi_prev'), 'floor': ('det_T1_new', 'det_floor', 'umi_T1_new', 'umi_floor')}.items():
        b = pd.cut(P[r], [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1]); ok = P[s] > 0
        t = pd.DataFrame({'umi_per_det_real_T1': (P[ur] / P[r]).groupby(b, observed=True).median(), 'umi_per_det_src': (P[us] / P[s].where(ok)).groupby(b, observed=True).median()})
        t.insert(0, 'cohort', c); t.insert(1, 'source', nm); t.index = t.index.astype(str); rows.append(t)
T_umi = pd.concat(rows); T_umi.index.name = 'real_T1_detection_bin'; T_umi.to_csv(f'{W}/results/b3_umi_per_detection_within_type.csv'); T_umi.round(3)""")
md("""## M4 — gene-level CP10k ratio along the real ranking""")
code("""RB = [0, 500, 1000, 1500, 2000, 3000, 4000, 7000]; RL = ['1–500', '501–1,000', '1,001–1,500', '1,501–2,000', '2,001–3,000', '3,001–4,000', '>4,000']
rows = []; CP = {}
for c in COH:
    A0, A1 = D[c]; mu = lv(A1, 'realU', 'mean_umi'); ok = mu > 0
    rank = pd.Series(mu).rank(ascending=False, method='first').to_numpy(); b = pd.cut(rank, RB, labels=RL)
    cp = {'realU': lv(A1, 'realU', 'X'), 'prev': lv(A0, 'k4000', 'X'), 'floor': lv(A1, 'floor', 'X'), 'T1prev': lv(A0, 'realT1', 'X'), 'T1new': lv(A1, 'realT1', 'X')}
    dt = {'realU': lv(A1, 'realU', 'det'), 'prev': lv(A0, 'k4000', 'det'), 'floor': lv(A1, 'floor', 'det'), 'T1prev': lv(A0, 'realT1', 'det'), 'T1new': lv(A1, 'realT1', 'det')}
    CP[c] = dict(cp=cp, det=dt, rank=rank, ok=ok, mu=mu)
    for nm, s in [('prev K4000', 'prev'), ('floor', 'floor'), ('realT1 (floor depth)', 'T1new')]:
        r = pd.Series(cp[s] / np.where(cp['realU'] > 0, cp['realU'], np.nan))[ok]; rd = pd.Series(dt[s] / np.where(dt['realU'] > 0, dt['realU'], np.nan))[ok]
        t = pd.DataFrame({'n_genes': pd.Series(b[ok]).value_counts(), 'cp10k_ratio_median': r.groupby(b[ok], observed=True).median(), 'detection_ratio_median': rd.groupby(b[ok], observed=True).median(),
                          'cp10k_ratio_p25': r.groupby(b[ok], observed=True).quantile(0.25), 'cp10k_ratio_p75': r.groupby(b[ok], observed=True).quantile(0.75)})
        t.insert(0, 'cohort', c); t.insert(1, 'source', nm); rows.append(t)
T_rank = pd.concat(rows); T_rank.index.name = 'real_rank_bin'; T_rank.to_csv(f'{W}/results/b4_cp10k_ratio_by_rank.csv'); T_rank.round(3)""")
code("""rows = []
for c in COH:
    z = CP[c]; ok = z['ok']; rank = z['rank']
    for nm, s in [('prev K4000', 'prev'), ('floor', 'floor')]:
        r = pd.Series(z['cp'][s] / np.where(z['cp']['realU'] > 0, z['cp']['realU'], np.nan))
        m5 = ok & (rank > 500) & (rank <= 1500); late = [r[ok & (rank > lo) & (rank <= hi)].median() for lo, hi in [(2000, 3000), (3000, 4000), (4000, 7000)]]
        rows.append(dict(cohort=c, source=nm, M4_ratio_rank500_1500=r[m5].median(), M4_ratio_2001_3000=late[0], M4_ratio_3001_4000=late[1], M4_ratio_gt4000=late[2],
                         M4_ratio_gt2000_all=r[ok & (rank > 2000)].median()))
G2 = pd.DataFrame(rows).set_index(['cohort', 'source']); G2.round(3)""")
md("""## Gate summary""")
code("""R = G1.join(G2)
def verdict(row):
    return pd.Series({'M1 (>=0.7)': 'PASS' if row.M1_ratio_2to5pct >= 0.7 else 'FAIL', 'M2 (>=0.5)': 'PASS' if row.M2_ratio_0p5to2pct >= 0.5 else 'FAIL',
                      'M3 (1.0-1.4)': 'PASS' if 1.0 <= row.M3_umi_per_det_src <= 1.4 else 'FAIL',
                      'M4 (0.85-1.15 and >=0.5)': 'PASS' if (0.85 <= row.M4_ratio_rank500_1500 <= 1.15) and min(row.M4_ratio_2001_3000, row.M4_ratio_3001_4000, row.M4_ratio_gt4000) >= 0.5 else 'FAIL'})
V = R.apply(verdict, axis=1); GS = R[['M1_ratio_2to5pct', 'M2_ratio_0p5to2pct', 'M3_umi_per_det_src', 'M4_ratio_rank500_1500', 'M4_ratio_2001_3000', 'M4_ratio_3001_4000', 'M4_ratio_gt4000']].join(V)
GS.to_csv(f'{W}/results/m_gate_summary.csv'); R.to_csv(f'{W}/results/m_gate_values.csv'); GS.round(3)""")
md("""## Over-shoot check (the gates only bound the under-generation side)

A deterministic floor that turns any continuous value ≥ 1e-4 into 1 can create detections that real cells do not have. Detection ratio (source / depth-matched real) for gene × type pairs in the lowest real-detection bins, and the share of pairs where the source detects more than 1.3× real, are in Table B1 (`share_pairs_ratio_gt_1.3`). Below, the cell-level consequence: genes detected per cell and the share of single-UMI entries.""")
code("""rows = []
for c in COH:
    A0, A1 = D[c]; d0 = pd.read_csv(f'{PREV}/results/{c}_depth.csv').set_index('source'); d1 = pd.read_csv(f'{W}/results/{c}_floor_depth.csv').set_index('source')
    for nm, d, s in [('real (full depth)', d1, 'realU'), ('real depth-matched to floor', d1, 'realT1'), ('floor', d1, 'floor'), ('real depth-matched to prev K4000', d0, 'realT1'), ('prev K4000', d0, 'k4000')]:
        rows.append(dict(cohort=c, source=nm, umi_mean=d.loc[s, 'umi_mean'], detected_genes_median=d.loc[s, 'detected_median'], share_entries_equal_1=d.loc[s, 'frac_ones']))
pd.DataFrame(rows).set_index(['cohort', 'source']).round(3).to_csv(f'{W}/results/b5_cell_level_shape.csv'); pd.read_csv(f'{W}/results/b5_cell_level_shape.csv', index_col=[0, 1])""")
md("""## M5 proxy — fine-state marker genes and the best-fine-cluster step

(a) Gene-level detection in all cells (source / realT1) for marker sets of the states that were absent (macrophage, ILC, cycling) and for low-abundance NK receptors. This is only a proxy: the generated composition differs and the frozen scANVI labels have no macrophage/ILC class. (b) Analysis 2 as planned: detection ratio against the gene's highest detection in any real author fine cluster (full-depth real, from the 2026-10-20 stats).""")
code("""SETS = {'macrophage': ['C1QA', 'C1QB', 'C1QC', 'VSIG4', 'CD163', 'MSR1', 'NR1H3'], 'ILC': ['KIT', 'RORC', 'IL1RL1', 'HPGDS', 'CD200R1'],
        'cycling': ['MKI67', 'TOP2A', 'TYMS', 'BIRC5', 'CENPF'], 'NK receptors': ['NCAM1', 'NCR1', 'KIR2DL1', 'KIR3DL1']}
rows = []
for c in COH:
    A0, A1 = D[c]; z = CP[c]; gi = {g: i for i, g in enumerate(A1.var_names)}
    for st, gl in SETS.items():
        for g in gl:
            if g not in gi: rows.append(dict(cohort=c, set=st, gene=g, in_panel=False)); continue
            i = gi[g]; d = {k: float(v[i]) for k, v in z['det'].items()}
            rows.append(dict(cohort=c, set=st, gene=g, in_panel=True, real_full_depth=d['realU'], real_T1_floor_depth=d['T1new'], prev_K4000=d['prev'], floor=d['floor'],
                             prev_over_T1=d['prev'] / d['T1prev'] if d['T1prev'] > 0 else np.nan, floor_over_T1=d['floor'] / d['T1new'] if d['T1new'] > 0 else np.nan))
M5 = pd.DataFrame(rows); M5.to_csv(f'{W}/results/b6_marker_gene_detection_all_cells.csv', index=False)
print(M5.groupby(['cohort', 'set'])[['prev_over_T1', 'floor_over_T1']].median().round(3).to_string()); M5.round(4)""")
code("""rows = []
EF = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
for c in COH:
    A0, A1 = D[c]; fi = A0[((A0.obs.level == 'fine') & (A0.obs.source == 'realU')).values]; xmax = pd.Series(fi.layers['det'].max(0), index=A0.var_names)
    z = CP[c]; ok = (z['det']['T1new'] > 0) & (z['det']['T1prev'] > 0) & (z['det']['realU'] > 0)
    b = pd.cut(xmax.to_numpy(), EF)
    for nm, s, t1 in [('prev K4000', 'prev', 'T1prev'), ('floor', 'floor', 'T1new')]:
        r = pd.Series(z['det'][s] / np.where(z['det'][t1] > 0, z['det'][t1], np.nan))[ok]
        t = pd.DataFrame({'n_genes': pd.Series(b[ok]).value_counts().sort_index(), 'detection_ratio_median': r.groupby(b[ok], observed=True).median()}); t.insert(0, 'cohort', c); t.insert(1, 'source', nm); t.index = t.index.astype(str); rows.append(t)
T_fine = pd.concat(rows); T_fine.index.name = 'max_real_fine_cluster_detection'; T_fine.to_csv(f'{W}/results/b7_detection_ratio_by_best_fine_cluster.csv'); T_fine.round(3)""")
md("""## Specificity check — is the extra detection on-target?

For genes that are type-specific in real cells (best six-type detection ≥ 20% and mean detection in the other types ≤ 25% of it), the table gives detection inside the best type (on-target) and in the other types (off-target) relative to depth-matched real. Over-detection that sits **off-target** is diffuse noise; over-detection **on-target** is a calibration (threshold) problem. Marker rows list canonical genes.""")
code("""TYPES = lambda A, s: [t for t in ['B cell', 'DC', 'Monocyte', 'NK cell', 'T cell'] if ((A.obs.level == 'six_type') & (A.obs.source == s) & (A.obs.label == t)).any()]
def dmat(A, s, T): return np.vstack([A.layers['det'][np.where(((A.obs.level == 'six_type') & (A.obs.source == s) & (A.obs.label == t)).values)[0][0]] for t in T])
MK = {'NK cell': ['NCAM1', 'NCR1', 'KLRF1', 'KIR2DL1', 'KIR3DL1'], 'T cell': ['CD3E', 'CD3D', 'IL7R'], 'B cell': ['MS4A1', 'CD79A'], 'Monocyte': ['CD14', 'LYZ', 'S100A8']}
rows, mrows = [], []
for c in COH:
    A0, A1 = D[c]; T = sorted(set(TYPES(A0, 'realT1')) & set(TYPES(A0, 'k4000')) & set(TYPES(A1, 'realT1')) & set(TYPES(A1, 'floor')))
    R = {'prev': (dmat(A0, 'realT1', T), dmat(A0, 'k4000', T)), 'floor': (dmat(A1, 'realT1', T), dmat(A1, 'floor', T))}
    for nm, (r, s) in R.items():
        best = r.argmax(0); mx = r.max(0); oth = (r.sum(0) - mx) / (len(T) - 1); q = (mx >= 0.2) & (oth <= 0.25 * mx); idx = np.where(q)[0]
        on_r, on_s = r[best[idx], idx], s[best[idx], idx]
        mo_r = np.array([(r[:, i].sum() - r[best[i], i]) / (len(T) - 1) for i in idx]); mo_s = np.array([(s[:, i].sum() - s[best[i], i]) / (len(T) - 1) for i in idx])
        rows.append(dict(cohort=c, source=nm, n_type_specific_genes=len(idx), on_target_det_ratio=np.median(on_s / on_r), off_target_det_ratio_of_sums=mo_s.sum() / mo_r.sum(),
                         specificity_real=np.median(on_r / np.maximum(mo_r, 1e-4)), specificity_src=np.median(on_s / np.maximum(mo_s, 1e-4))))
        gi = {g: i for i, g in enumerate(A1.var_names)}
        for tt, gl in MK.items():
            for g in gl:
                if g in gi and tt in T:
                    i = gi[g]; j = T.index(tt); o = [k for k in range(len(T)) if k != j]
                    mrows.append(dict(cohort=c, source=nm, gene=g, target=tt, on_target_real=r[j, i], on_target_src=s[j, i], off_target_real=r[o, i].mean(), off_target_src=s[o, i].mean()))
SP = pd.DataFrame(rows).set_index(['cohort', 'source']); SP.to_csv(f'{W}/results/b8_type_specificity.csv'); MKt = pd.DataFrame(mrows); MKt.to_csv(f'{W}/results/b8_marker_on_off_target.csv', index=False)
print(SP.round(3).to_string()); MKt.round(4)""")
md("""## Figure""")
code("""def run_med(x, w=151): return pd.Series(x).rolling(w, center=True, min_periods=40).median().to_numpy()
def p_rank(ax, c):
    z = CP[c]; ok = z['ok']; order = np.argsort(-z['mu'])[:int(ok.sum())]
    for s in ['prev', 'floor']:
        r = np.clip(z['cp'][s][order] / z['cp']['realU'][order], 1e-3, None); ax.plot(np.arange(1, len(order) + 1), run_med(r), color=COL[s], lw=2, ls='--' if s == 'prev' else '-', label=LAB[s])
    r = np.clip(z['cp']['T1new'][order] / z['cp']['realU'][order], 1e-3, None); ax.plot(np.arange(1, len(order) + 1), run_med(r), color=COL['realT1'], lw=1.4, label=LAB['realT1'])
    ax.axhline(1, color=INK2, lw=0.8); ax.axhspan(0.85, 1.15, color=GRID, alpha=0.6, lw=0); ax.axvspan(500, 1500, color=GRID, alpha=0.3, lw=0); ax.axvline(2000, color=MUTED, lw=0.8, ls=':')
    ax.set_yscale('log'); ax.set_ylim(8e-4, 4); ax.set_yticks([1e-3, 1e-2, 1e-1, 1], ['0 (≤0.001)', '0.01', '0.1', '1']); ax.set_xlim(0, len(order))
    ax.set_xlabel('genes ordered by real mean UMI'); ax.set_ylabel('mean CP10k, source / real (running median)'); ax.set_title('A  M4: gene-level ratio along the real ranking', loc='left', fontsize=10, color=INK)
def p_cal(ax, c):
    P = PAIRS[c]; ed = np.r_[0.0007, np.geomspace(0.001, 1, 25)]
    for nm, (r, s) in {'prev': ('det_T1_prev', 'det_prev'), 'floor': ('det_T1_new', 'det_floor')}.items():
        g = P.groupby(pd.cut(P[r], ed), observed=True); xm = g[r].median(); ax.plot(xm, (g[s].median()).clip(lower=1e-4), color=COL[nm], lw=2, ls='--' if nm == 'prev' else '-', marker='o', ms=3)
    ax.plot([1e-3, 1], [1e-3, 1], color=INK2, lw=0.8); ax.axvspan(1e-3, 0.05, color=GRID, alpha=0.5, lw=0)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-3, 1.05); ax.set_ylim(1e-4, 1.5); ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1], ['0 (≤0.0001)', '0.001', '0.01', '0.1', '1'])
    ax.set_xlabel('real detection inside the type (depth-matched)'); ax.set_ylabel('source detection inside the same type'); ax.set_title('B  M1/M2: detection inside six-types (bin median)', loc='left', fontsize=10, color=INK)
def p_fine(ax, c):
    sub = T_fine[T_fine.cohort == c]
    for nm, k in [('prev K4000', 'prev'), ('floor', 'floor')]:
        t = sub[sub.source == nm]; x = [float(i.split(',')[1].strip(' ]')) - 0.025 for i in t.index]; ax.plot(x, t.detection_ratio_median.clip(lower=1e-3), color=COL[k], lw=2, ls='--' if k == 'prev' else '-', marker='o', ms=3.5)
    ax.axhline(1, color=INK2, lw=0.8); ax.set_yscale('log'); ax.set_ylim(8e-4, 4); ax.set_yticks([1e-3, 1e-2, 1e-1, 1], ['0 (≤0.001)', '0.01', '0.1', '1'])
    ax.set_xlabel('highest real detection in any author fine cluster (bin upper edge)'); ax.set_ylabel('gene detection, source / depth-matched real'); ax.set_title('C  Fine-state genes: step vs best fine cluster', loc='left', fontsize=10, color=INK)
def p_umi(ax, c):
    P = PAIRS[c]; ed = np.geomspace(0.003, 1, 18)
    for nm, (r, s, ur, us) in {'prev': ('det_T1_prev', 'det_prev', 'umi_T1_prev', 'umi_prev'), 'floor': ('det_T1_new', 'det_floor', 'umi_T1_new', 'umi_floor')}.items():
        b = pd.cut(P[r], ed); xm = P[r].groupby(b, observed=True).median(); ax.plot(xm, (P[us] / P[s].where(P[s] > 0)).groupby(b, observed=True).median(), color=COL[nm], lw=2, ls='--' if nm == 'prev' else '-', marker='o', ms=3)
    b = pd.cut(P['det_T1_new'], ed); ax.plot(P['det_T1_new'].groupby(b, observed=True).median(), (P['umi_T1_new'] / P['det_T1_new']).groupby(b, observed=True).median(), color=COL['realT1'], lw=1.4)
    ax.axhspan(1.0, 1.4, color=GRID, alpha=0.6, lw=0); ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.003, 1.05); ax.set_ylim(0.9, 20); ax.set_yticks([1, 2, 5, 10, 20], ['1', '2', '5', '10', '20'])
    ax.set_xlabel('real detection inside the type (depth-matched)'); ax.set_ylabel('UMI per detected cell (median)'); ax.set_title('D  M3: counts carried by each detection', loc='left', fontsize=10, color=INK)
os.makedirs(f'{W}/figs', exist_ok=True)
for c in COH:
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.2)); p_rank(axs[0, 0], c); p_cal(axs[0, 1], c); p_fine(axs[1, 0], c); p_umi(axs[1, 1], c)
    for ax in axs.ravel(): ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    fig.legend(handles=[mpl.lines.Line2D([], [], color=COL[k], lw=2, ls='--' if k == 'prev' else '-', label=LAB[k]) for k in ['prev', 'floor', 'realT1']], loc='upper center', ncol=3, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle(f'{TITLE[c]}: floor read-out vs previous K4000', x=0.01, ha='left', fontsize=12.5, color=INK, y=0.995); fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ('png', 'svg', 'pdf'): fig.savefig(f'{W}/figs/lowexpr_floor_{c}.{ext}', dpi=200 if ext == 'png' else None)
    plt.show()""")
cf = os.path.join(W, 'code', 'conclusion.md')
if os.path.exists(cf): md(open(cf, encoding='utf-8').read())
nbf.write(nb, f'{W}/lowexpr_floor.ipynb')
