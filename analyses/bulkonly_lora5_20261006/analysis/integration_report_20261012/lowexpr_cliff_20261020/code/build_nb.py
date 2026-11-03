"""Build lowexpr_cliff.ipynb (analyses 1 and 2). Run: python build_nb.py <out_dir>; then execute with nbclient."""
import sys, nbformat as nbf
W = sys.argv[1]; nb = nbf.v4.new_notebook(); C = nb.cells
md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Low-expression cliff in PBISC-generated cells — analyses 1 and 2 (2026-10-20)

**Question.** After normalisation, PBISC gene counts fall sharply below real from about the 1,000th–2,000th gene and for genes with real mean ≤ 0.1 UMI per cell. Is this the read-out format, or the model?

**Analysis 1 — PBISC-like read-out on real cells.** Real held-out cells are passed through read-outs and compared with PBISC on the same genes:
| source | what it is |
|---|---|
| `realU` | real held-out cells, untouched (reference) |
| `realT1` | real, binomially thinned to a depth drawn from PBISC K = 1,000 cells of the same six-type (depth only, nothing modelled) |
| `realT2` | the report's "PBISC-like read-out" (`e_checks.py`): kNN(15)-smoothed profile on the cell's own detected genes, scaled to PBISC depth, one Poisson draw |
| `realT2b` | as realT2 but support re-drawn as Bernoulli(p), p = share of the 15 neighbours detecting the gene (calibrated support) |
| `realT2m` | as realT2 but support = genes with p ≥ 0.5 (a mode-seeking support channel) |
| `k1000`, `k4000` | PBISC generated cells (Harmony-set files, all 5 draws) |

**Analysis 2 — detection rate inside cell clusters.** The x-axis is real detection inside clusters instead of the overall mean: (a) the maximum detection over real author fine clusters (n ≥ 100 cells), as planned, and (b) detection inside each of the six PBISC types, depth-matched (realT1), which compares like with like (PBISC types come from the frozen scANVI annotator).

**Data.** 30,000 random cells per source and cohort (seed 20261020). Real = `cell_transcriptome_validation_20261014/inputs/*_real_heldout_panel_counts.h5ad` (sepsis: adults only, 240,211 scoped cells; COVID: 333,503). Generated = Harmony-set files (sepsis 302,500 adult cells; COVID 410,000). Gene panel = genes present in both K files and the real file (sepsis 5,961; COVID 5,745). Per-gene statistics are computed by `code/prep.py` and stored in `data/<cohort>_gene_stats.h5ad`; this notebook only reads them. Zero-count checks use the full files, not the sample.""")
code("""import numpy as np, pandas as pd, anndata as ad, matplotlib as mpl, matplotlib.pyplot as plt, os
pd.set_option('display.width', 220); pd.set_option('display.max_columns', 30)
W = '%s'
COH = ['sepsis', 'covid']; TITLE = {'sepsis': 'Sepsis (Ye 2026)', 'covid': 'COVID (Stephenson 2025)'}
INK, INK2, MUTED, GRID = '#1f1f1e', '#5f5e5a', '#8a8984', '#e6e5e0'
COL = {'k1000': '#eb6834', 'k4000': '#eb6834', 'realT2': '#2a78d6', 'realT2m': '#1baf7a', 'realT1': '#8a8984', 'realU': '#1f1f1e', 'realT2b': '#2a78d6'}
LAB = {'k1000': 'PBISC K = 1,000', 'k4000': 'PBISC K = 4,000', 'realT2': 'real, PBISC-like read-out (T2)', 'realT2m': 'real, T2 with mode support (p ≥ 0.5)',
       'realT1': 'real, depth-matched (T1)', 'realU': 'real', 'realT2b': 'real, T2 with Bernoulli support'}
LS = {'k4000': (0, (4, 2))}
mpl.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK2, 'xtick.color': INK2,
                     'ytick.color': INK2, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
D = {}
for c in COH:
    A = ad.read_h5ad(f'{W}/data/{c}_gene_stats.h5ad'); F = pd.read_csv(f'{W}/results/{c}_fullfile_gene_totals.csv').set_index('gene')
    D[c] = (A, F)
    print(c, A.shape, A.obs.level.value_counts().to_dict())
pd.concat({c: pd.read_csv(f'{W}/results/{c}_depth.csv').set_index('source') for c in COH}).round(2)""" % W)
md("""Depth after read-out: realT2/T2b/T2m are scaled to the PBISC depth by construction; realT1 can only thin, so it ends a little shallower. Note that PBISC cells detect **more** genes per cell than realT2 (245 vs 177 in sepsis) — the PBISC deficit is not a general shortage of detected genes.""")
code("""def tables(c):
    A, F = D[c]; al = A[(A.obs.level == 'all').values]; S = list(al.obs.source)
    cp = pd.DataFrame(al.X.T, index=A.var_names, columns=S); det = pd.DataFrame(al.layers['det'].T, index=A.var_names, columns=S)
    mu = pd.DataFrame(al.layers['mean_umi'].T, index=A.var_names, columns=S)
    return cp, det, mu
SRC = ['realT1', 'realT2', 'realT2b', 'realT2m', 'k1000', 'k4000']
out1 = []
for c in COH:
    cp, det, mu = tables(c); ok = mu.realU > 0; r = cp[SRC].div(cp.realU, axis=0)[ok]
    rank = mu.realU.rank(ascending=False, method='first')[ok]
    rb = pd.cut(rank, [0, 500, 1000, 1500, 2000, 3000, 4000, 6000], labels=['1–500', '501–1,000', '1,001–1,500', '1,501–2,000', '2,001–3,000', '3,001–4,000', '>4,000'])
    mb = pd.cut(mu.realU[ok], [0, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 1e3], labels=['≤0.003', '0.003–0.01', '0.01–0.03', '0.03–0.1', '0.1–0.3', '0.3–1', '>1'])
    for name, b in [('real rank', rb), ('real mean UMI', mb)]:
        t = r.groupby(b, observed=True).median(); t.insert(0, 'n_genes', b.value_counts()); t.index = pd.MultiIndex.from_product([[c], [name], t.index])
        out1.append(t)
T1 = pd.concat(out1); T1.to_csv(f'{W}/results/a1_ratio_cp10k_by_rank_and_mean.csv'); T1.round(3)""")
md("""**Table A1.** Median per-gene ratio of mean CP10k (source / real), by the gene's real rank and real mean UMI.

- **realT1 ≈ 1 everywhere**, so depth alone does not create the drop.
- **realT2 declines gradually and monotonically** (0.8 → 0.35 → 0.25 → 0.2 → 0.16 → 0.07). This decline is built into T2: it multiplies the cell's own detections (already one sampling event) by a smoothed mean and then samples again, so sporadic genes are sampled twice. T2b (Bernoulli support) behaves the same.
- **PBISC has a different shape.** Genes ranked 500–1,500 (mean 0.1–1 UMI) are **over**-represented 1.35–1.65×. The fall comes between ranks 1,500 and 2,000 (0.2–0.3). Beyond rank 4,000 the median is exactly 0.
- **Mode support (T2m) erases genes far earlier** (≈ 0 from rank 1,500) and does not over-shoot. PBISC is neither T2 nor T2m.""")
code("""zz = []
for c in COH:
    A, F = D[c]; z = (F.real_n_det > 0) & (F.k1000_n_det == 0) & (F.k4000_n_det == 0)
    rd = F.real_n_det / F.real_n_cells
    zz.append(dict(cohort=c, real_cells=int(F.real_n_cells.iloc[0]), gen_cells_per_K=int(F.k1000_n_cells.iloc[0]),
                   genes_zero_in_both_K_full_files=int(z.sum()), of_which_measured_in_cohort=int(F.measured_gen[z].sum()),
                   real_cells_detecting_median=float(F.real_n_det[z].median()), real_cells_detecting_max=int(F.real_n_det[z].max()),
                   max_real_detection_rate=float(rd[z].max()), genes_real_det_ge_1pct_but_gen_lt_5pct_of_real=int(((rd >= 0.01) & ((F.k1000_n_det + F.k4000_n_det) / (F.k1000_n_cells + F.k4000_n_cells) < 0.05 * rd)).sum())))
    dr = A[(A.obs.level == 'draw').values]; dd = pd.DataFrame(dr.layers['det'].T, index=A.var_names)
    zz[-1].update(sample_zero_in_all_10_draws=int((dd == 0).all(axis=1).sum()), sample_zero_in_some_draws_only=int(((dd == 0).any(axis=1) & ~(dd == 0).all(axis=1)).sum()))
Z = pd.DataFrame(zz).set_index('cohort'); Z.to_csv(f'{W}/results/a1_zero_genes_fullfile.csv'); Z.T""")
md("""**Table A2 — is any gene switched off by the pipeline?** Over the full files (sepsis 302,500 and COVID 410,000 generated cells per K), about 400 genes are never generated by either K although real cells have them. All of these are very rare in real cells: at most 0.06% (sepsis) or 0.09% (COVID) of cells, typically 7–14 cells in 240,000–333,000. No gene with meaningful real expression is exactly zero, and genes that are zero in one draw are mostly not zero in others. So there is **no sign of hard masking in the pipeline**. The earlier figure of "~1,030 absent genes" came from a 20,000-cell sample and counted any count; it mixes these rare genes with sampling noise. The real loss is graded and concerns many more genes: 491 (sepsis) / 665 (COVID) genes are detected in ≥ 1% of real cells but at < 5% of that rate in PBISC.""")
code("""def six(c, layer, src):
    A, _ = D[c]; st = A[((A.obs.level == 'six_type') & (A.obs.label != 'Other')).values]
    types = sorted(set(st.obs.label[st.obs.source == 'realT1']) & set(st.obs.label[st.obs.source == 'k1000']))
    get = lambda s: np.concatenate([(st[((st.obs.source == s) & (st.obs.label == t)).values].X if layer == 'X' else st[((st.obs.source == s) & (st.obs.label == t)).values].layers[layer])[0] for t in types])
    return pd.DataFrame({s: get(s) for s in src}), types
EDGES = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9, 1.0]
out2, out4 = [], []
for c in COH:
    d, types = six(c, 'det', ['realT1', 'realT2', 'realT2b', 'realT2m', 'k1000', 'k4000'])
    b = pd.cut(d.realT1, EDGES); g = d.groupby(b, observed=True)
    t = g.median(); t.insert(0, 'n_gene_type_pairs', g.size()); t['k1000 / realT1'] = (g.k1000.median() / g.realT1.median())
    t['k4000 / realT1'] = g.k4000.median() / g.realT1.median(); t.index = pd.MultiIndex.from_product([[c], t.index.astype(str)]); out2.append(t)
    u, _ = six(c, 'mean_umi', ['realT1', 'k1000', 'k4000']); cp, _ = six(c, 'X', ['realT1', 'realT2', 'k1000', 'k4000'])
    e = pd.DataFrame({'umi_per_detection_realT1': u.realT1 / d.realT1, 'umi_per_detection_k1000': u.k1000 / d.k1000.replace(0, np.nan),
                      'umi_per_detection_k4000': u.k4000 / d.k4000.replace(0, np.nan), 'cp10k k1000 / realT1': cp.k1000 / cp.realT1,
                      'cp10k k4000 / realT1': cp.k4000 / cp.realT1, 'cp10k realT2 / realT1': cp.realT2 / cp.realT1})
    t4 = e.groupby(pd.cut(d.realT1, [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1]), observed=True).median()
    t4.index = pd.MultiIndex.from_product([[c], t4.index.astype(str)]); out4.append(t4)
T2_ = pd.concat(out2); T2_.to_csv(f'{W}/results/a2_within_type_detection_calibration.csv')
T4_ = pd.concat(out4); T4_.to_csv(f'{W}/results/a2_umi_per_detection_within_type.csv')
print('types used:', types); T2_.round(3)""")
md("""**Table A3 — detection inside each six-type, PBISC vs depth-matched real (main result of analysis 2).** Unit = gene × cell type (B, DC, monocyte, NK, T), binned by realT1 detection.

- **Detection ≥ ~5–10%:** PBISC detection is close to calibrated (0.84–0.99 × realT1 at every level up to 97%). There is no step at 0.5, so a mode-seeking support (T2m, zero below 0.1) is ruled out.
- **Detection < 5%:** a floor appears. At 2–5% PBISC detects 0.18×; at 1–2%, 0.07×; below 1%, essentially never.
- **realT2 is not calibrated anywhere below about 50%.** Its deficit is the double-sampling artefact described under Table A1, so it is a pessimistic emulator. At the very bottom (≤ 5%) it happens to agree with PBISC.""")
code("""T4_.round(3)""")
md("""**Table A4 — counts per detection inside each six-type.** Real low-detection genes are sporadic single molecules: when detected, they carry about 1.0–1.2 UMI (realT1). PBISC puts **1.7–2.1 UMI** on each detection of the same genes. Detection itself is calibrated in this range (A3), so the mean per cell is inflated. Hence mean CP10k is 1.1–1.46× real for gene × type pairs with 5–30% detection, which is the over-shoot of ranks 500–1,500 in A1.

This fits the read-out design (Bernoulli support, then a Poisson draw on the expected count). To keep detection near the real rate after the Poisson zero (e^{−μ}), the expected count μ on support must be well above 1. If μ ≈ 1–2, the mean count per detection is μ / (1 − e^{−μ}) ≈ 1.6–2.3. In real data, a low-expression detection is itself the single Poisson event. PBISC applies the support and then a second Poisson, which is the same double-sampling problem T2 exaggerates.""")
code("""out3 = []
for c in COH:
    A, F = D[c]; fi = A[(A.obs.level == 'fine').values]; xmax = pd.Series(fi.layers['det'].max(0), index=A.var_names)
    cp, det, mu = tables(c); ok = det.realU > 0
    r = det[SRC].div(det.realT1, axis=0)[ok & (det.realT1 > 0)]
    b = pd.cut(xmax[r.index], [0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]); t = r.groupby(b, observed=True).median()
    t.insert(0, 'n_genes', b.value_counts()); t.index = pd.MultiIndex.from_product([[c], t.index.astype(str)]); out3.append(t)
T3_ = pd.concat(out3); T3_.to_csv(f'{W}/results/a2_detection_ratio_by_max_fine_cluster_detection.csv'); T3_.round(3)""")
md("""**Table A5 — analysis 2 as planned.** x = the gene's highest detection in any real author fine cluster (full depth); y = gene-level detection relative to depth-matched real (source / realT1).

PBISC shows a **step between 0.2 and 0.4**: 0.15–0.18 at 0.1–0.2, 0.54 / 0.24 (sepsis / COVID) at 0.2–0.3, 0.88 / 0.73 at 0.3–0.4, and close to 1 above. The step sits lower than 0.5, and fine clusters inherit it from the six-type floor. A gene detected in 20–30% of one small fine cluster is diluted below 5% inside its six-type at PBISC depth, which puts it under the floor. So fine-state genes are the first casualties: macrophage C1Q/VSIG4, ILC KIT/RORC, cycling MKI67, low-abundance NK receptors.""")
code("""def p_rank(ax, c):
    cp, det, mu = tables(c); ok = mu.realU > 0; order = mu.realU[ok].sort_values(ascending=False).index
    x = np.arange(1, len(order) + 1)
    for s in ['realT1', 'realT2', 'realT2m', 'k1000', 'k4000']:
        r = (cp.loc[order, s] / cp.loc[order, 'realU']).clip(lower=1e-3)
        y = r.rolling(151, center=True, min_periods=40).median()
        ax.plot(x, y, color=COL[s], lw=2 if s.startswith('k') else 1.6, ls=LS.get(s, '-'), label=LAB[s])
    i01 = int((mu.realU[order] > 0.1).sum())
    ax.axvspan(i01, len(order), color=GRID, alpha=0.5, lw=0); ax.text(i01 + 40, 2.6, 'real mean ≤ 0.1 UMI', color=INK2, fontsize=8, va='top')
    for v in (1000, 2000): ax.axvline(v, color=MUTED, lw=0.8, ls=':')
    ax.axhline(1, color=INK2, lw=0.8)
    ax.set_yscale('log'); ax.set_ylim(8e-4, 3.5); ax.set_xlim(0, len(order))
    ax.set_yticks([1e-3, 1e-2, 1e-1, 1], ['0 (≤0.001)', '0.01', '0.1', '1'])
    ax.set_xlabel('genes ordered by real mean UMI (left = highest)'); ax.set_ylabel('mean CP10k, source / real (running median)')
    ax.set_title('A  Gene-level ratio along the real ranking', loc='left', fontsize=10, color=INK)

def p_cal(ax, c):
    d, _ = six(c, 'det', ['realT1', 'realT2', 'realT2m', 'k1000', 'k4000'])
    ed = np.r_[0.0007, np.geomspace(0.001, 1, 25)]; b = pd.cut(d.realT1, ed); g = d.groupby(b, observed=True)
    xm = g.realT1.median()
    q = g.k1000.quantile([0.25, 0.75]).unstack(); ax.fill_between(xm, q[0.25].clip(lower=1e-4), q[0.75].clip(lower=1e-4), color=COL['k1000'], alpha=0.15, lw=0)
    for s in ['realT2', 'realT2m', 'k1000', 'k4000']:
        ax.plot(xm, g[s].median().clip(lower=1e-4), color=COL[s], lw=2 if s.startswith('k') else 1.6, ls=LS.get(s, '-'), marker='o' if s == 'k1000' else None, ms=3.5)
    ax.plot([1e-3, 1], [1e-3, 1], color=INK2, lw=0.8); ax.axvspan(1e-3, 0.05, color=GRID, alpha=0.5, lw=0)
    ax.text(0.0012, 0.6, 'floor zone\\n(real < 5%)', color=INK2, fontsize=8, va='top')
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(1e-3, 1.05); ax.set_ylim(1e-4, 1.2)
    ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1], ['0 (≤0.0001)', '0.001', '0.01', '0.1', '1'])
    ax.set_xlabel('real detection inside the cell type (depth-matched, T1)'); ax.set_ylabel('source detection inside the same type')
    ax.set_title('B  Detection inside six-types (gene × type; line = bin median)', loc='left', fontsize=10, color=INK)

def p_fine(ax, c):
    A, F = D[c]; fi = A[(A.obs.level == 'fine').values]; xmax = pd.Series(fi.layers['det'].max(0), index=A.var_names)
    cp, det, mu = tables(c); keep = det.realT1 > 0
    ed = np.linspace(0, 1, 21); b = pd.cut(xmax[keep], ed)
    for s in ['realT2', 'realT2m', 'k1000', 'k4000']:
        r = (det[s] / det.realT1)[keep]; g = r.groupby(b, observed=True).median(); xm = xmax[keep].groupby(b, observed=True).median()
        ax.plot(xm, g.clip(lower=1e-3), color=COL[s], lw=2 if s.startswith('k') else 1.6, ls=LS.get(s, '-'), marker='o' if s == 'k1000' else None, ms=3.5)
    ax.axhline(1, color=INK2, lw=0.8); ax.set_yscale('log'); ax.set_ylim(8e-4, 2); ax.set_xlim(0, 1)
    ax.set_yticks([1e-3, 1e-2, 1e-1, 1], ['0 (≤0.001)', '0.01', '0.1', '1'])
    ax.set_xlabel('highest real detection in any author fine cluster'); ax.set_ylabel('gene detection, source / depth-matched real')
    ax.set_title('C  Analysis 2 as planned: step vs best fine cluster', loc='left', fontsize=10, color=INK)

def p_umi(ax, c):
    d, _ = six(c, 'det', ['realT1', 'k1000', 'k4000']); u, _ = six(c, 'mean_umi', ['realT1', 'k1000', 'k4000'])
    ed = np.geomspace(0.003, 1, 18); b = pd.cut(d.realT1, ed); xm = d.realT1.groupby(b, observed=True).median()
    for s in ['realT1', 'k1000', 'k4000']:
        m = d[s] > 0; y = (u[s] / d[s]).where(m).groupby(b, observed=True).median()
        ax.plot(xm, y, color=COL[s], lw=2 if s.startswith('k') else 1.6, ls=LS.get(s, '-'), marker='o' if s == 'k1000' else None, ms=3.5)
    ax.axvspan(0.003, 0.05, color=GRID, alpha=0.5, lw=0)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.003, 1.05); ax.set_ylim(0.9, 20)
    ax.set_yticks([1, 2, 5, 10, 20], ['1', '2', '5', '10', '20'])
    ax.set_xlabel('real detection inside the cell type (depth-matched, T1)'); ax.set_ylabel('UMI per detected cell (median)')
    ax.set_title('D  Counts carried by each detection', loc='left', fontsize=10, color=INK)

os.makedirs(f'{W}/figs', exist_ok=True)
for c in COH:
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.2))
    p_rank(axs[0, 0], c); p_cal(axs[0, 1], c); p_fine(axs[1, 0], c); p_umi(axs[1, 1], c)
    for ax in axs.ravel(): ax.yaxis.grid(True, color=GRID, lw=0.6); ax.set_axisbelow(True)
    hs = [mpl.lines.Line2D([], [], color=COL[s], lw=2, ls=LS.get(s, '-'), label=LAB[s]) for s in ['k1000', 'k4000', 'realT1', 'realT2', 'realT2m']]
    fig.legend(handles=hs, loc='upper center', ncol=5, frameon=False, fontsize=9, bbox_to_anchor=(0.5, 0.955))
    fig.suptitle(f'{TITLE[c]}: why low-expression genes fall away in PBISC', x=0.01, ha='left', fontsize=12.5, color=INK, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    for ext in ('png', 'svg', 'pdf'): fig.savefig(f'{W}/figs/lowexpr_cliff_{c}.{ext}', dpi=200 if ext == 'png' else None)
    plt.show()""")
md("""## Conclusion

**Analysis 1 (realT2).** The PBISC-like read-out on real cells also loses low-expression genes, but **gradually**, through double sampling (its own detections are sampled again). It does **not** reproduce PBISC's shape: over-shoot of mid genes, a fall around rank 1,500–2,000, exact zeros. The cliff is therefore not explained by "real biology plus the read-out format". It comes from how PBISC generates support and counts. Depth alone (T1) explains nothing.

**Analysis 2 (cluster detection).** There is a step, but not at 0.5:
1. **A detection floor.** Inside a cell type, PBISC is close to calibrated whenever real detection is ≥ ~5–10%. Below ~5% the gene is largely not generated (2–5% → 0.18×, 1–2% → 0.07×, < 1% → ~0). Seen from fine clusters, the step sits at a best-cluster detection of about 0.2–0.4, because fine states are diluted inside their six-type.
2. **Counts per detection are too high** (1.7–2.1 vs 1.0–1.2 UMI). Mid-detection genes are therefore inflated 1.1–1.5× in CP10k.

Both are signatures of the read-out design: Bernoulli support, then a Poisson draw on the expected count. Real low-expression detections are single Poisson events; PBISC treats them as a support state plus a count ≥ 1. It keeps them where support probability is substantial (with excess counts), but does not generate them where the per-type probability is only a few percent.

**Correction to the earlier hypothesis.** A mode-seeking support (zeroing everything below p = 0.5) is ruled out. T2m behaves completely differently, and PBISC is calibrated from 0.05 to 0.97. What remains is a floor at low probability. Whether it lives in the support logit (p̂ ≈ 0) or in μ̂ needs the saved p̂ and μ̂ before sampling (analysis 3, server checkpoint).

**Why mean 0.1 / rank ~1,500 is the boundary.** Among genes with real mean 0.03–0.1 UMI, 80% never reach 10% detection in any six-type (32–34% never reach 5%). At 0.01–0.03, 87% never reach 5%. At 0.1–0.3, almost none (0–1%) are below 5%.""")
code("""rows = []
for c in COH:
    A, F = D[c]; cp, det, mu = tables(c)
    st = A[((A.obs.level == 'six_type') & (A.obs.source == 'realT1') & (A.obs.label != 'Other')).values]; mx = pd.Series(st.layers['det'].max(0), index=A.var_names)
    for lo, hi in [(0, 0.01), (0.01, 0.03), (0.03, 0.1), (0.1, 0.3), (0.3, 1), (1, 1e3)]:
        m = (mu.realU > lo) & (mu.realU <= hi)
        rows.append(dict(cohort=c, real_mean_umi=f'{lo}–{hi}', n_genes=int(m.sum()), share_max_type_det_lt_5pct=(mx[m] < 0.05).mean(), share_max_type_det_lt_10pct=(mx[m] < 0.10).mean()))
T6 = pd.DataFrame(rows); T6.to_csv(f'{W}/results/a2_floor_share_by_real_mean.csv', index=False); T6.round(2)""")
md("""## Added 2026-10-20: is the low-gene deficit a normalisation artefact?

Hypothesis: the genes ranked just below the top 500 are over-predicted, which raises each cell's total; CP10k normalisation would then make lower genes look under-expressed. Test: compare raw UMI (depth-corrected), CP10k and detection (normalisation-free) against depth-matched real (T1), by real rank.""")
code("""rows = []
for c in COH:
    cp, det, mu = tables(c); ok = mu.realU > 0
    f = mu.k1000.sum() / mu.realT1.sum()          # PBISC total UMI / T1 total UMI
    rank = mu.realU.rank(ascending=False, method='first')
    b = pd.cut(rank, [0, 500, 1500, 2000, 3000, 4000, 7000], labels=['1-500', '501-1,500', '1,501-2,000', '2,001-3,000', '3,001-4,000', '>4,000'])
    t = pd.DataFrame({'raw UMI K1k/T1 (depth-corrected)': (mu.k1000 / mu.realT1 / f)[ok].groupby(b[ok], observed=True).median(),
                      'CP10k K1k/T1': (cp.k1000 / cp.realT1)[ok].groupby(b[ok], observed=True).median(),
                      'detection K1k/T1': (det.k1000 / det.realT1)[ok].groupby(b[ok], observed=True).median(),
                      'UMI share T1 %': mu.realT1.groupby(b, observed=True).sum() / mu.realT1.sum() * 100,
                      'UMI share K1k %': mu.k1000.groupby(b, observed=True).sum() / mu.k1000.sum() * 100})
    t['depth factor K1k/T1'] = f; t.index = pd.MultiIndex.from_product([[c], t.index]); rows.append(t)
T7 = pd.concat(rows); T7.to_csv(f'{W}/results/a3_raw_vs_cp10k_vs_detection_by_rank.csv'); T7.round(3)""")
md("""**Table A7.** Raw UMI (divided by the 1.27-1.29x depth factor) and CP10k agree to within 2-3% at every rank, and detection, which normalisation cannot touch, is 10-25x lower beyond rank 2,000. The over-predicted genes (rank 501-1,500) hold 7.2 -> 9.0% (sepsis) and 8.0 -> 10.9% (COVID) of all UMI, so normalisation shrinks other genes by only 2-3%. **The low-gene deficit is real under-generation, not a normalisation artefact; the over-prediction and the floor are separate defects.**""")
nbf.write(nb, f'{W}/lowexpr_cliff.ipynb')
