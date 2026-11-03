import nbformat as nbf
C = []; md = lambda s: C.append(nbf.v4.new_markdown_cell(s)); code = lambda s: C.append(nbf.v4.new_code_cell(s))
md("""# Why does PBISC K=4,000 integrate poorly with real cells? The 2,000 Harmony HVGs (2026-10-22)

The 6,000-gene count curve shows PBISC K=4,000 matching real cells over the top ~1,000 genes. Harmony does not use those 6,000 genes, though.
It uses 2,000 seurat_v3 HVGs picked on real half A, log1p CP10k, scaled by real-A mean / SD (clip 10), PCA 30 fit on real A, then Harmony (batch = source).
This notebook follows those exact genes through each step.

1. Where the 2,000 HVGs sit on the real expression ranking.
2. Gene-count and detection curves for the HVGs only (real A, real B, PBISC K=4,000, PBISC floor).
3. What Harmony actually sees: the per-HVG shift of the generated cells in real-A SD units, and their SD ratio.
4. Which genes create the generated-vs-real offset in PCA space (attribution of the centroid offset by gene rank band).
5. Causal check: re-run the same Harmony + iLISI pipeline using only the HVGs within the top 1,000 real-ranked genes, or only those beyond.

Data: ~60,000-cell bundles per source (`ilisi_floor_20261021/data/big`). Real A = reference patients; real B = other patients (real-vs-real control).""")
code(r'''import os, sys, numpy as np, pandas as pd, anndata as ad, matplotlib as mpl, matplotlib.pyplot as plt, harmonypy
from sklearn.decomposition import PCA
sys.path.insert(0, "/home/claude/work/s3/code"); import ilisi_lib as L
IN = "/home/claude/work/s3/big"; OUT = "/home/claude/work/hvg"; FIG, RES = f"{OUT}/figs", f"{OUT}/results"
for d in (FIG, RES): os.makedirs(d, exist_ok=True)
COH = ["sepsis", "covid"]; SRC = ["ref", "realU", "prev4k", "floor"]
LAB = {"ref": "Real A (reference)", "realU": "Real B (other patients)", "prev4k": "PBISC K=4,000", "floor": "PBISC floor"}
COL = {"ref": "#4A4F5C", "realU": "#9AA3B2", "prev4k": "#eb6834", "floor": "#2a78d6"}
INK, INK2, MUTED, GRID = "#1F2430", "#4A5162", "#8A90A0", "#E4E6EB"
mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#B8BDC7", "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none", "pdf.fonttype": 42})
D = {(c, s): ad.read_h5ad(f"{IN}/{c}_{s}.h5ad") for c in COH for s in SRC}
{k: v.shape for k, v in D.items()}''')
md("## 1–3. Per-gene statistics: raw counts, then the scaled values Harmony sees")
code(r'''G = {}; HV = {}; SC = {}
for c in COH:
    R = D[(c, "ref")]; genes = np.array(R.var_names); XR = R.X.tocsr()
    hv = L.hvg_seurat_v3([XR]); HV[c] = hv                                  # identical call to big_joint.py
    YR = L.lognorm(XR)[:, hv].toarray().astype(np.float32); mu, sd = YR.mean(0), YR.std(0) + 1e-6; SC[c] = (mu, sd)
    t = pd.DataFrame({"gene": genes}); t["is_hvg"] = False; t.loc[hv, "is_hvg"] = True
    for s in SRC:
        X = D[(c, s)].X.tocsr(); n = X.shape[0]
        t[f"{s}_mean_umi"] = np.asarray(X.sum(0)).ravel() / n; t[f"{s}_det"] = np.asarray((X > 0).sum(0)).ravel() / n
    t["real_rank"] = t.ref_mean_umi.rank(ascending=False, method="first").astype(int)
    for s in SRC:                                                           # scaled log1p CP10k on the HVGs (clip 10), as fed to PCA
        Z = np.clip((L.lognorm(D[(c, s)].X.tocsr())[:, hv].toarray().astype(np.float32) - mu) / sd, -10, 10)
        t.loc[hv, f"{s}_z_mean"] = Z.mean(0); t.loc[hv, f"{s}_z_sd"] = Z.std(0)
        SC[(c, s)] = Z[np.sort(np.random.default_rng(1).choice(len(Z), min(25000, len(Z)), replace=False))]; del Z   # keep a 25k-cell subsample (memory)
    G[c] = t; t.to_csv(f"{RES}/{c}_hvg_gene_stats.csv", index=False)
BANDS = [(1, 500), (501, 1000), (1001, 2000), (2001, 3000), (3001, 6000)]
def band(r): return next(f"{a:,}–{b:,}" for a, b in BANDS if a <= r <= b)
rows = []
for c in COH:
    h = G[c][G[c].is_hvg].copy(); h["band"] = [band(r) for r in h.real_rank]
    for b, d in h.groupby("band", sort=False):
        rows.append({"cohort": c, "real_rank_band": b, "n_hvg": len(d), "share_of_hvg": len(d) / len(h),
                     **{f"{s}_mean_umi_ratio_to_realA": float(np.median(d[f"{s}_mean_umi"] / d.ref_mean_umi)) for s in SRC[1:]},
                     **{f"{s}_det_ratio_to_realA": float(np.median(d[f"{s}_det"] / d.ref_det)) for s in SRC[1:]},
                     **{f"{s}_median_z_shift": float(np.median(d[f"{s}_z_mean"])) for s in SRC[1:]},
                     **{f"{s}_median_sd_ratio": float(np.median(d[f"{s}_z_sd"] / d.ref_z_sd)) for s in SRC[1:]}})
B = pd.DataFrame(rows); B["order"] = B.real_rank_band.map({f"{a:,}–{b:,}": i for i, (a, b) in enumerate(BANDS)})
B = B.sort_values(["cohort", "order"]).drop(columns="order"); B.to_csv(f"{RES}/hvg_by_rank_band.csv", index=False); B.round(3)''')
md("## 4. Which genes make the generated cells sit away from real cells in PCA space?\nPCA 30 fit on real A (scaled HVGs, as in the pipeline; 25,000-cell subsample for the fit, full-data gene means for the offset). Offset = centroid of the source minus centroid of real A, in PC space. Each gene's share is its projection onto that offset, so the shares sum to 1. Offset size is given relative to the real-A spread: √(Σ PC variance).")
code(r'''arows = []; PCAS = {}
for c in COH:
    pca = PCA(30, random_state=0).fit(SC[(c, "ref")]); PCAS[c] = pca; Wt = pca.components_            # 30 x 2000
    spread = np.sqrt(pca.explained_variance_.sum()); rank = G[c].real_rank.to_numpy()[HV[c]]
    for s in SRC[1:]:
        dz = (G[c].loc[HV[c], f"{s}_z_mean"] - G[c].loc[HV[c], "ref_z_mean"]).to_numpy(); d = Wt @ dz   # offset in PC space (full-data means)
        share = (dz[None, :] * Wt * d[:, None]).sum(0) / (d @ d)                                       # per-gene share of ||d||^2
        r = {"cohort": c, "source": s, "offset_over_realA_spread": float(np.linalg.norm(d) / spread)}
        for a, b in BANDS: r[f"share_rank_{a}-{b}"] = float(share[(rank >= a) & (rank <= b)].sum())
        r["top10_genes"] = ", ".join(np.array(G[c].gene)[HV[c]][np.argsort(-share)[:10]]); arows.append(r)
A = pd.DataFrame(arows); A.to_csv(f"{RES}/pca_offset_attribution.csv", index=False); A.round(3)''')
md("## 5. Causal check: Harmony + iLISI on HVG subsets\nSame V0 pipeline (scale by real-A mean / SD, PCA 30 fit on real A, Harmony batch = source, iLISI perplexity 30, k = 90, median of 10,000 cells), 20,000 + 20,000 cells. HVG subsets: all 2,000; only those within the top 1,000 real-ranked genes; only those beyond rank 1,000. If K=4,000 mixes well on the top-1,000 HVGs, the poor integration comes from the genes it under-generates.")
code(r'''def ilisi(c, s, cols, seed=3):
    rng = np.random.default_rng(seed); ZR, ZS = SC[(c, "ref")], SC[(c, s)]
    ia = rng.choice(len(ZR), 20000, replace=False); ib = rng.choice(len(ZS), 20000, replace=False)
    A_, B_ = ZR[ia][:, cols], ZS[ib][:, cols]; pca = PCA(30, random_state=0).fit(A_)
    P = np.vstack([pca.transform(A_), pca.transform(B_)]).astype(np.float64); src = np.r_[np.zeros(20000), np.ones(20000)].astype(int).astype(str)
    ho = harmonypy.run_harmony(P, pd.DataFrame({"s": src}), ["s"], random_state=0, max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.Z_corr); Z = Z.T if Z.shape[0] == 30 else Z
    return L.lisi_stats(P, src, seed=0)[0], L.lisi_stats(Z, src, seed=0)[0]
irows = []
for c in COH:
    rank = G[c].real_rank.to_numpy()[HV[c]]
    for name, cols in (("all 2,000 HVGs", np.arange(len(rank))), ("HVGs within real top 1,000", np.where(rank <= 1000)[0]), ("HVGs beyond real rank 1,000", np.where(rank > 1000)[0])):
        for s in SRC[1:]:
            pre, post = ilisi(c, s, cols); irows.append({"cohort": c, "hvg_set": name, "n_genes": len(cols), "source": s, "iLISI_pre": pre, "iLISI_post": post})
            print(c, name, s, round(pre, 3), round(post, 3), flush=True)
I = pd.DataFrame(irows); I.to_csv(f"{RES}/ilisi_by_hvg_subset.csv", index=False)
I.pivot_table(index=["cohort", "hvg_set", "n_genes"], columns="source", values="iLISI_post").round(3)''')
code(r'''rows = []
for c in COH:
    t = G[c]
    for lo, hi in BANDS[:4]:
        for hv_ in (True, False):
            d = t[(t.real_rank >= lo) & (t.real_rank <= hi) & (t.is_hvg == hv_)]
            rows.append({"cohort": c, "real_rank_band": f"{lo:,}–{hi:,}", "genes": "HVG" if hv_ else "non-HVG", "n": len(d),
                         **{f"{s}_mean_umi_ratio": float(np.median(d[f"{s}_mean_umi"] / d.ref_mean_umi)) for s in SRC[1:]},
                         **{f"{s}_det_ratio": float(np.median(d[f"{s}_det"] / d.ref_det)) for s in SRC[1:]}})
H = pd.DataFrame(rows); H.to_csv(f"{RES}/hvg_vs_nonhvg_same_rank.csv", index=False); H.round(2)''')
md("## Figure")
code(r'''def T(ax, s): ax.set_title(s, loc="left", fontsize=11, fontweight="bold", color=INK, pad=6)
def figure(c):
    t = G[c]; h = t[t.is_hvg].sort_values("real_rank"); fl = 2e-5
    fig = plt.figure(figsize=(16, 10.5)); gs = fig.add_gridspec(3, 3, height_ratios=[1, 1.05, 1.0], width_ratios=[1, 1.25, 1.0], hspace=0.62, wspace=0.32,
                                                                 left=0.06, right=0.985, top=0.9, bottom=0.07)
    name = {"sepsis": "Sepsis — Ye 2026", "covid": "COVID-19 — Imperial ESKD"}[c]
    fig.text(0.06, 0.955, f"{name}: the 2,000 HVGs Harmony uses — real vs PBISC K=4,000 vs PBISC floor", fontsize=17, fontweight="bold", color=INK)
    fig.text(0.06, 0.928, "HVGs = seurat_v3 on real half A (as in the integration pipeline) · ~60,000 cells per source · real rank = order of real-A mean UMI over the 6,000 panel genes",
             fontsize=9.5, color=INK2)
    # A: where HVGs sit
    ax = fig.add_subplot(gs[0, 0]); bins = np.arange(0, 6001, 250)
    ax.hist(h.real_rank, bins=bins, color="#4A4F5C", edgecolor="white"); ax.axvline(1000, color=COL["prev4k"], ls="--", lw=1.2)
    f1 = (h.real_rank <= 1000).mean(); ax.text(1080, ax.get_ylim()[1] * 0.92, f"{100 * f1:.0f}% of HVGs in real top 1,000\n{100 * (1 - f1):.0f}% beyond rank 1,000",
                                               fontsize=9, color=INK, va="top")
    ax.set_xlabel("real rank among 6,000 panel genes (left = highest mean)"); ax.set_ylabel("HVGs per 250 ranks"); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    T(ax, "A  Where the 2,000 HVGs sit")
    # B: HVG gene-count curve
    ax = fig.add_subplot(gs[0, 1:]); x = h.real_rank.to_numpy()
    for s in ("prev4k", "floor", "realU"):
        y = h[f"{s}_mean_umi"].clip(lower=fl).to_numpy(); ax.scatter(x, y, s=2, color=COL[s], alpha=0.35, lw=0, rasterized=True)
        ax.plot(x, 10 ** pd.Series(np.log10(y)).rolling(61, center=True, min_periods=10).median(), color=COL[s], lw=1.8, label=LAB[s])
    ax.plot(x, h.ref_mean_umi.clip(lower=fl), color=COL["ref"], lw=2.2, label=LAB["ref"]); ax.axvline(1000, color=MUTED, ls=":", lw=1)
    ax.set_yscale("log"); ax.set_ylim(fl * 0.8, 40); ax.set_xlim(0, 6000); ax.set_yticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10], ["0.0001", "0.001", "0.01", "0.1", "1", "10"])
    ax.set_xlabel("HVGs placed at their real rank (x = real rank among 6,000 panel genes)"); ax.set_ylabel("mean UMI per cell"); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left", ncol=4); T(ax, "B  Gene-count curve, HVGs only (dots = each HVG; line = running median over 61 HVGs)")
    # C: detection ratio
    ax = fig.add_subplot(gs[1, 0])
    for s in ("realU", "prev4k", "floor"):
        r = (h[f"{s}_det"] / h.ref_det).clip(1e-3, 10); ax.scatter(x, r, s=2, color=COL[s], alpha=0.3, lw=0, rasterized=True)
        ax.plot(x, 10 ** pd.Series(np.log10(r.to_numpy())).rolling(61, center=True, min_periods=10).median(), color=COL[s], lw=1.8)
    ax.axhline(1, color=INK, lw=0.8); ax.set_yscale("log"); ax.set_ylim(5e-3, 5); ax.set_xlim(0, 6000)
    ax.set_xlabel("real rank"); ax.set_ylabel("detection fraction / real A"); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True); T(ax, "C  Detection, HVGs only")
    # D: z-shift (what Harmony sees)
    ax = fig.add_subplot(gs[1, 1])
    for s in ("realU", "prev4k", "floor"):
        ax.scatter(x, h[f"{s}_z_mean"], s=2, color=COL[s], alpha=0.3, lw=0, rasterized=True)
        ax.plot(x, pd.Series(h[f"{s}_z_mean"].to_numpy()).rolling(61, center=True, min_periods=10).median(), color=COL[s], lw=1.8, label=LAB[s])
    ax.axhline(0, color=INK, lw=0.8); ax.set_ylim(-1.2, 1.2); ax.set_xlim(0, 6000)
    ax.set_xlabel("real rank"); ax.set_ylabel("mean scaled value (real-A SD units)"); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    T(ax, "D  Input to PCA: per-HVG shift vs real A"); ax.legend(frameon=False, fontsize=8, loc="upper right")
    # E: SD ratio
    ax = fig.add_subplot(gs[1, 2])
    for s in ("realU", "prev4k", "floor"):
        r = (h[f"{s}_z_sd"] / h.ref_z_sd).clip(1e-2, 5); ax.scatter(x, r, s=2, color=COL[s], alpha=0.3, lw=0, rasterized=True)
        ax.plot(x, pd.Series(r.to_numpy()).rolling(61, center=True, min_periods=10).median(), color=COL[s], lw=1.8)
    ax.axhline(1, color=INK, lw=0.8); ax.set_ylim(0, 2); ax.set_xlim(0, 6000)
    ax.set_xlabel("real rank"); ax.set_ylabel("SD of scaled value / real A"); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True); T(ax, "E  Input to PCA: per-HVG spread")
    # F: offset attribution
    ax = fig.add_subplot(gs[2, 0]); a = A[A.cohort == c].set_index("source"); cols = [f"share_rank_{lo}-{hi}" for lo, hi in BANDS]
    labs = [f"{lo:,}–{hi:,}" for lo, hi in BANDS]; xx = np.arange(len(BANDS))
    for j, s in enumerate(("prev4k", "floor")):
        ax.bar(xx + (j - 0.5) * 0.38, a.loc[s, cols].to_numpy(), 0.36, color=COL[s], label=f"{LAB[s]} (offset {a.loc[s, 'offset_over_realA_spread']:.2f}× spread)")
    nh = [((h.real_rank >= lo) & (h.real_rank <= hi)).sum() for lo, hi in BANDS]
    ax.set_xticks(xx, [f"{l}\n({n} HVGs)" for l, n in zip(labs, nh)], fontsize=8); ax.set_ylabel("share of PCA centroid offset"); ax.axhline(0, color=INK, lw=0.8)
    ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=8, loc="upper right"); ax.set_xlabel("real rank band")
    ax.set_ylim(0, ax.get_ylim()[1] * 1.35); T(ax, "F  Which HVGs push generated cells away (PCA 30)")
    # G: causal iLISI
    ax = fig.add_subplot(gs[2, 1]); ii = I[I.cohort == c]; sets = ["all 2,000 HVGs", "HVGs within real top 1,000", "HVGs beyond real rank 1,000"]
    for j, s in enumerate(("realU", "prev4k", "floor")):
        v = [ii[(ii.hvg_set == st) & (ii.source == s)].iLISI_post.item() for st in sets]
        ax.bar(np.arange(3) + (j - 1) * 0.27, v, 0.25, color=COL[s], label=LAB[s])
        for k, vv in enumerate(v): ax.text(k + (j - 1) * 0.27, vv + 0.01, f"{vv:.2f}", ha="center", fontsize=7.5, color=INK)
    n_ = [ii[ii.hvg_set == st].n_genes.iloc[0] for st in sets]
    short = {"all 2,000 HVGs": "all HVGs", "HVGs within real top 1,000": "HVGs in real top 1,000", "HVGs beyond real rank 1,000": "HVGs beyond rank 1,000"}
    ax.set_xticks(range(3), [f"{short[st]}\n(n = {n:,})" for st, n in zip(sets, n_)], fontsize=8); ax.set_ylim(1, 2.15); ax.set_ylabel("iLISI after Harmony (2 = perfect)")
    ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=8, loc="upper center", ncol=3); T(ax, "G  Causal check: Harmony on HVG subsets")
    # H: HVG vs non-HVG at the same real rank
    ax = fig.add_subplot(gs[2, 2]); bb = BANDS[:4]; xx = np.arange(len(bb))
    for j, (s, hv_, alpha, lab) in enumerate((("prev4k", True, 1.0, "K=4,000, HVGs"), ("prev4k", False, 0.35, "K=4,000, non-HVGs"),
                                               ("floor", True, 1.0, "floor, HVGs"), ("floor", False, 0.35, "floor, non-HVGs"))):
        v = []
        for lo, hi in bb:
            d = t[(t.real_rank >= lo) & (t.real_rank <= hi) & (t.is_hvg == hv_)]; v.append(np.median(d[f"{s}_mean_umi"] / d.ref_mean_umi))
        ax.bar(xx + (j - 1.5) * 0.2, v, 0.19, color=COL[s], alpha=alpha, label=lab)
    ax.axhline(1, color=INK, lw=0.8); ax.set_xticks(xx, [f"{lo:,}–{hi:,}" for lo, hi in bb], fontsize=8); ax.set_xlabel("real rank band")
    ax.set_ylabel("mean UMI / real A (median gene)"); ax.set_ylim(0, 1.9); ax.yaxis.grid(True, color=GRID); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper right"); T(ax, "H  Same rank, HVG vs non-HVG")
    for ext in ("png", "pdf", "svg"): fig.savefig(f"{FIG}/hvg_harmony_genes_{c}.{ext}", dpi=250)
    plt.show()
figure("sepsis")''')
code('figure("covid")')
md("## Save the HVG statistics as h5ad (genes × statistics) for checking")
code(r'''for c in COH:
    t = G[c].set_index("gene"); num = t.select_dtypes("number")
    a = ad.AnnData(X=num.to_numpy(np.float32).T, obs=pd.DataFrame(index=num.columns), var=t[["is_hvg", "real_rank"]].astype({"is_hvg": str}))
    a.uns["note"] = f"{c}: rows = per-gene statistics (mean UMI, detection, scaled mean / SD on HVGs) for real A, real B, PBISC K=4,000, PBISC floor; columns = 6,000 panel genes."
    a.write_h5ad(f"{OUT}/hvg_gene_stats_{c}.h5ad")
sorted(os.listdir(OUT))''')
nb = nbf.v4.new_notebook(); nb.cells = C; nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
nbf.write(nb, "/home/claude/work/hvg/harmony_hvg_diagnosis.ipynb"); print("ok")
