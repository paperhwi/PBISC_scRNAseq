"""Gene-level mean-expression curves: generated vs real.

Definitions (PBISC notes 2026-10-20 / 2026-10-21 / Figure 6b-c):
* mean UMI and mean CP10k per gene, detection = fraction of cells with >= 1 UMI;
* genes are ranked by the *real* mean CP10k; the generated/real ratio is reported per rank band
  (default 1-500, 501-1,500, 1,501-2,000, 2,001-3,000, 3,001-4,000, >4,000) and per real-mean bin;
* correlation of log1p mean expression (Pearson, Spearman) over the shared genes;
* an optional depth-matched real control (binomial thinning to the generated depth).
"""
import numpy as np
import pandas as pd

from ._utils import align_genes, cp10k, depth, thin_to_depth, to_csr

DEFAULT_RANK_BANDS = (0, 500, 1500, 2000, 3000, 4000, np.inf)
DEFAULT_MEAN_BINS = (0, 0.03, 0.1, 0.2, 0.3, 0.4, 1, 3, 10, np.inf)


def gene_stats(X, genes):
    """Per-gene mean UMI, mean CP10k and detection of a cells x genes count matrix."""
    X = to_csr(X)
    n = X.shape[0]
    mean_umi = np.asarray(X.mean(0)).ravel()
    mean_cp10k = np.asarray(cp10k(X).mean(0)).ravel()
    det = np.asarray((X > 0).sum(0)).ravel() / n
    return pd.DataFrame({"gene": np.asarray(genes).astype(str), "mean_umi": mean_umi, "mean_cp10k": mean_cp10k, "detection": det})


def compare_gene_means(real_X, real_genes, gen_X, gen_genes, depth_match=True, seed=0):
    """Gene table with real and generated statistics side by side.

    depth_match=True also adds the 'real_dm_*' columns: real cells binomially thinned to the median generated depth, the control that
    separates a genuine under/over-generation from a pure depth difference.
    Returns a DataFrame sorted by real mean CP10k (rank 1 = most expressed) with columns rank, real_*, gen_*, ratio_cp10k.
    """
    ia, ib, names = align_genes(real_genes, gen_genes)
    R = to_csr(real_X)[:, ia]
    G = to_csr(gen_X)[:, ib]
    r, g = gene_stats(R, names), gene_stats(G, names)
    t = pd.DataFrame({"gene": names})
    for c in ("mean_umi", "mean_cp10k", "detection"):
        t[f"real_{c}"], t[f"gen_{c}"] = r[c].to_numpy(), g[c].to_numpy()
    if depth_match:
        Rt = thin_to_depth(R, float(np.median(depth(G))), seed=seed)
        d = gene_stats(Rt, names)
        for c in ("mean_umi", "mean_cp10k", "detection"):
            t[f"real_dm_{c}"] = d[c].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        t["ratio_cp10k"] = np.where(t.real_mean_cp10k > 0, t.gen_mean_cp10k / t.real_mean_cp10k, np.nan)
    t = t.sort_values("real_mean_cp10k", ascending=False, kind="mergesort").reset_index(drop=True)
    t.insert(1, "rank", np.arange(1, len(t) + 1))
    return t


def band_table(t, by="rank", edges=None):
    """Per rank-band (or real-mean-bin) summary of generated vs real.

    ratio = sum(generated CP10k) / sum(real CP10k) in the band; median_ratio = median per-gene ratio (genes with real > 0);
    det_ratio = sum(generated detection) / sum(real detection); never_generated = genes with zero generated UMI.
    """
    if by == "rank":
        edges = DEFAULT_RANK_BANDS if edges is None else edges
        key, lo_is_open = t["rank"].to_numpy(), True
    elif by == "mean":
        edges = DEFAULT_MEAN_BINS if edges is None else edges
        key, lo_is_open = t["real_mean_umi"].to_numpy(), False
    else:
        raise ValueError("by must be 'rank' or 'mean'")
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (key > lo) & (key <= hi) if lo_is_open else (key >= lo) & (key < hi)
        s = t[m]
        if len(s) == 0:
            continue
        pos = s.real_mean_cp10k > 0
        rows.append({
            "band": f"{lo:g}-{hi:g}", "n_genes": len(s),
            "real_mean_cp10k": s.real_mean_cp10k.mean(), "gen_mean_cp10k": s.gen_mean_cp10k.mean(),
            "ratio": s.gen_mean_cp10k.sum() / max(s.real_mean_cp10k.sum(), 1e-12),
            "median_ratio": float(np.median(s.gen_mean_cp10k[pos] / s.real_mean_cp10k[pos])) if pos.any() else np.nan,
            "det_ratio": s.gen_detection.sum() / max(s.real_detection.sum(), 1e-12),
            "never_generated": int((s.gen_mean_umi == 0).sum()),
            **({"dm_real_ratio": s.real_dm_mean_cp10k.sum() / max(s.real_mean_cp10k.sum(), 1e-12)} if "real_dm_mean_cp10k" in s else {}),
        })
    return pd.DataFrame(rows)


def correlation_summary(t, min_real_mean=0.0):
    """Pearson / Spearman r of log1p mean CP10k (generated vs real) and of detection, over genes with real mean > min_real_mean."""
    from scipy.stats import pearsonr, spearmanr
    s = t[t.real_mean_cp10k > min_real_mean]
    x, y = np.log1p(s.real_mean_cp10k), np.log1p(s.gen_mean_cp10k)
    out = {"n_genes": len(s), "pearson_log1p_cp10k": float(pearsonr(x, y)[0]), "spearman_cp10k": float(spearmanr(s.real_mean_cp10k, s.gen_mean_cp10k)[0]),
           "pearson_detection": float(pearsonr(s.real_detection, s.gen_detection)[0]),
           "genes_never_generated": int((t.gen_mean_umi == 0).sum()), "genes_never_in_real": int((t.real_mean_umi == 0).sum())}
    if "real_dm_mean_cp10k" in t:
        out["pearson_log1p_cp10k_depthmatched_real"] = float(pearsonr(np.log1p(s.real_mean_cp10k), np.log1p(s.real_dm_mean_cp10k))[0])
    return out


def plot_gene_mean_scatter(ax, t, label="generated", color="#c0392b"):
    """Figure 6b-c style: log10 mean UMI per gene, real (x) vs generated (y), with r in the title."""
    s = t[(t.real_mean_umi > 0) & (t.gen_mean_umi > 0)]
    ax.scatter(np.log10(s.real_mean_umi), np.log10(s.gen_mean_umi), s=3, alpha=0.35, color=color, linewidths=0, rasterized=True)
    lim = [min(np.log10(s.real_mean_umi).min(), np.log10(s.gen_mean_umi).min()), max(np.log10(s.real_mean_umi).max(), np.log10(s.gen_mean_umi).max())]
    ax.plot(lim, lim, color="k", lw=0.6, ls="--")
    r = np.corrcoef(np.log10(s.real_mean_umi), np.log10(s.gen_mean_umi))[0, 1]
    ax.set_xlabel("Real: log10 mean UMI / cell"); ax.set_ylabel(f"{label}: log10 mean UMI / cell"); ax.set_title(f"r = {r:.3f} (n = {len(s):,} genes)")
    return ax


def plot_gene_mean_curve(ax, t, label="generated", color="#c0392b", smooth=50):
    """Rank curve: log10 mean CP10k vs real rank for real (grey) and generated (colour; rolling median over `smooth` genes)."""
    x = t["rank"].to_numpy()
    eps = 1e-4
    ax.plot(x, np.log10(t.real_mean_cp10k + eps), color="0.4", lw=1.2, label="Real")
    if "real_dm_mean_cp10k" in t:
        ax.plot(x, np.log10(t.real_dm_mean_cp10k.rolling(smooth, center=True, min_periods=1).median() + eps), color="0.7", lw=1, ls="--", label="Real (depth-matched)")
    ax.plot(x, np.log10(t.gen_mean_cp10k.rolling(smooth, center=True, min_periods=1).median() + eps), color=color, lw=1.2, label=label)
    ax.set_xlabel("Gene rank by real mean CP10k"); ax.set_ylabel("log10 mean CP10k (+1e-4)"); ax.legend(frameon=False, fontsize=7)
    return ax
