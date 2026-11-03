"""Figures for the lab-meeting deck. White background, Arial (Liberation Sans is the metric-identical fallback)."""
import sys, json, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import logging; logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
R = "/home/claude/work/report"; U = "/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; F = f"{R}/figs"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "axes.spines.top": False,
                     "axes.spines.right": False, "legend.frameon": False, "savefig.dpi": 220})
DRAW = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3"]
SRC = {"real": "#7f7f7f", "k1000": "#2b83ba", "k4000": "#d7191c"}
SRCN = {"real": "Real (public)", "k1000": "K = 1,000", "k4000": "K = 4,000"}
COH = {"sepsis": "Sepsis (Ye 2026)", "covid": "ESKD COVID (Stephenson 2025)"}
AB = {"T cell": "T", "Monocyte": "Mono", "B cell": "B", "NK cell": "NK", "Other": "Other", "DC": "DC"}
TAB = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)


def scat(ax, U2, lab, title, order=None, colors=None, s=0.3, legend=True, label_on=False, ncol=1, fs=6.5):
    lab = np.asarray(lab).astype(str); cats = order or sorted(pd.unique(lab))
    colors = colors or {c: TAB[i % len(TAB)] for i, c in enumerate(cats)}
    idx = np.random.default_rng(0).permutation(len(U2))
    ax.scatter(U2[idx, 0], U2[idx, 1], c=[colors.get(x, "#cccccc") for x in lab[idx]], s=s, linewidths=0, rasterized=True)
    if label_on:
        for c in cats:
            m = lab == c
            if m.sum() > 100: ax.text(*np.median(U2[m], 0), c, fontsize=fs, ha="center", va="center", weight="bold")
    if legend:
        h = [plt.Line2D([], [], marker="o", ls="", ms=4, color=colors[c], label=c) for c in cats if c in set(lab)]
        ax.legend(handles=h, loc="center left", bbox_to_anchor=(1, 0.5), fontsize=fs, ncol=ncol, handletextpad=0.1, borderaxespad=0)
    ax.set_title(title); ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(False) for sp in ax.spines.values()]


def gen_umaps(coh):
    fig, axs = plt.subplots(2, 4, figsize=(15, 6.6))
    l1 = ["T cell", "Monocyte", "B cell", "NK cell", "DC", "Other"]; l1c = dict(zip(l1, ["#4C72B0", "#C44E52", "#55A868", "#8172B3", "#DD8452", "#999999"]))
    for i, K in enumerate(["k1000", "k4000"]):
        name = f"{coh}_{K}"; rows = np.load(f"{U}/{name}_counts.npz")["rows"]; o = pd.read_csv(f"{U}/{name}_obs.csv.gz").iloc[rows]
        pre, post = np.load(f"{R}/{name}_umap_pre.npy"), np.load(f"{R}/{name}_umap_post.npy")
        lei = np.load(f"{R}/{name}_leiden_r0.5.npy")[rows].astype(str)
        dc = dict(zip([f"draw{j}" for j in range(5)], DRAW))
        scat(axs[i, 0], pre, o.draw, f"{SRCN[K]} — before Harmony (draw)", colors=dc, order=list(dc), legend=(i == 0) and False)
        scat(axs[i, 1], post, o.draw, f"{SRCN[K]} — after Harmony (draw)", colors=dc, order=list(dc), legend=True)
        scat(axs[i, 2], post, o.celltype_l1_pred, "after Harmony — predicted L1", colors=l1c, order=l1, legend=True)
        lo = [str(x) for x in sorted(set(lei.astype(int)))]
        scat(axs[i, 3], post, lei, "after Harmony — Leiden r = 0.5", order=lo, legend=False, label_on=True, fs=7)
    fig.suptitle(f"{COH[coh]}: generated cells, Harmony with batch = draw (100 cells / sample / draw shown)", x=0.01, ha="left", fontsize=11)
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_gen_umap.png"); plt.close(fig)


def batch_bars(A):
    names = ["sepsis_k1000", "sepsis_k4000", "covid_k1000", "covid_k4000"]; names = [n for n in names if n in A]
    lab = [n.replace("sepsis_", "Sepsis\n").replace("covid_", "COVID\n").replace("k1000", "K=1,000").replace("k4000", "K=4,000") for n in names]
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6)); x = np.arange(len(names)); w = 0.36
    ax = axs[0]
    ax.bar(x - w / 2, [A[n]["kBET_pre"] for n in names], w, color="#bdbdbd", label="before Harmony")
    ax.bar(x + w / 2, [A[n]["kBET_post"] for n in names], w, color="#2b83ba", label="after Harmony")
    ax.plot([x[0] - 0.5, x[-1] + 0.5], [np.mean([A[n]["kBET_null_post"] for n in names])] * 2, ls="--", c="k", lw=0.8, label="permutation null")
    for xi, n in zip(x, names): ax.text(xi + w / 2, A[n]["kBET_post"] + 0.01, f"{A[n]['kBET_post']:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x, lab); ax.set_ylabel("kBET rejection rate (k = 30)\nlower = better draw mixing"); ax.set_ylim(0, 1); ax.legend(fontsize=8, loc="upper right")
    ax.set_title("kBET: residual draw effect")
    ax = axs[1]
    ax.bar(x - w / 2, [A[n]["iLISI_draw_pre"] for n in names], w, color="#bdbdbd", label="before Harmony")
    ax.bar(x + w / 2, [A[n]["iLISI_draw_post"] for n in names], w, color="#2b83ba", label="after Harmony")
    ax.axhline(5, ls="--", c="k", lw=0.8); ax.text(x[-1] + 0.5, 5.05, "max = 5 draws", ha="right", fontsize=8)
    for xi, n in zip(x, names): ax.text(xi + w / 2, A[n]["iLISI_draw_post"] + 0.05, f"{A[n]['iLISI_draw_post']:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x, lab); ax.set_ylabel("median iLISI (draw)\nhigher = better mixing"); ax.set_ylim(0, 5.5); ax.set_title("iLISI: local draw diversity")
    plt.tight_layout(); fig.savefig(f"{F}/batch_metrics.png"); plt.close(fig)


def clust_bars(A):
    names = [n for n in ["sepsis_k1000", "sepsis_k4000", "covid_k1000", "covid_k4000"] if n in A]
    lab = [n.replace("sepsis_", "Sepsis ").replace("covid_", "COVID ").replace("k1000", "K1k").replace("k4000", "K4k") for n in names]
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.5)); x = np.arange(len(names)); w = 0.2
    for ax, met in zip(axs[:2], ["ARI", "NMI"]):
        for j, (key, c, l) in enumerate([("l1_r0.5", "#2b83ba", "L1, r0.5"), ("l2_r0.5", "#abd9e9", "L2, r0.5"), ("l1_r1.0", "#d7191c", "L1, r1.0"), ("l2_r1.0", "#fdae61", "L2, r1.0")]):
            ax.bar(x + (j - 1.5) * w, [A[n][f"{met}_{key}"] for n in names], w, color=c, label=l)
        ax.set_xticks(x, lab); ax.set_ylabel(f"{met} vs reference labels"); ax.set_title(f"{met}: Leiden clusters vs cell-type labels"); ax.set_ylim(0, 0.75)
    axs[0].legend(fontsize=7, ncol=2, loc="upper left")
    ax = axs[2]
    ax.bar(x - 0.18, [A[n]["clusters_r0.5"] for n in names], 0.36, color="#2b83ba", label="r = 0.5")
    ax.bar(x + 0.18, [A[n]["clusters_r1.0"] for n in names], 0.36, color="#d7191c", label="r = 1.0")
    for xi, n in zip(x, names):
        for dx, r in [(-0.18, "0.5"), (0.18, "1.0")]:
            ax.text(xi + dx, A[n][f"clusters_r{r}"] + 0.4, f"{A[n][f'draw_dominated_clusters_r{r}']}", ha="center", fontsize=7, color="k")
    ax.set_xticks(x, lab); ax.set_ylabel("Leiden clusters"); ax.set_ylim(0, 30); ax.legend(fontsize=8, ncol=2, loc="upper left"); ax.set_title("clusters (number above bar = draw-dominated)")
    plt.tight_layout(); fig.savefig(f"{F}/clustering_metrics.png"); plt.close(fig)


def cluster_draw_comp(coh):
    fig, axs = plt.subplots(1, 2, figsize=(13, 3.6), sharey=True)
    for ax, K in zip(axs, ["k1000", "k4000"]):
        o = pd.read_csv(f"{U}/{coh}_{K}_obs.csv.gz"); lab = np.load(f"{R}/{coh}_{K}_leiden_r0.5.npy")
        ct = pd.crosstab(lab, o.draw, normalize="index"); l1 = pd.crosstab(lab, o.celltype_l1_pred, normalize="index").idxmax(1)
        size = pd.Series(lab).value_counts().reindex(ct.index)
        bottom = np.zeros(len(ct))
        for j, d in enumerate(ct.columns):
            ax.bar(np.arange(len(ct)), ct[d], bottom=bottom, color=DRAW[j], label=d, width=0.85); bottom += ct[d].to_numpy()
        ax.set_xticks(np.arange(len(ct)), [f"{c}\n{AB.get(l1[c], l1[c])}\n{size[c]/1e3:.0f}k" for c in ct.index], fontsize=7)
        ax.axhline(0.2, ls=":", c="k", lw=0.6); ax.set_title(f"{SRCN[K]}: draw share per Leiden cluster (r = 0.5)")
    axs[0].set_ylabel("fraction of cluster"); axs[1].legend(fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_cluster_draw.png"); plt.close(fig)


def joint_umap(coh):
    z = np.load(f"{R}/{coh}_joint_umap.npz", allow_pickle=True); hp = os.path.exists(f"{R}/{coh}_joint_umap_pre.npy")
    src = z["src"].astype(str); maj = z["major"].astype(str)
    fig, axs = plt.subplots(1, 4 if hp else 3, figsize=(18, 4.4) if hp else (14, 4.4))
    if not hp: axs = [None] + list(axs)
    if hp: scat(axs[0], np.load(f"{R}/{coh}_joint_umap_pre.npy"), np.array([SRCN[s] for s in src]), "PCA projection (before Harmony): source", colors={SRCN[k]: v for k, v in SRC.items()}, order=list(SRCN.values()), legend=False)
    scat(axs[1], z["umap"], np.array([SRCN[s] for s in src]), "Harmony (batch = source): source", colors={SRCN[k]: v for k, v in SRC.items()}, order=list(SRCN.values()))
    cats = sorted(pd.unique(maj)); cc = {c: TAB[i % 20] for i, c in enumerate(cats)}
    m = src == "real"; scat(axs[2], z["umap"][m], maj[m], "real cells: published labels (major)", colors=cc, order=cats, legend=False, label_on=True, fs=6)
    m = src != "real"; scat(axs[3], z["umap"][m], maj[m], "generated cells: transferred labels", colors=cc, order=cats, ncol=1, fs=6)
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_joint_umap.png"); plt.close(fig)


def mixing_bars():
    fig, ax = plt.subplots(figsize=(5.6, 3.4)); x = 0; ticks = []
    for coh in ["sepsis", "covid"]:
        J = json.load(open(f"{R}/{coh}_joint.json"))["mix"]
        for K in ["k1000", "k4000"]:
            ax.bar(x - 0.18, J[f"iLISI_real_vs_{K}_pre"], 0.36, color="#bdbdbd", label="PCA projection" if x == 0 else None)
            ax.bar(x + 0.18, J[f"iLISI_real_vs_{K}_post"], 0.36, color=SRC[K], label=None)
            ax.text(x + 0.18, J[f"iLISI_real_vs_{K}_post"] + 0.02, f"{J[f'iLISI_real_vs_{K}_post']:.2f}", ha="center", fontsize=8)
            ticks.append(f"{coh.title() if coh == 'sepsis' else 'COVID'}\nreal vs {SRCN[K]}"); x += 1
    ax.axhline(2, ls="--", c="k", lw=0.8); ax.text(3.5, 2.03, "perfect mixing = 2", ha="right", fontsize=8)
    ax.axhline(1, ls=":", c="k", lw=0.6)
    ax.set_xticks(range(4), ticks, fontsize=8); ax.set_ylim(0, 2.6); ax.set_ylabel("median iLISI (real vs generated)")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#bdbdbd"), plt.Rectangle((0, 0), 1, 1, color="#555555")], labels=["PCA projection", "after Harmony"], fontsize=8, loc="upper left", ncol=2)
    ax.set_title("Generated vs real cells barely co-mix"); plt.tight_layout(); fig.savefig(f"{F}/joint_mixing.png"); plt.close(fig)


def lineage(x):
    x = str(x)
    if x.startswith("B_ASC"): return "Plasmablasts"
    if x.startswith("B_"): return "B cells"
    if x.startswith("CD14mono"): return "CD14 monocytes"
    if x.startswith(("CD16mono", "Int.mono")): return "CD16 / int. monocytes"
    if x.startswith(("DC", "ASDC", "pDC")): return "DC"
    if x.startswith("NK ") or x == "NK": return "NK cells"
    if x.startswith(("T CD4", "T reg")): return "CD4+ T cells"
    if x.startswith("T CD8"): return "CD8+ T cells"
    if x in ("MAIT", "T g/d", "NKT", "ILC"): return "MAIT / gdT / NKT / ILC"
    if x in ("Platelet", "RBC"): return "Platelet / RBC"
    return "HSPC" if x == "HSPC" else x


def composition(coh, grps):
    T = {s: pd.read_csv(f"{R}/{coh}_comp_{s}_major.csv", index_col=0) for s in ["real", "k1000", "k4000"]}
    if coh == "covid": T = {s: v.groupby(v.index.map(lineage)).sum() for s, v in T.items()}
    cats = list(T["real"].index); cc = {c: TAB[i % 20] for i, c in enumerate(sorted(cats))}
    fig, ax = plt.subplots(figsize=(8.5, 3.8)); x = 0; ticks = []; pos = []
    for g in grps:
        for s in ["real", "k1000", "k4000"]:
            b = 0
            for c in cats:
                v = T[s].loc[c, g] if c in T[s].index else 0; ax.bar(x, v, bottom=b, color=cc[c], width=0.8); b += v
            ticks.append({"real": "real", "k1000": "K1k", "k4000": "K4k"}[s]); pos.append(x); x += 1
        ax.text(x - 2, 1.03, g.replace("_", " "), ha="center", fontsize=9, weight="bold"); x += 0.6
    ax.set_xticks(pos, ticks, fontsize=8); ax.set_ylim(0, 1.1); ax.set_ylabel("fraction of PBMC")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=cc[c]) for c in cats], labels=cats, fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_title(f"{COH[coh]}: lineage composition (generated cells labelled by kNN transfer)", loc="left")
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_composition.png"); plt.close(fig)


def draw_comp(coh, case):
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
    real = pd.read_csv(f"{R}/{coh}_comp_real_major.csv", index_col=0)[case]
    if coh == "covid": real = real.groupby(real.index.map(lineage)).sum()
    cats = list(real.index); cc = {c: TAB[i % 20] for i, c in enumerate(sorted(cats))}
    for ax, K in zip(axs, ["k1000", "k4000"]):
        t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False); t = t[t.group == case]
        if coh == "covid": t = t.assign(major=t.major.map(lineage))
        ct = pd.crosstab(t.major, t.draw, normalize="columns").reindex(cats).fillna(0); ct.insert(0, "real", real)
        b = np.zeros(ct.shape[1])
        for c in cats: ax.bar(range(ct.shape[1]), ct.loc[c], bottom=b, color=cc[c], width=0.8); b += ct.loc[c].to_numpy()
        ax.set_xticks(range(ct.shape[1]), ct.columns, fontsize=8); ax.set_title(f"{SRCN[K]}: {case.replace('_', ' ')} by draw")
    axs[0].set_ylabel("fraction"); axs[1].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=cc[c]) for c in cats], labels=cats, fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_draw_comp.png"); plt.close(fig)


def shift(coh, C):
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.2), sharex=True, sharey=True)
    for ax, K in zip(axs, ["k1000", "k4000"]):
        d = pd.read_csv(f"{R}/{coh}_{K}_shift.csv", index_col=0); x, y = d.real_log2FC, d[f"{K}_log2FC"]
        ax.scatter(x, y, s=18, color=SRC[K], alpha=0.8, linewidths=0)
        lim = max(np.abs(x).max(), np.abs(y).max()) * 1.08; ax.plot([-lim, lim], [-lim, lim], ls="--", c="grey", lw=0.7)
        ax.axhline(0, c="k", lw=0.4); ax.axvline(0, c="k", lw=0.4)
        top = (np.abs(x) + np.abs(y)).sort_values(ascending=False).index[:9]
        for i in top: ax.text(x[i], y[i], i.split("-")[0] if len(i) > 14 else i, fontsize=6.5)
        ax.set_title(f"{SRCN[K]}: Spearman ρ = {C[f'{K}_label_shift_spearman']:.2f}, sign agreement {C[f'{K}_label_shift_sign_agree']:.0%} (n = {C[f'{K}_label_shift_n']})", fontsize=9)
        ax.set_xlabel("real: log2 (case / control) subset fraction")
    axs[0].set_ylabel("generated: log2 (case / control)")
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_shift.png"); plt.close(fig)


def rare_table(coh, rare, case, ctrl):
    real = pd.read_csv(f"{R}/{coh}_real_labels.csv.gz", keep_default_na=False)
    if coh == "sepsis": real = real[real.grp.str.startswith("adult")]
    rows = []
    for lab in rare:
        r = {"subset": lab, "real_case_%": 100 * (real[real.grp == case].label == lab).mean(), "real_ctrl_%": 100 * (real[real.grp == ctrl].label == lab).mean()}
        for K in ["k1000", "k4000"]:
            t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False)
            r[f"{K}_case_%"] = 100 * (t[t.group == case].label == lab).mean(); r[f"{K}_ctrl_%"] = 100 * (t[t.group == ctrl].label == lab).mean()
            r[f"{K}_conf"] = t[t.label == lab].conf.median() if (t.label == lab).any() else np.nan
        rows.append(r)
    df = pd.DataFrame(rows); df.to_csv(f"{R}/{coh}_rare.csv", index=False); return df


RARE = {"sepsis": (["Mono_c06-CD14-MKI67", "B_c08-MKI67", "T_CD4_c07-MKI67", "T_CD8_c06-MKI67", "NK_c04-MKI67", "T_CD4_c02-GPR183-NR4A2",
                     "B_c07-MZB1", "pDC", "cDC", "Macrophage", "HSC", "ILC", "Megakaryocyte", "gdT"], "adult_sepsis", "adult_HC"),
        "covid": (["CD14mono_dex", "CD14mono_IFN", "B_ASC_IgG", "B_ASC_dividing", "T CD4 cycling", "ASDC", "DC1", "pDC", "HSPC", "ILC",
                   "T g/d", "MAIT", "NK CD56 bright", "Platelet"], "covid_positive", "covid_negative")}


def rare_fig(coh):
    rare, case, ctrl = RARE[coh]; df = rare_table(coh, rare, case, ctrl)
    fig, ax = plt.subplots(figsize=(11, 3.6)); x = np.arange(len(df)); w = 0.27; fl = 0.001
    for j, (s, c) in enumerate([("real", SRC["real"]), ("k1000", SRC["k1000"]), ("k4000", SRC["k4000"])]):
        v = df[f"{s}_case_%"].clip(lower=fl); ax.bar(x + (j - 1) * w, v, w, color=c, label=SRCN[s])
        z = df[f"{s}_case_%"] == 0
        for xi in x[z.to_numpy()]: ax.text(xi + (j - 1) * w, fl * 1.2, "0", ha="center", fontsize=6.5, color=c)
    ax.set_yscale("log"); ax.set_ylim(fl, 30); ax.set_ylabel(f"% of cells in {case.replace('_', ' ')}")
    ax.set_xticks(x, df.subset, rotation=35, ha="right", fontsize=8); ax.legend(fontsize=8, ncol=3, loc="upper right")
    ax.set_title(f"{COH[coh]}: rare / disease-associated subsets after label transfer", loc="left")
    plt.tight_layout(); fig.savefig(f"{F}/{coh}_rare.png"); plt.close(fig)


if __name__ == "__main__":
    what = sys.argv[1:]
    A = json.load(open(f"{R}/metrics_A.json")) if os.path.exists(f"{R}/metrics_A.json") else {}
    for coh in [c for c in ["sepsis", "covid"] if c in what]:
        if os.path.exists(f"{R}/{coh}_k4000_umap_post.npy"): gen_umaps(coh); cluster_draw_comp(coh)
        joint_umap(coh)
        rare_fig(coh)
        C = json.load(open(f"{R}/{coh}_composition.json"))
        if coh == "sepsis":
            composition(coh, ["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]); draw_comp(coh, "adult_sepsis"); shift(coh, C)
        else:
            composition(coh, ["covid_negative", "covid_positive"]); draw_comp(coh, "covid_positive"); shift(coh, C)
    if "summary" in what:
        batch_bars(A); clust_bars(A); mixing_bars()


def real_ref():
    import h5py
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.6), gridspec_kw={"width_ratios": [1.25, 1, 1, 1]})
    so = pd.concat([pd.read_csv(f"/mnt/user-data/outputs/sepsis_ann/sepsis_obs_part{i}.csv.gz", keep_default_na=False, usecols=["celltype_major_ye2026", "in_paper_scheme"]) for i in range(3)], ignore_index=True)
    Us = np.load("/home/claude/work/sepsis/run_major/umap.npy"); keep = so.in_paper_scheme.astype(str).eq("True").to_numpy()
    sub = np.where(keep)[0]; sub = np.random.default_rng(0).permutation(sub)[:250000]
    scat(axs[0], Us[sub], so.celltype_major_ye2026.to_numpy()[sub], f"Sepsis real: {keep.sum():,} cells, 133 samples (Harmony by sample)", s=0.15, legend=False, label_on=True, fs=7)
    L = pd.read_csv("/mnt/user-data/outputs/covid_ann/covid_obs_recap.csv.gz", index_col=0, keep_default_na=False, usecols=[0, 1, 2] if False else None)
    lab = L.cell_type_detailed_recap.to_numpy()
    for ax, f, t in zip(axs[1:], ["myeloid", "bcell", "tnk"], ["COVID real: myeloid", "COVID real: B / plasma", "COVID real: T / NK"]):
        with h5py.File(f"/home/claude/work/covid/{f}.h5ad") as h: Uc = h["obsm/X_umap"][:]; idx = h["obs/cell_index"][:]
        sub = np.random.default_rng(0).permutation(len(Uc))[:150000]
        scat(ax, Uc[sub], lab[idx][sub], f"{t} ({len(Uc):,} cells)", s=0.2, legend=False, label_on=True, fs=6)
    plt.tight_layout(); fig.savefig(f"{F}/real_ref_umap.png"); plt.close(fig)
