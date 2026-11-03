"""Figures for the two calibration checks (real-vs-real baseline; real cells through a PBISC-like read-out)."""
import json, os, logging, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
O = "/home/claude/work/report/check"; F = "/home/claude/work/report/figs"
plt.rcParams.update({"font.family": ["Arial", "Liberation Sans"], "font.size": 9, "axes.titlesize": 10, "figure.facecolor": "white",
                     "axes.facecolor": "white", "savefig.facecolor": "white", "axes.spines.top": False, "axes.spines.right": False,
                     "legend.frameon": False, "savefig.dpi": 220})
COH = {"sepsis": "Sepsis (reference: Ye 2026 cells)", "covid": "ESKD COVID (reference: Stephenson 2025 cells)"}
def order(coh):
    oth = "covid" if coh == "sepsis" else "sepsis"
    return [("realU", "Real, held-out samples (same study)", "#6a51a3"), ("realT1", "Real, depth-matched", "#9e9ac8"),
            ("realT2", "Real, PBISC-like read-out", "#e69f00"), (oth, f"Real, other study ({'COVID' if oth == 'covid' else 'sepsis'})", "#5aa469"),
            ("k1000", "Generated K = 1,000", "#2b83ba"), ("k4000", "Generated K = 4,000", "#d7191c")]


def mixing():
    cohs = [c for c in ["sepsis", "covid"] if os.path.exists(f"{O}/{c}_checks.json")]
    fig, axs = plt.subplots(1, len(cohs), figsize=(6.2 * len(cohs), 3.6), sharex=True)
    axs = np.atleast_1d(axs)
    for ax, coh in zip(axs, cohs):
        r = json.load(open(f"{O}/{coh}_checks.json")); od = order(coh)[::-1]; y = np.arange(len(od))
        pre = [r[k]["iLISI_pre"] for k, _, _ in od]; post = [r[k]["iLISI_post"] for k, _, _ in od]
        ax.barh(y + 0.19, pre, 0.36, color="#e3e3e3", label="before Harmony (projection)")
        ax.barh(y - 0.19, post, 0.36, color=[c for _, _, c in od], label="after Harmony")
        for yi, v in zip(y, post): ax.text(v + 0.02, yi - 0.19, f"{v:.2f}", va="center", fontsize=8)
        for yi, v in zip(y, pre): ax.text(v + 0.02, yi + 0.19, f"{v:.2f}", va="center", fontsize=7, color="#7f7f7f")
        ax.set_yticks(y, [l for _, l, _ in od]); ax.axvline(1, c="k", lw=0.5, ls=":"); ax.axvline(2, c="k", lw=0.8, ls="--")
        ax.set_xlim(0.8, 2.15); ax.set_xlabel("median iLISI vs reference (1 = separate, 2 = fully mixed)"); ax.set_title(COH[coh], loc="left")
    axs[0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#e3e3e3"), plt.Rectangle((0, 0), 1, 1, color="#6a51a3")],
                  labels=["before Harmony (projection)", "after Harmony (pairwise)"], fontsize=8, loc="lower right")
    plt.tight_layout(); fig.savefig(f"{F}/calib_mixing.png"); plt.close(fig)


def lineage():
    cohs = [c for c in ["sepsis", "covid"] if os.path.exists(f"{O}/{c}_checks.json")]
    fig, axs = plt.subplots(1, len(cohs), figsize=(6.2 * len(cohs), 3.5), sharey=True); axs = np.atleast_1d(axs)
    L = ["T", "NK", "B", "Myeloid"]
    for ax, coh in zip(axs, cohs):
        r = json.load(open(f"{O}/{coh}_checks.json")); od = [o for o in order(coh) if o[0] != "realT1"]; w = 0.8 / len(od)
        for j, (k, lab, c) in enumerate(od):
            ax.bar(np.arange(4) + (j - (len(od) - 1) / 2) * w, [r[k].get(f"iLISI_post_{l}", np.nan) for l in L], w, color=c, label=lab)
        ax.set_xticks(range(4), ["T cells", "NK cells", "B cells", "Myeloid"]); ax.set_ylim(1, 2); ax.axhline(2, c="k", lw=0.8, ls="--")
        ax.set_title(COH[coh], loc="left"); ax.set_ylabel("median iLISI within lineage, after Harmony")
    h, l = axs[0].get_legend_handles_labels(); l = [x.replace(" (COVID)", "") for x in l]
    fig.legend(h, l, loc="lower center", ncol=5, fontsize=8, bbox_to_anchor=(0.5, 0.0))
    plt.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(f"{F}/calib_lineage.png"); plt.close(fig)


def umaps():
    cohs = [c for c in ["sepsis", "covid"] if os.path.exists(f"{O}/{c}_umap_k1000.npz")]
    fig, axs = plt.subplots(len(cohs), 4, figsize=(16, 4.0 * len(cohs))); axs = np.atleast_2d(axs)
    for i, coh in enumerate(cohs):
        oth = "covid" if coh == "sepsis" else "sepsis"; lab = dict((k, (l, c)) for k, l, c in order(coh))
        for ax, k in zip(axs[i], ["realU", "realT2", oth, "k1000"]):
            z = np.load(f"{O}/{coh}_umap_{k}.npz", allow_pickle=True); U2, s = z["umap"], z["src"]
            p = np.random.default_rng(0).permutation(len(U2))
            ax.scatter(U2[p, 0], U2[p, 1], c=np.where(s[p] == 0, "#cfcfcf", lab[k][1]), s=0.3, linewidths=0, rasterized=True)
            ax.set_title(f"{'Sepsis' if coh == 'sepsis' else 'COVID'} reference (grey) vs\n{lab[k][0]}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([]); [sp.set_visible(False) for sp in ax.spines.values()]
    plt.tight_layout(); fig.savefig(f"{F}/calib_umap.png"); plt.close(fig)


def table():
    rows = []
    for coh in ["sepsis", "covid"]:
        if not os.path.exists(f"{O}/{coh}_checks.json"): continue
        r = json.load(open(f"{O}/{coh}_checks.json"))
        for k, lab, _ in [("ref", "Reference (real)", "")] + order(coh):
            st = r["stats"][k]; d = {"cohort": coh, "source": lab, **{f"{a}": st[a] for a in st}}
            if k != "ref":
                d.update({x: r[k][x] for x in r[k] if not x.startswith("iLISI_pre_") and not x.startswith("iLISI_post_")})
            rows.append(d)
    df = pd.DataFrame(rows); df.to_csv(f"{O}/calibration_summary.csv", index=False); return df


if __name__ == "__main__":
    mixing(); lineage(); umaps(); print(table().round(3).to_string())
