"""Summaries of the panel-matched re-scoring: real vs generated, adult healthy control vs adult sepsis."""
import numpy as np, pandas as pd
O = "/home/claude/work/report/rescore"; rng = np.random.default_rng(20261012); NB = 2000
real = pd.read_csv(f"{O}/real_by_sample.csv"); gen = pd.read_csv(f"{O}/generated_by_sample_draw.csv")
cv = pd.read_csv(f"{O}/thresholds_cv_panel_matched.csv")
real = real[real.n_parent >= 20]; gen = gen[gen.n_parent >= 20]
gs = gen.groupby(["claim_id", "K", "sample", "group"])[["posA_frac", "posB_frac", "score_mean"]].mean().reset_index()   # average over 5 draws


def contrast(df, col, case="adult_sepsis", ctrl="adult_HC"):
    a = df[df.group == case][col].to_numpy(); b = df[df.group == ctrl][col].to_numpy()
    if len(a) < 2 or len(b) < 2: return np.nan, np.nan, np.nan, len(a), len(b)
    bs = [rng.choice(a, len(a)).mean() - rng.choice(b, len(b)).mean() for _ in range(NB)]
    return a.mean() - b.mean(), *np.percentile(bs, [2.5, 97.5]), len(a), len(b)


rows = []
for c in sorted(real.claim_id.unique()):
    r = real[real.claim_id == c]
    for v, col in [("A_panel", "posA_frac"), ("B_panel_depth", "posB_frac")]:
        row = {"claim_id": c, "variant": v, "cv_auc": float(cv[(cv.claim_id == c) & (cv.variant == v)].cv_auc.iloc[0])}
        for src, d, cc in [("real_label", r, "label_frac"), ("real_classifier", r, col), ("real_score", r, "score_mean")] + \
                          [(f"k{K//1000}k_classifier", gs[(gs.claim_id == c) & (gs.K == K)], col) for K in (1000, 4000)] + \
                          [(f"k{K//1000}k_score", gs[(gs.claim_id == c) & (gs.K == K)], "score_mean") for K in (1000, 4000)]:
            row[f"{src}_HC_median"] = d[d.group == "adult_HC"][cc].median(); row[f"{src}_sepsis_median"] = d[d.group == "adult_sepsis"][cc].median()
            m, lo, hi, na, nb = contrast(d, cc); row.update({f"{src}_diff": m, f"{src}_lo": lo, f"{src}_hi": hi, f"{src}_n": f"{na}/{nb}"})
        # paired K4 - K1 contrast (bootstrap samples jointly)
        w = gs[gs.claim_id == c].pivot_table(index=["sample", "group"], columns="K", values=col).dropna().reset_index()
        a = w[w.group == "adult_sepsis"]; b = w[w.group == "adult_HC"]
        dd = lambda A, B: (A[4000] - A[1000]).mean() - (B[4000] - B[1000]).mean()
        bs = [dd(a.sample(len(a), replace=True, random_state=int(rng.integers(1e9))), b.sample(len(b), replace=True, random_state=int(rng.integers(1e9)))) for _ in range(500)]
        row.update({"K4_minus_K1_diff": dd(a, b), "K4_minus_K1_lo": np.percentile(bs, 2.5), "K4_minus_K1_hi": np.percentile(bs, 97.5)})
        # per-draw contrasts
        for K in (1000, 4000):
            pdw = [contrast(gen[(gen.claim_id == c) & (gen.K == K) & (gen.draw == dr)], col)[0] for dr in range(5)]
            row[f"k{K//1000}k_per_draw_diff"] = ";".join(f"{x:+.3f}" for x in pdw)
        rows.append(row)
out = pd.DataFrame(rows); out.to_csv(f"{O}/summary_panel_matched.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 50)
A = out[out.variant == "A_panel"]
show = ["claim_id", "cv_auc", "real_label_HC_median", "real_label_sepsis_median", "real_classifier_HC_median", "real_classifier_sepsis_median",
        "real_classifier_diff", "real_classifier_lo", "real_classifier_hi", "k1k_classifier_HC_median", "k1k_classifier_sepsis_median", "k1k_classifier_diff",
        "k1k_classifier_lo", "k1k_classifier_hi", "k4k_classifier_HC_median", "k4k_classifier_sepsis_median", "k4k_classifier_diff", "K4_minus_K1_diff", "K4_minus_K1_lo", "K4_minus_K1_hi"]
print(A[show].round(3).T.to_string())
print(A[["claim_id", "real_score_diff", "real_score_lo", "real_score_hi", "k1k_score_diff", "k1k_score_lo", "k1k_score_hi", "k4k_score_diff", "k4k_score_lo", "k4k_score_hi",
         "real_score_HC_median", "k1k_score_HC_median", "k1k_per_draw_diff", "k4k_per_draw_diff"]].round(3).to_string())
B = out[out.variant == "B_panel_depth"]
print(B[["claim_id", "cv_auc", "real_classifier_HC_median", "real_classifier_sepsis_median", "k1k_classifier_HC_median", "k1k_classifier_sepsis_median", "k4k_classifier_HC_median", "k4k_classifier_sepsis_median"]].round(3).to_string())
