"""D: composition of generated cells (labels transferred from real cells) vs the real reference, per clinical group."""
import sys, json, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr
coh = sys.argv[1]; R = "/home/claude/work/report"
real = pd.read_csv(f"{R}/{coh}_real_labels.csv.gz", keep_default_na=False)
if coh == "sepsis":
    real = real[real.grp.str.startswith("adult")]
    grps = ["adult_HC", "adult_ICU_ctrl", "adult_sepsis"]; case, ctrl = "adult_sepsis", "adult_HC"
else:
    grps = ["covid_negative", "covid_positive"]; case, ctrl = "covid_positive", "covid_negative"
def comp(df, col):
    return pd.crosstab(df[col], df.grp, normalize="columns").reindex(columns=grps).fillna(0)
out = {"groups": grps}
tables = {}
for col in ["major", "label"]:
    Cr = comp(real, col); tables[("real", col)] = Cr
    for K in ["k1000", "k4000"]:
        t = pd.read_csv(f"{R}/{coh}_{K}_transfer.csv.gz", keep_default_na=False).rename(columns={"group": "grp"})
        Cg = comp(t, col).reindex(Cr.index).fillna(0); tables[(K, col)] = Cg
        for g in grps:
            out[f"{K}_{col}_{g}_r"] = float(pearsonr(Cr[g], Cg[g])[0]); out[f"{K}_{col}_{g}_L1"] = float(np.abs(Cr[g] - Cg[g]).sum())
        # disease shift concordance (log2 case/control of subset proportions; subsets >= 0.2% in real)
        keep = Cr.max(1) >= 0.002; eps = 1e-3
        lr = np.log2((Cr.loc[keep, case] + eps) / (Cr.loc[keep, ctrl] + eps)); lg = np.log2((Cg.loc[keep, case] + eps) / (Cg.loc[keep, ctrl] + eps))
        out[f"{K}_{col}_shift_spearman"] = float(spearmanr(lr, lg)[0]); out[f"{K}_{col}_shift_sign_agree"] = float((np.sign(lr) == np.sign(lg)).mean())
        out[f"{K}_{col}_shift_n"] = int(keep.sum())
        if col == "label":
            pd.DataFrame({"real_log2FC": lr, f"{K}_log2FC": lg}).to_csv(f"{R}/{coh}_{K}_shift.csv")
        # per-draw variability of the case-group composition r
        rs = []
        for d, td in t.groupby("draw"):
            Cd = comp(td, col).reindex(Cr.index).fillna(0); rs.append(pearsonr(Cr[case], Cd[case])[0])
        out[f"{K}_{col}_{case}_r_per_draw"] = [round(float(x), 3) for x in rs]
        out[f"{K}_median_conf"] = float(t.conf.median())
for (s, col), T in tables.items():
    T.to_csv(f"{R}/{coh}_comp_{s}_{col}.csv")
json.dump(out, open(f"{R}/{coh}_composition.json", "w"), indent=1)
print(json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in out.items()}, indent=0))
