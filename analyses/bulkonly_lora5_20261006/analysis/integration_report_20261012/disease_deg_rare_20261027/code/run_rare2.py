import numpy as np, pandas as pd
R=pd.read_csv("out/rare_expression.csv")
R=R[R.r_floor.notna()].copy()
for nm in ["realB","floor"]:
    R[f"mkRatio_{nm}"]=R[f"mkScore_{nm}"]/R["mkScore_refA"]
    R[f"mkDetRatio_{nm}"]=R[f"mkDet_{nm}"]/R["mkDet_refA"]
R["nDiff_excess"]=R.nDiff_floor-R.nDiff_realB
pd.set_option("display.width",260)
print("=== marker-programme retention within the SAME subset label ===")
print(R[["cohort","subset","n_refA","n_realB_tr","n_floor_tr","mkScore_refA","mkRatio_realB","mkRatio_floor","mkDetRatio_realB","mkDetRatio_floor","r_realB","r_floor","nDiff_realB","nDiff_floor"]].round(3).to_string(index=False))
print()
print("median marker-score ratio: realB %.2f  floor %.2f"%(R.mkRatio_realB.median(),R.mkRatio_floor.median()))
print("median marker-detection ratio: realB %.2f  floor %.2f"%(R.mkDetRatio_realB.median(),R.mkDetRatio_floor.median()))
print("median profile r: realB %.4f  floor %.4f"%(R.r_realB.median(),R.r_floor.median()))
print("median n genes |delta|>=0.5: realB %.0f  floor %.0f"%(R.nDiff_realB.median(),R.nDiff_floor.median()))
R.to_csv("out/rare_expression_ratios.csv",index=False)

G=pd.read_csv("out/rare_gene_diff.csv.gz")
print("\n=== genes most depleted / enriched in floor within each rare subset (|d|, vs real A; real B shown as control) ===")
for (coh,sub),d in G.groupby(["cohort","subset"]):
    d=d.dropna(subset=["d_floor"])
    if len(d)==0: continue
    lo=d.nsmallest(8,"d_floor"); hi=d.nlargest(8,"d_floor")
    print(f"\n-- {coh} / {sub}")
    print("   down in floor:", ", ".join(f"{r.gene}({r.d_floor:+.2f}|B{r.d_realB:+.2f})" for r in lo.itertuples()))
    print("   up   in floor:", ", ".join(f"{r.gene}({r.d_floor:+.2f}|B{r.d_realB:+.2f})" for r in hi.itertuples()))
