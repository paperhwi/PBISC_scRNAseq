import numpy as np, pandas as pd, pickle
from scipy import stats
from deg_lib import *
allde=pickle.load(open("out/allde.pkl","rb"))
rows=[]; tops=[]
for (coh,ct),tabs in allde.items():
    R=tabs["real"]; rexp=(R.mean_case+R.mean_ctrl)/2
    rank=stats.rankdata(-rexp); hi=rank>2000
    realDEG=set(R.gene[R.deg])
    for s in ["floor","prev"]:
        if s not in tabs: continue
        F=tabs[s].set_index("gene").loc[R.gene].reset_index()
        r=dict(cohort=coh,celltype=ct,source=s)
        r["sign_rank2000plus"]=float((np.sign(R.logFC[hi])==np.sign(F.logFC[hi])).mean())
        r["slope_rank2000plus"]=float(np.polyfit(R.logFC[hi],F.logFC[hi],1)[0])
        # precision of generated top-100 DEGs (|logFC| ranked among q<0.05) vs real DEG set
        t=F[F.q<0.05]
        top=t.reindex(t.logFC.abs().sort_values(ascending=False).index).gene.head(100).tolist()
        r["gen_top100_n"]=len(top)
        r["gen_top100_precision"]=float(np.mean([g in realDEG for g in top])) if top else np.nan
        # same-direction among those
        lr=R.set_index("gene").logFC; lf=F.set_index("gene").logFC
        r["gen_top100_signmatch"]=float(np.mean([np.sign(lf[g])==np.sign(lr[g]) for g in top])) if top else np.nan
        rows.append(r)
    # top DEG lists
    def tl(T,k=30):
        t=T[T.q<0.05]
        t=t.reindex(t.logFC.abs().sort_values(ascending=False).index).head(k)
        return list(zip(t.gene,t.logFC.round(2)))
    tops.append(dict(cohort=coh,celltype=ct,
                     real=tl(R), floor=tl(tabs["floor"]) if "floor" in tabs else [],
                     prev=tl(tabs["prev"]) if "prev" in tabs else []))
S=pd.DataFrame(rows); S.to_csv("out/deg_rank2000_summary.csv",index=False)
pd.set_option("display.width",250)
print(S.round(3).to_string(index=False))
print()
print("=== aggregate (median over 5 six-types; and excluding DC) ===")
D=pd.read_csv("out/deg_summary.csv")
for coh in ["sepsis","covid"]:
    for s in ["prev","floor"]:
        d=D[(D.cohort==coh)&(D.source==s)]; dn=d[d.celltype!="DC"]
        q=S[(S.cohort==coh)&(S.source==s)]; qn=q[q.celltype!="DC"]
        print(f"{coh:7s} {s:5s} logFC_r med {d.logFC_r.median():.2f} / exDC mean {dn.logFC_r.mean():.2f} | "
              f"rank2000+ r exDC {dn.logFC_r_rank2000plus.mean():.2f} | slope exDC {qn.slope_rank2000plus.mean():.2f} | "
              f"sign rank2000+ exDC {qn.sign_rank2000plus.mean():.2f} | gen top100 prec exDC {qn.gen_top100_precision.mean():.2f}")
pickle.dump(tops,open("out/tops.pkl","wb"))
for t in tops:
    print(f"\n--- {t['cohort']} / {t['celltype']}")
    print("  real :", ", ".join(f"{g}({v:+.1f})" for g,v in t["real"][:15]))
    print("  floor:", ", ".join(f"{g}({v:+.1f})" for g,v in t["floor"][:15]))
