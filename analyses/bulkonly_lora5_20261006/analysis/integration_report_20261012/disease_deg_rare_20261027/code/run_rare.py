import numpy as np, pandas as pd, pickle
from scipy import stats
from rare_lib import *
rng=np.random.default_rng(20261027)
MIN=30
rows=[]; genediff=[]
for coh in ["sepsis","covid"]:
    X,V,SRC,genes,(DEF,PAR)=load_coh(coh)
    mA,fA=V["refA"]
    # expressed-gene filter from real A
    detA=np.asarray((X["refA"]>0).mean(0)).ravel()
    gok=detA>=0.01
    for sub,fn in DEF.items():
        parA=fn(mA,fA); tgA=parA&(fA==sub)
        if tgA.sum()<MIN: continue
        A=X["refA"]
        # top-10 markers from real A (subset vs rest of parent)
        d=np.asarray(A[tgA].mean(0)).ravel()-np.asarray(A[parA&~tgA].mean(0)).ravel()
        det=np.asarray((A[tgA]>0).mean(0)).ravel(); d2=d.copy(); d2[det<0.1]=-np.inf
        top=np.argsort(-d2)[:10]; mk=genes[top]
        prof={"refA":np.asarray(A[tgA].mean(0)).ravel()}
        ncell={"refA":int(tgA.sum())}
        for lab in ["realB_tr","floor_tr"]:
            m,f=V[lab]; p=fn(m,f); t=p&(f==sub)
            ncell[lab]=int(t.sum())
            prof[lab]=np.asarray(X[SRC[lab]][t].mean(0)).ravel() if t.sum()>0 else None
        if prof["floor_tr"] is None or ncell["floor_tr"]<MIN:
            rows.append(dict(cohort=coh,subset=sub,parent=PAR[sub],n_refA=ncell["refA"],n_realB_tr=ncell["realB_tr"],n_floor_tr=ncell["floor_tr"],note="floor n<30")); continue
        r=dict(cohort=coh,subset=sub,parent=PAR[sub],n_refA=ncell["refA"],n_realB_tr=ncell["realB_tr"],n_floor_tr=ncell["floor_tr"],markers=",".join(mk))
        a=prof["refA"][gok]
        for lab,nm in [("realB_tr","realB"),("floor_tr","floor")]:
            if prof[lab] is None or ncell[lab]<10: continue
            b=prof[lab][gok]
            r[f"r_{nm}"]=float(np.corrcoef(a,b)[0,1])
            r[f"rho_{nm}"]=float(stats.spearmanr(a,b).statistic)
            diff=b-a
            r[f"nDiff_{nm}"]=int((np.abs(diff)>=0.5).sum())
            r[f"meanAbsDiff_{nm}"]=float(np.abs(diff).mean())
        # marker score level within subset, and detection
        for lab,nm in [("refA","refA"),("realB_tr","realB"),("floor_tr","floor")]:
            m,f=(V[lab] if lab!="refA" else V["refA"]); p=fn(m,f); t=p&(f==sub) if lab!="refA" else tgA
            if t.sum()<10: continue
            Z=X[SRC.get(lab,"refA")] if lab!="refA" else X["refA"]
            r[f"mkScore_{nm}"]=float(np.asarray(Z[t][:,top].mean(1)).ravel().mean())
            r[f"mkDet_{nm}"]=float(np.asarray((Z[t][:,top]>0).mean(0)).ravel().mean())
            # identity AUROC within parent
            sc=np.asarray(Z[:,top].mean(1)).ravel()
            if 0<t[p].sum()<p.sum():
                from sklearn.metrics import roc_auc_score
                r[f"auroc_{nm}"]=float(roc_auc_score(t[p],sc[p]))
        rows.append(r)
        # per-gene differences for the biggest subsets
        if ncell["floor_tr"]>=100:
            dd=pd.DataFrame(dict(gene=genes[gok],refA=a,floor=prof["floor_tr"][gok],
                                 realB=prof["realB_tr"][gok] if prof["realB_tr"] is not None else np.nan))
            dd["d_floor"]=dd.floor-dd.refA; dd["d_realB"]=dd.realB-dd.refA
            dd["cohort"]=coh; dd["subset"]=sub
            genediff.append(dd)
R=pd.DataFrame(rows); R.to_csv("out/rare_expression.csv",index=False)
G=pd.concat(genediff,ignore_index=True) if genediff else pd.DataFrame()
G.to_csv("out/rare_gene_diff.csv.gz",index=False,compression="gzip")
pd.set_option("display.width",280); pd.set_option("display.max_columns",40)
cols=["cohort","subset","parent","n_refA","n_realB_tr","n_floor_tr","r_realB","r_floor","nDiff_realB","nDiff_floor","mkScore_refA","mkScore_realB","mkScore_floor","auroc_refA","auroc_realB","auroc_floor"]
print(R.reindex(columns=cols).round(3).to_string(index=False))
