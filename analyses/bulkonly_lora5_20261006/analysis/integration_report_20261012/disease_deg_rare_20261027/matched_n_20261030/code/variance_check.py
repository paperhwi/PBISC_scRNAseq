"""Is the generated between-sample variance deflated? Per-gene SD of unit-level log2CPM within each group, real vs floor."""
import numpy as np, pandas as pd
from deg_lib import *
rows=[]
for coh in ["sepsis","covid"]:
    srcs={s:load(f"{s}_{coh}") for s in ["real","floor"]}; genes=list(srcs["real"]["genes"])
    meas=[g for g,ok in zip(srcs["real"]["genes"],(srcs["real"]["S"].sum(0)>0)&(srcs["floor"]["S"].sum(0)>0)) if ok]
    for ct in ["T cell","B cell","NK cell","Monocyte"]:
        g={s:cpm_matrix(srcs[s],design(f"{s}_{coh}",srcs[s]),ct) for s in srcs}
        for grp in ["case","ctrl"]:
            sd={}; mu={}
            for s in g:
                x=g[s][g[s]._c==grp][meas]
                sd[s]=x.std(0); mu[s]=x.mean(0)
            ok=(mu["real"]>1)&(mu["floor"]>1)      # expressed genes only
            ratio=(sd["floor"][ok]/sd["real"][ok]).replace([np.inf,-np.inf],np.nan).dropna()
            rows.append(dict(cohort=coh,celltype=ct,group=grp,n_genes=int(ok.sum()),
                             median_SD_real=float(sd["real"][ok].median()),median_SD_floor=float(sd["floor"][ok].median()),
                             median_SD_ratio=float(ratio.median()),pct_genes_floorSD_lower=100*float((ratio<1).mean())))
V=pd.DataFrame(rows); V.to_csv("matched/variance_check.csv",index=False)
pd.set_option("display.width",200); print(V.round(3).to_string(index=False))
