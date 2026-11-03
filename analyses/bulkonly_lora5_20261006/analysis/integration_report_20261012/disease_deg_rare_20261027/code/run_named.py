import numpy as np, pandas as pd, pickle
allde=pickle.load(open("out/allde.pkl","rb"))
PANEL={
 ("sepsis","Monocyte"):["IL1R2","MCEMP1","AQP9","TNFAIP6","ACSL1","S100A8","S100A9","HLA-DRA","HLA-DRB1","HLA-DPA1","CD74","IFI27","C1QA","C1QB","C1QC","VSIG4","CLU","RNASE1","PLD4","CXCL1","CCL20","ADAMTS2","RETN","SERPINB2"],
 ("covid","Monocyte"):["IFI27","ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","OAS1","SIGLEC1","LY6E","HLA-DRA","CD74","S100A8","S100A9"],
 ("covid","T cell"):["IFI27","ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","OAS1","LY6E","S100A8","S100A9"],
 ("covid","NK cell"):["IFI27","ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","OAS1","LY6E","S100A8","S100A9"],
 ("covid","B cell"):["IFI27","ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","OAS1","LY6E","S100A8","S100A9"],
 ("covid","DC"):["IFI27","ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","OAS1","LY6E","S100A8","S100A9"],
 ("sepsis","T cell"):["NR4A1","NR4A2","NR4A3","ATF3","EGR1","SOCS3","DUSP1","CDKN1A","HBB","SLC4A10","ZNF683","S100A8","S100A9"],
 ("sepsis","NK cell"):["NR4A1","ATF3","SOCS3","EGR1","HBB","GNLY","NKG7","PRF1","S100A8","S100A9"],
 ("sepsis","B cell"):["NR4A1","SOCS3","TNF","MZB1","JCHAIN","S100A8","S100A9"],
 ("sepsis","DC"):["FCER1A","LILRA4","GZMB","PTGDS","CLIC3","CD1C","IRF7","IL1R2"],
}
rows=[]
for (coh,ct),genes in PANEL.items():
    if (coh,ct) not in allde: continue
    tabs=allde[(coh,ct)]
    R=tabs["real"].set_index("gene")
    for s in ["floor","prev"]:
        if s not in tabs: continue
        F=tabs[s].set_index("gene")
        for g in genes:
            if g not in R.index or g not in F.index:
                rows.append(dict(cohort=coh,celltype=ct,source=s,gene=g,note="not in panel")); continue
            rows.append(dict(cohort=coh,celltype=ct,source=s,gene=g,
                real_logFC=round(float(R.logFC[g]),2), real_q=float(R.q[g]), real_DEG=bool(R.deg[g]),
                gen_logFC=round(float(F.logFC[g]),2), gen_q=float(F.q[g]), gen_DEG=bool(F.deg[g]),
                ratio=round(float(F.logFC[g]/R.logFC[g]),2) if abs(R.logFC[g])>0.2 else np.nan,
                sign_ok=bool(np.sign(R.logFC[g])==np.sign(F.logFC[g]))))
N=pd.DataFrame(rows); N.to_csv("out/named_genes.csv",index=False)
pd.set_option("display.width",260); pd.set_option("display.max_rows",400)
f=N[N.source=="floor"]
for coh,ct in [("sepsis","Monocyte"),("sepsis","DC"),("covid","Monocyte"),("covid","T cell")]:
    d=f[(f.cohort==coh)&(f.celltype==ct)]
    print(f"\n### {coh} / {ct}  (real vs floor)")
    print(d[["gene","real_logFC","real_DEG","gen_logFC","gen_DEG","ratio","sign_ok"]].to_string(index=False))
print("\n=== IFN panel magnitude ratio (floor/real), COVID, by cell type ===")
ifn=["IFI27","ISG15","IFI6","IFI44L","IFIT1","IFIT3","MX1","OAS1","LY6E"]
d=f[(f.cohort=="covid")&(f.gene.isin(ifn))]
print(d.pivot_table(index="gene",columns="celltype",values="ratio").round(2).to_string())
print("\n=== S100A8/A9 leakage into lymphoid lineages (COVID) ===")
d=f[(f.cohort=="covid")&(f.gene.isin(["S100A8","S100A9"]))]
print(d[["celltype","gene","real_logFC","real_DEG","gen_logFC","gen_DEG"]].to_string(index=False))
