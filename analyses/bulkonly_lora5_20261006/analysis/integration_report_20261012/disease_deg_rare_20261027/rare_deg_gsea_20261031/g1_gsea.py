"""Preranked GSEA of the disease contrast per cell population, real vs PBISC floor (same pseudobulk DEG tables).
Rank = sign(log2FC) * -log10(p) (Welch). Gene sets: MSigDB Hallmark 2023.2 and Reactome 2023.2, restricted to the ~5,700-gene panel.
gseapy prerank, 1,000 gene-set permutations, size 10-500 in panel, seed 0."""
import numpy as np, pandas as pd, pickle, gseapy as gp, os, sys, warnings; warnings.filterwarnings("ignore")
A=pickle.load(open("out/allde.pkl","rb")); SIX=["T cell","B cell","NK cell","Monocyte","DC"]
GS={"hallmark":"heatgsea/../gsea/genesets/h.all.v2023.2.Hs.symbols.gmt","reactome":"gsea/genesets/c2.cp.reactome.gmt"}
out=[]
for lib,path in GS.items():
    for coh in ["sepsis","covid"]:
        for ct in SIX:
            for src in ["real","floor"]:
                t=A[(coh,ct)][src].dropna(subset=["p"]); p=np.clip(t.p.values,1e-300,1)
                rnk=pd.Series(np.sign(t.logFC.values)*-np.log10(p),index=t.gene.values).sort_values(ascending=False)
                rnk=rnk[~rnk.index.duplicated()]+np.random.default_rng(0).normal(0,1e-6,len(rnk))   # break ties
                res=gp.prerank(rnk=rnk,gene_sets=path,min_size=10,max_size=500,permutation_num=1000,seed=0,threads=2,outdir=None,verbose=False)
                r=res.res2d[["Term","ES","NES","NOM p-val","FDR q-val","Tag %","Lead_genes"]].copy()
                r.columns=["term","ES","NES","p","fdr","tag","lead_genes"]; r["library"]=lib; r["cohort"]=coh; r["celltype"]=ct; r["source"]=src
                out.append(r); print(lib,coh,ct,src,len(r),flush=True)
G=pd.concat(out,ignore_index=True)
for c in ["ES","NES","p","fdr"]: G[c]=pd.to_numeric(G[c],errors="coerce")
G.to_csv("heatgsea/out/gsea_all.csv.gz",index=False,compression="gzip")
