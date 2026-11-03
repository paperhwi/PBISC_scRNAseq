import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
IN="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/data/big/"
TR="/mnt/user-data/uploads/integration_report_20261012/ilisi_floor_20261021/results/big/"
st=lambda m,*p: np.any([np.char.startswith(np.asarray(m,dtype=str),q) for q in p],axis=0)
SEPSIS={"Mono_c06-CD14-MKI67":lambda m,l:m=="Monocytes","B_c08-MKI67":lambda m,l:np.isin(m,["B cells","Plasma cells"]),
 "T_CD4_c07-MKI67":lambda m,l:m=="CD4+ T cells","T_CD8_c06-MKI67":lambda m,l:m=="CD8+ T cells","NK_c04-MKI67":lambda m,l:m=="NK cells",
 "B_c07-MZB1":lambda m,l:np.isin(m,["B cells","Plasma cells"]),"pDC":lambda m,l:m=="DC","cDC":lambda m,l:m=="DC",
 "Macrophage":lambda m,l:np.isin(m,["Monocytes","Macrophages"]),"gdT":lambda m,l:np.isin(m,["CD4+ T cells","CD8+ T cells","NKT cells"])|(l=="gdT"),
 "ILC":lambda m,l:(m=="NK cells")|(l=="ILC"),"HSC":lambda m,l:np.ones(len(m),bool),
 "Mono_c01-CD14-HLA-DRA":lambda m,l:m=="Monocytes","NKT_c05-IFNG-CCL4hi":lambda m,l:m=="NKT cells"}
PAR_S={"Mono_c06-CD14-MKI67":"Monocytes","B_c08-MKI67":"B + plasma","T_CD4_c07-MKI67":"CD4 T","T_CD8_c06-MKI67":"CD8 T","NK_c04-MKI67":"NK",
 "B_c07-MZB1":"B + plasma","pDC":"DC","cDC":"DC","Macrophage":"Mono + macrophage","gdT":"all T","ILC":"NK + ILC","HSC":"all cells",
 "Mono_c01-CD14-HLA-DRA":"Monocytes","NKT_c05-IFNG-CCL4hi":"NKT"}
COVID={"CD14mono_dex":lambda m,l:st(m,"CD14mono"),"CD14mono_IFN":lambda m,l:st(m,"CD14mono"),"B_ASC_IgG":lambda m,l:st(m,"B_ASC"),
 "B_ASC_dividing":lambda m,l:st(m,"B_ASC"),"T CD4 cycling":lambda m,l:st(m,"T CD4","T reg"),"ASDC":lambda m,l:st(m,"DC","ASDC","pDC"),
 "DC1":lambda m,l:st(m,"DC","ASDC","pDC"),"pDC":lambda m,l:st(m,"DC","ASDC","pDC"),"HSPC":lambda m,l:np.ones(len(m),bool),
 "ILC":lambda m,l:st(m,"NK","ILC"),"T g/d":lambda m,l:st(m,"T ","MAIT","NKT"),"MAIT":lambda m,l:st(m,"T ","MAIT","NKT"),
 "NK CD56 bright":lambda m,l:st(m,"NK")}
PAR_C={"CD14mono_dex":"CD14 mono","CD14mono_IFN":"CD14 mono","B_ASC_IgG":"plasmablasts","B_ASC_dividing":"plasmablasts","T CD4 cycling":"CD4 T",
 "ASDC":"DC","DC1":"DC","pDC":"DC","HSPC":"all cells","ILC":"NK + ILC","T g/d":"all T","MAIT":"all T","NK CD56 bright":"NK"}
def lognorm(X):
    X=sp.csr_matrix(X,copy=True).astype(np.float64)
    t=np.asarray(X.sum(1)).ravel(); t[t==0]=1
    X=sp.diags(1e4/t)@X; X.data=np.log1p(X.data); return sp.csr_matrix(X)
def load_coh(coh):
    R=ad.read_h5ad(f"{IN}{coh}_ref.h5ad"); U=ad.read_h5ad(f"{IN}{coh}_realU.h5ad"); F=ad.read_h5ad(f"{IN}{coh}_floor.h5ad")
    tr={s:pd.read_csv(f"{TR}{coh}_{s}_transfer.csv.gz",index_col=0,keep_default_na=False) for s in ("realU","floor")}
    genes=np.array(R.var_names)
    X={"refA":lognorm(R.X),"realB":lognorm(U.X),"floor":lognorm(F.X)}
    V={"refA":(R.obs.major.astype(str).to_numpy(),R.obs.fine.astype(str).to_numpy()),
       "realB_auth":(tr["realU"].major.to_numpy(),tr["realU"].fine.to_numpy()),
       "realB_tr":(tr["realU"].major_tr.to_numpy(),tr["realU"].fine_tr.to_numpy()),
       "floor_tr":(tr["floor"].major_tr.to_numpy(),tr["floor"].fine_tr.to_numpy())}
    SRC={"refA":"refA","realB_auth":"realB","realB_tr":"realB","floor_tr":"floor"}
    return X,V,SRC,genes,(SEPSIS,PAR_S) if coh=="sepsis" else (COVID,PAR_C)
