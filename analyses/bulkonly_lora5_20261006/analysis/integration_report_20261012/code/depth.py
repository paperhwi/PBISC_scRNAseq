import sys, time, numpy as np, pandas as pd, scipy.sparse as sp
sys.path.insert(0,'/home/claude/work'); import store
U="/mnt/user-data/uploads/PBISC/.claude_pylib/extract"; t0=time.time()
for coh, sym_col in [("sepsis","gene_symbol"),("covid","gene")]:
    st=store.Store(f"/home/claude/work/{coh}/store"); sym=st.var[sym_col].astype(str).to_numpy()
    z=np.load(f"{U}/{coh}_k1000_counts.npz",allow_pickle=True); gvar=z["var"].astype(str)
    pos=pd.Series(np.arange(len(sym)),index=sym); pos=pos[~pos.index.duplicated()]
    gpos=pd.Series(np.arange(len(gvar)),index=gvar)
    panel=[g for g in gvar if g in pos.index]; gi=gpos.loc[panel].to_numpy(); ri=pos.loc[panel].to_numpy()
    print(coh,"setup",time.time()-t0,flush=True)
    rng=np.random.default_rng(0); idx=np.sort(rng.choice(st.n_cells,20000,replace=False))
    R=st.rows(idx); print("read",time.time()-t0,flush=True)
    Rp=R[:,ri].tocsr(); tot_all=np.asarray(R.sum(1)).ravel(); tot=np.asarray(Rp.sum(1)).ravel(); det=np.diff(Rp.indptr)
    for K in ["k1000","k4000"]:
        z=np.load(f"{U}/{coh}_{K}_counts.npz",allow_pickle=True)
        X=sp.csr_matrix((z["data"].astype(np.float32),z["indices"].astype(np.int32),z["indptr"]),shape=(len(z["rows"]),len(z["var"])))[:,gi].tocsr()
        gt=np.asarray(X.sum(1)).ravel(); gd=np.diff(X.indptr)
        print(coh,K,"gen panel UMI",np.percentile(gt,[25,50,75]).round(0),"detected",np.percentile(gd,[25,50,75]).round(0),"frac ones",round(float((X.data==1).mean()),3),flush=True)
    print(coh,"real panel UMI",np.percentile(tot,[25,50,75]).round(0),"all-gene UMI",np.percentile(tot_all,[25,50,75]).round(0),"detected",np.percentile(det,[25,50,75]).round(0),"frac ones",round(float((Rp.data==1).mean()),3),"panel share",round(float(np.median(tot/tot_all)),3),flush=True)
