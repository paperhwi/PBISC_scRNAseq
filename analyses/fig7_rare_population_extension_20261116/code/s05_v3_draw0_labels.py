"""Label all cells of an independently trained adapter draw (final model, draw 0; 656,000 cells) for adapter robustness.
Reference = draw-1 cells with the official 30-NN transfer labels (annotated h5ad). Classifier: 15-NN majority vote in a
50-PC space (log CP10k, 2,000 HVGs chosen on draw 1; PCA fitted on 120,000 draw-1 cells). Checks: (i) held-out draw-1 cells
vs official labels; (ii) draw-0 per-sample fine composition of the published 60k subset vs the official draw-0 60k labels."""
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, warnings; warnings.filterwarnings('ignore')
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import KNeighborsClassifier
UP='/mnt/user-data/uploads/PBISC/'; OUT='/home/claude/run/cache/'
rng=np.random.default_rng(20261115)
def lognorm(X):
    X=sp.csr_matrix(X,dtype=np.float32); s=np.asarray(X.sum(1)).ravel(); s[s==0]=1
    X=sp.diags(1e4/s)@X; X.data=np.log1p(X.data); return X
b=ad.read_h5ad(UP+'deliverables/annotated_h5ad_20261112/pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad',backed='r')
lab=b.obs.celltype_fine.astype(str).values; n=b.n_obs
idx=np.sort(rng.choice(n,160000,replace=False)); ref=idx[:120000]; ho=idx[120000:]
Xb=lognorm(b[idx].to_memory().X);
mu=np.asarray(Xb.mean(0)).ravel(); var=np.asarray(Xb.multiply(Xb).mean(0)).ravel()-mu**2
hv=np.argsort(-(var/(mu+1e-3)))[:2000]
svd=TruncatedSVD(50,random_state=0).fit(Xb[:120000][:,hv])
Zr=svd.transform(Xb[:120000][:,hv]); Zh=svd.transform(Xb[120000:][:,hv])
knn=KNeighborsClassifier(15,n_jobs=2).fit(Zr,lab[idx[:120000]])
ph=knn.predict(Zh); acc=(ph==lab[idx[120000:]]).mean()
rows=[]
for c in ['NK CD56 bright','NK','NK IFN stim','B_ASC_dividing','HSPC','DC1','DC2','pDC','CD14mono_IFN']:
    t=lab[idx[120000:]]==c; p=ph==c; rows.append(dict(label=c,n_true=int(t.sum()),precision=(t&p).sum()/max(p.sum(),1),recall=(t&p).sum()/max(t.sum(),1)))
V=pd.DataFrame(rows); V.insert(0,'overall_fine_accuracy',acc); print(V.round(3).to_string())
a=ad.read_h5ad(UP+'NewLora/imperial_eskd_draw0_K4000_ceiling.h5ad',backed='r')
pred=[]
for s in range(0,a.n_obs,50000):
    X=lognorm(a[s:s+50000].to_memory().X); pred.append(knn.predict(svd.transform(X[:,hv]))); print(s,flush=True)
pred=np.concatenate(pred)
o=pd.DataFrame({'sample_id':a.obs.sample_id.astype(str).values,'fine':pred}); o.to_parquet(OUT+'draw0_knn_labels.parquet')
# check against official draw-0 60k composition
L0=pd.read_csv(UP+'NewLora/analysis_20261102/covid/ceiling_cell_labels.csv.gz')
c_off=pd.crosstab(L0['sample'],L0.fine_tr,normalize='index'); c_new=pd.crosstab(o.sample_id,o.fine,normalize='index')
chk=[]
for c in ['NK CD56 bright','NK','B_ASC_dividing','HSPC','CD14mono','T CD8 CTL']:
    x=c_off.get(c).reindex(c_new.index); y=c_new.get(c); chk.append(dict(label=c,sample_r=np.corrcoef(x.fillna(0),y.fillna(0))[0,1],mean_official=x.mean(),mean_knn=y.mean()))
C=pd.DataFrame(chk); print(C.round(4).to_string())
V.to_csv('data_v3/draw0_label_validation_heldout.csv',index=False); C.to_csv('data_v3/draw0_label_validation_vs_official60k.csv',index=False)
