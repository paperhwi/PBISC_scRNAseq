import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, warnings, json
warnings.filterwarnings("ignore"); ad.settings.allow_write_nullable_strings=False
U='/mnt/user-data/uploads/PBISC/'; D=U+'eskd_covid_clinical_assoc_20261115/data/'
RARE13=['ASDC','B_ASC_IgG','B_ASC_dividing','CD14mono_IFN','CD14mono_dex','DC1','HSPC','ILC','MAIT','NK CD56 bright','T CD4 cycling','T g/d','pDC']
RARE_EXT=RARE13+['B_ASC_IgA','B_ASC_IgM','B_ASC_other']
lm=pd.read_csv(U+'deliverables/annotated_h5ad_20261112/tables/label_mapping_covid.csv').set_index('fine')
# ---------- Real ----------
r=pd.read_csv(D+'real_obs_full.csv.gz',low_memory=False)
r['fine']=r.cell_type_IFNstim_recap; r['broad6']=r.fine.map(lm.six_type_broad)
assert r.broad6.notna().all()
r[['_index','sample_id','patient_id','fine','broad6']].to_parquet('cache/real_cells.parquet')
z=np.load(D+'real_pseudobulk_sample_x_fine_panel.npz')
genes=z['genes']; grp=pd.Series(z['groups']).str.split('|',expand=True); grp.columns=['sample_id','fine']; grp['n']=z['ncells']
pbR=z['pb']
# rare cells
parts=[np.load(D+'real_rare_cells_counts_part%d.npz'%k,allow_pickle=True) for k in (0,1)]
var=pd.read_csv(D+'real_var.csv.gz')['_index'].astype(str).values
mats=[sp.csr_matrix((p['data'],p['indices'],p['indptr']),shape=(len(p['rows']),len(var))) for p in parts]
X=sp.vstack(mats).tocsr(); rows=np.concatenate([p['rows'] for p in parts])
assert (r['_index'].values[rows]==np.concatenate([p['obs_names'] for p in parts])).all()
keepc=['sample_id','patient_id','fine','broad6','case_control','WHO_severity','WHO_temp_severity','fatal_disease','centre','days_since_posswab','pool','percentMito','isg_score']
aR=ad.AnnData(X.astype(np.float32),obs=r.iloc[rows][['_index']+keepc].set_index('_index'),var=pd.DataFrame(index=var))
for c in aR.obs.columns:
    if aR.obs[c].dtype==object: aR.obs[c]=aR.obs[c].astype(str)
aR.var_names_make_unique()
# check: pseudobulk vs rare matrix (independent paths)
vi=pd.Series(np.arange(len(var)),index=var); pc=vi[genes].values
chk=[]
for (s,f),sub in aR.obs.groupby(['sample_id','fine'],observed=True):
    if len(chk)>=400: break
    i=np.where((grp.sample_id==s)&(grp.fine==f))[0][0]
    idx=aR.obs_names.get_indexer(sub.index)
    chk.append(np.abs(np.asarray(X[idx][:,pc].sum(0)).ravel()-pbR[i]).max())
print('pseudobulk check max abs diff over',len(chk),'groups:',max(chk))
assert grp.n.sum()==len(r)
aR.write_h5ad('cache/real_rare_cells.h5ad',compression='gzip')
np.savez_compressed('cache/real_pb.npz',pb=pbR,genes=genes); grp.to_csv('cache/real_pb_groups.csv',index=False)
# ---------- PBISC ----------
p=ad.read_h5ad(U+'deliverables/annotated_h5ad_20261112/pbisc_covid_eskd_newref_draw1_K4000_ceiling_annotated.h5ad',backed='r')
o=p.obs.copy(); o['fine']=o.celltype_fine.astype(str); o['broad6']=o.celltype_broad.astype(str)
o.to_parquet('cache/pbisc_obs.parquet')
gk=pd.Categorical(o.sample_id.astype(str)+'|'+o.fine); code=gk.codes; G=len(gk.categories)
pbP=np.zeros((G,p.n_vars)); rare=o.fine.isin(RARE_EXT).values; rm=[]
B=40000
for a in range(0,p.n_obs,B):
    b=min(a+B,p.n_obs); m=p.X[a:b]; m=sp.csr_matrix(m)
    M=sp.csr_matrix((np.ones(b-a),(code[a:b],np.arange(b-a))),shape=(G,b-a))
    pbP+=np.asarray((M@m).todense())
    if rare[a:b].any(): rm.append(m[rare[a:b]])
gp=pd.Series(gk.categories).str.split('|',expand=True); gp.columns=['sample_id','fine']; gp['n']=np.bincount(code,minlength=G)
np.savez_compressed('cache/pbisc_pb.npz',pb=pbP.astype(np.float32),genes=p.var.gene_symbol.astype(str).values,gene_id=p.var_names.values)
gp.to_csv('cache/pbisc_pb_groups.csv',index=False)
keepp=['sample_id','donor_id','fine','broad6','condition','WHO_severity','WHO_temp_severity','fatal_disease','time_from_first_positive_swab','transfer_confidence','in_primary_scope','wave']
ob=o[rare][keepp].copy()
for c in ob.columns:
    if str(ob[c].dtype) in('object','category'): ob[c]=ob[c].astype(str)
aP=ad.AnnData(sp.vstack(rm).tocsr().astype(np.float32),obs=ob,var=p.var.copy())
aP.write_h5ad('cache/pbisc_rare_cells.h5ad',compression='gzip')
print(aR, aP)
print(aP.obs.fine.value_counts()); print(aR.obs.fine.value_counts())
