"""COVID mapping validation with patient-disjoint, two-fold out-of-fold tests.

R2: reference PCA + kNN. R3: read-out-aware logistic regression (2,000 HVGs,
C=.05). R1: original Harmony+kNN route audited on real held-out controls in
fold 0 of each split protocol. Fixed generated cells are scored by R2/R3;
only folds that exclude their real-reference patient are retained.
"""
from common import *
import json,time,warnings
import scanpy as sc
import harmonypy
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import balanced_accuracy_score

start=time.time();rng=np.random.default_rng(7)
def log(*x):print(round(time.time()-start),*x,flush=True)
REAL=ROOT/'bulkonly_lora5_20261006/analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad'
with h5py.File(REAL,'r') as f:
    obs=ad.io.read_elem(f['obs']);real_genes=ad.io.read_elem(f['var']).index.astype(str)
ro=pd.DataFrame({'sample':obs.sample_id.astype(str).values,'patient':obs.patient_id.astype(str).values,
                 'major':obs.cell_type_detailed_recap.astype(str).values,'fine':obs.cell_type_IFNstim_recap.astype(str).values,
                 'group':obs.case_control.astype(str).str.lower().map({'positive':'covid_positive','negative':'covid_negative','recovery':'recovery'}).values})
ro['lin6']=ro.major.map(lin6)
roster=set(pd.read_csv(REPORT/'results/covid_roster187_local_ids.csv').sample_id)
ok=ro['sample'].isin(roster)&~ro.major.isin(['RBC','Platelet'])
samples=sorted(ro.loc[ok,'sample'].unique())
selection=np.sort(np.concatenate([rng.choice(np.where(ok&(ro['sample']==s))[0],size=min(400,int((ok&(ro['sample']==s)).sum())),replace=False) for s in samples]))
rm=ro.iloc[selection].reset_index(drop=True)
gen={};gvar=None
for K,branch in [('k1000','bulkonly_lora5_20261006'),('k4000','bulkonly_lora5_k4000_20261006')]:
    path=ROOT/branch/f'analysis/harmony/covid_{K}/covid_{K}_harmony_embedding.h5ad'
    with h5py.File(path,'r') as f:
        m=ad.io.read_elem(f['obs']);v=ad.io.read_elem(f['var']);assert 'counts' in f['layers']
        genes=v.gene_symbol.astype(str).to_numpy()
    if gvar is None:gvar=genes
    else:assert np.array_equal(gvar,genes)
    gi=np.sort(np.concatenate([rng.choice(np.asarray(ix),100,replace=False) for ix in m.groupby(['sample_id','draw'],observed=True).indices.values()]))
    gm=m.iloc[gi].reset_index(drop=True)
    gm['sample']=gm.sample_id.astype(str);gm['patient']=gm.donor_id.astype(str)
    gm['lin6']=gm.celltype_l1_pred.astype(str).map({'T cell':'T','B cell':'B','Monocyte':'Myeloid','NK cell':'NK','DC':'Myeloid','Other':'Other'})
    gen[K]=(path,gi,gm)
p1=pd.Series(np.arange(len(real_genes)),index=real_genes);p1=p1[~p1.index.duplicated()]
pg=pd.Series(np.arange(len(gvar)),index=gvar);pg=pg[~pg.index.duplicated()]
panel=[g for g in pg.index if g in p1.index];log('reading panel',len(panel),'real cells',len(selection))
cache=OUT/'inputs/mapping_cache';cache.mkdir(exist_ok=True)
if (cache/'real_counts.npz').exists():
    assert np.array_equal(np.load(cache/'real_rows.npy'),selection)
    X=sp.load_npz(cache/'real_counts.npz')
else:
    X=read_csr_rows(REAL,selection,p1.loc[panel].to_numpy());sp.save_npz(cache/'real_counts.npz',X);np.save(cache/'real_rows.npy',selection)
G={}
for K,(path,gi,gm) in gen.items():
    if (cache/f'{K}_counts.npz').exists():
        assert np.array_equal(np.load(cache/f'{K}_rows.npy'),gi);G[K]=sp.load_npz(cache/f'{K}_counts.npz')
    else:
        G[K]=read_csr_rows(path,gi,pg.loc[panel].to_numpy(),node='layers/counts');sp.save_npz(cache/f'{K}_counts.npz',G[K]);np.save(cache/f'{K}_rows.npy',gi)
log('counts loaded')
patients=np.array(sorted(rm.patient.unique()));rng_split=np.random.default_rng(7);rng_split.shuffle(patients);half=len(patients)//2
old=json.load(open(REPORT/'proportion_qualification/results/covid_mapping_meta.json'))
protocols=[('sample_overlap',0,np.isin(rm['sample'],old['ref_samples'])),
           ('patient_disjoint',0,np.isin(rm.patient,patients[:half])),('patient_disjoint',1,np.isin(rm.patient,patients[half:]))]
summary=[];confusion=[];cells=[];splits=[];generated=[];fitlog=[]
levels=['lin6','major','fine']
depth=np.concatenate([np.asarray(G[k].sum(1)).ravel() for k in G])
depth_by_lineage={}
for k in G:
    lineage=gen[k][2].lin6.to_numpy().copy()
    lineage[gen[k][2].aifi_l2.astype(str).eq('Plasma cell').to_numpy()]='Plasma'
    total=np.asarray(G[k].sum(1)).ravel()
    for c in np.unique(lineage):depth_by_lineage.setdefault(c,[]).append(total[lineage==c])
depth_by_lineage={c:np.concatenate(v) for c,v in depth_by_lineage.items()}

def evaluate(protocol,fold,source,route,level,truth,pred,md):
    df=pd.DataFrame({'patient':md.patient.to_numpy(),'sample':md['sample'].to_numpy(),'true':truth,'pred':pred})
    df['correct']=df.true==df.pred
    pacc=df.groupby('patient').correct.mean()
    est,lo,hi=boot_mean(pacc,rng)
    # Per-sample composition error -> patient averages, matching bootstrap unit.
    cats=sorted(set(truth)|set(pred));tp=sample_props(df,['patient','sample'],'true',cats);pp=sample_props(df,['patient','sample'],'pred',cats)
    l1=np.abs(tp-pp).sum(axis=1).groupby(level='patient').mean();lmean,llo,lhi=boot_mean(l1,rng)
    row=dict(protocol=protocol,fold=fold,source=source,route=route,level=level,n_cells=len(df),n_patients=df.patient.nunique(),
             patient_accuracy=est,accuracy_lo=lo,accuracy_hi=hi,pooled_accuracy=float(df.correct.mean()),
             balanced_accuracy=balanced_accuracy_score(truth,pred),composition_L1=lmean,L1_lo=llo,L1_hi=lhi)
    summary.append(row)
    counts=df.groupby(['patient','sample','true','pred']).size().rename('n').reset_index()
    counts['protocol']=protocol;counts['fold']=fold;counts['source']=source;counts['route']=route;counts['level']=level
    confusion.append(counts)

for protocol,fold,train in protocols:
    tag=f'{protocol}_fold{fold}'
    complete=OUT/'results'/f'{tag}_completed.json'
    if complete.exists():
        prior=json.load(open(complete));summary+=prior['summary'];splits+=prior['splits'];fitlog+=prior['fitlog'];log('loaded completed',tag);continue
    idx0=len(summary);sp0=len(splits);fit0=len(fitlog)
    Xa,Xb=X[train],X[~train];ma,mb=rm[train].reset_index(drop=True),rm[~train].reset_index(drop=True)
    overlap=set(ma.patient)&set(mb.patient)
    if protocol=='patient_disjoint':assert not overlap
    splits.append(dict(protocol=protocol,fold=fold,ref_patients=sorted(ma.patient.unique()),test_patients=sorted(mb.patient.unique()),
                       overlap_patients=sorted(overlap),ref_samples=sorted(ma['sample'].unique()),test_samples=sorted(mb['sample'].unique())))
    for role,md in [('reference',ma),('test',mb)]:
        md[['patient','sample','group']].drop_duplicates().assign(role=role,protocol=protocol,fold=fold).to_csv(OUT/'results'/f'{tag}_{role}_roster.csv',index=False)
    a=ad.AnnData(Xa);sc.pp.highly_variable_genes(a,flavor='seurat_v3',n_top_genes=2000);hv=np.where(a.var.highly_variable)[0]
    Ya=norm(Xa);mu=np.asarray(Ya[:,hv].mean(0)).ravel();sd=np.sqrt(np.maximum(np.asarray(Ya[:,hv].multiply(Ya[:,hv]).mean(0)).ravel()-mu**2,0))+1e-6
    def scale(Y):return np.clip((Y[:,hv].toarray()-mu)/sd,-10,10).astype(np.float32)
    pca=PCA(30,random_state=0,svd_solver='randomized').fit(scale(Ya))
    def proj(x):return pca.transform(scale(norm(x))).astype(np.float32)
    Pa=proj(Xa)
    def t2(x,lineage):
        tot=np.asarray(x.sum(1)).ravel();target=np.array([rng.choice(depth_by_lineage.get(c,depth)) for c in lineage]);P=proj(x)
        _,nb=NearestNeighbors(n_neighbors=15,n_jobs=2).fit(P).kneighbors(P)
        w=sp.csr_matrix((np.full(nb.size,1/15,np.float32),nb.ravel(),np.arange(0,nb.size+1,15)),shape=(len(nb),len(nb)))
        prof=(sp.diags(1/np.maximum(tot,1))@x).tocsr();sup=x.copy();sup.data[:]=1
        e=sp.vstack([(w[i:i+2000]@prof).multiply(sup[i:i+2000]).tocsr() for i in range(0,len(tot),2000)])
        e=(sp.diags(target/np.maximum(np.asarray(e.sum(1)).ravel(),1e-12))@e).tocsr();e.data=rng.poisson(e.data).astype(np.float32);e.eliminate_zeros();return e
    T2a,T2b=t2(Xa,ma.lin6.to_numpy()),t2(Xb,mb.lin6.to_numpy());Sa=scale(norm(T2a));models={}
    for lv in levels:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always',ConvergenceWarning)
            models[lv]=LogisticRegression(C=.05,max_iter=1000,solver='lbfgs').fit(Sa,ma[lv].to_numpy())
        fitlog.append(dict(protocol=protocol,fold=fold,level=lv,iterations=int(models[lv].n_iter_.max()),warnings=[str(w.message) for w in ws]))
        log(tag,'classifier fitted',lv,'iterations',models[lv].n_iter_)
    del Sa
    cats={lv:pd.Categorical(ma[lv]).categories for lv in levels}
    def vote(ref,q):
        d,ix=NearestNeighbors(n_neighbors=30,n_jobs=2).fit(ref).kneighbors(q);w=1/(d+1e-3);out={}
        for lv in levels:
            c=pd.Categorical(ma[lv],categories=cats[lv]).codes[ix];scores=np.zeros((len(q),len(cats[lv])))
            for j in range(30):np.add.at(scores,(np.arange(len(q)),c[:,j]),w[:,j])
            out[lv]=np.asarray(cats[lv])[scores.argmax(1)]
        return out
    for name,x,md in [('realU',Xb,mb),('realT2',T2b,mb)]+[(k,G[k],gen[k][2]) for k in G if protocol=='patient_disjoint']:
        px=proj(x);pred={'R2':vote(Pa,px),'R3':{lv:models[lv].predict(scale(norm(x))) for lv in levels}}
        if name in ('realU','realT2') and fold==0 and name=='realT2':
            pre=np.vstack([Pa,px]);labels=np.r_[np.zeros(len(Pa)),np.ones(len(px))].astype(str)
            h=harmonypy.run_harmony(pre,pd.DataFrame({'source':labels}),['source'],random_state=0,max_iter_harmony=20,verbose=False)
            z=np.asarray(h.Z_corr).T.astype(np.float32);pred['R1']=vote(z[:len(Pa)],z[len(Pa):]);log(tag,'Harmony control done')
        for route,p in pred.items():
            for lv in levels:
                if name.startswith('real'):evaluate(protocol,fold,name,route,lv,md[lv].to_numpy(),p[lv],md)
                else:
                    # Cross-fitted transfer: no real-reference visit of this patient in the mapper.
                    keep=~md.patient.isin(set(ma.patient))
                    tab=md.loc[keep,['sample','patient','draw','group','wave','WHO_severity','sex_flag']].copy()
                    tab['label']=p[lv][keep];tab['K']=name;tab['route']=route;tab['level']=lv;tab['fold']=fold
                    tab['trained_adapter']=[bool(md.loc[i,'donor_trained_'+str(md.loc[i,'draw'])]) for i in md.index[keep]]
                    generated.append(tab)
        log(tag,'mapped',name)
    partc=pd.concat(confusion,ignore_index=True);partc=partc[(partc.protocol==protocol)&(partc.fold==fold)]
    partc.to_csv(OUT/'results'/f'{tag}_confusion_counts.csv.gz',index=False)
    partg=pd.concat(generated,ignore_index=True) if generated else pd.DataFrame()
    if not partg.empty:partg[partg.fold==fold].to_csv(OUT/'results'/f'{tag}_generated_labels.csv.gz',index=False)
    json.dump({'summary':summary[idx0:],'splits':splits[sp0:],'fitlog':fitlog[fit0:]},open(complete,'w'),indent=2)
    log('completed',tag)
pd.DataFrame(summary).to_csv(OUT/'results/covid_mapping_validation.csv',index=False)
json.dump({'panel_genes':len(panel),'real_cells_per_sample_cap':400,'generated_cells_per_sample_draw':100,
           'HVG_fit':'reference counts only','PCA_fit':'reference only, randomized solver seed 0',
           'readout_depth':'empirical generated panel depth by six lineages (native scANVI plus plasma); no clinical labels',
           'R1_scope':'PBISC-like real controls, fold 0; R2/R3 two-fold all-patient validation',
           'split_rosters':splits,'classifier_fit_log':fitlog},open(OUT/'results/covid_mapping_methods.json','w'),indent=2)
log('all mapping outputs saved')
