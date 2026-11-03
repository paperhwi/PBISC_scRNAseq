"""Disease effects after patient-held-out R2/R3 label transfer."""
from common import *
from scipy.stats import spearmanr
rng=np.random.default_rng(147);NB=2000
g=pd.concat([pd.read_csv(f) for f in (OUT/'results').glob('patient_disjoint_fold*_generated_labels.csv.gz')],ignore_index=True)
r=pd.read_csv(REPORT/'results/covid_real_labels.csv.gz').rename(columns={'label':'fine','grp':'group'})
cross=pd.read_csv(ROOT/'paper_claim_analysis/covid_sample_patient_id_reconciliation.csv').set_index('sample_id')
r['patient']=r['sample'].map(cross.patient_id);r['lin6']=r.major.map(lin6)
rows=[];effects=[]
for route in ['R2','R3']:
    for lv in ['lin6','major','fine']:
        cats=sorted(r[lv].unique());rr=r[r.group.isin(['covid_positive','covid_negative'])]
        tr=sample_props(rr,['patient','sample','group'],lv,cats).groupby(level=['patient','group']).mean()
        gp={}
        for K in ['k1000','k4000']:
            for clean in [False,True]:
                t=g[(g.K==K)&(g.route==route)&(g.level==lv)&g.group.isin(['covid_positive','covid_negative'])]
                if clean:t=t[~t.trained_adapter.astype(bool)&~t.sex_flag.astype(bool)]
                # Equal weighting of both out-of-reference mapper folds where both are eligible.
                p=sample_props(t,['patient','sample','group','fold'],'label',cats).groupby(level=['patient','sample','group']).mean().groupby(level=['patient','group']).mean()
                gp[K+('_clean' if clean else '')]=p
        allp={'real':tr,**gp};patients=sorted(set().union(*(set(p.index.get_level_values('patient')) for p in allp.values())))
        w=rng.multinomial(len(patients),np.full(len(patients),1/len(patients)),size=NB);e={};be={}
        for name,p in allp.items():
            es=[];bs=[]
            for state in ['covid_positive','covid_negative']:
                z=p.xs(state,level='group');lc=clr(z);pres=np.isin(patients,z.index)
                ar=pd.DataFrame(lc,index=z.index).reindex(patients,fill_value=0).to_numpy();den=w@pres
                b=(w@ar)/np.maximum(den[:,None],1);b[den==0]=np.nan;bs.append(b);es.append(lc.mean(0))
            e[name]=es[0]-es[1];be[name]=bs[0]-bs[1]
        use=tr.groupby(level='group').mean().max(axis=0).to_numpy()>=.005
        for name,p in gp.items():
            x=e['real'][use];y=e[name][use];xb=be['real'][:,use];yb=be[name][:,use]
            good=np.isfinite(xb).all(1)&np.isfinite(yb).all(1);xb=xb[good];yb=yb[good]
            xc=xb-xb.mean(1,keepdims=True);yc=yb-yb.mean(1,keepdims=True);b=(xc*yc).sum(1)/np.maximum((xc**2).sum(1),1e-12)
            lo,hi=np.quantile(b,[.025,.975])
            rows.append(dict(route=route,level=lv,source=name,slope=np.polyfit(x,y,1)[0],slope_lo=lo,slope_hi=hi,
                             spearman=spearmanr(x,y).statistic,sign_agree=np.mean(np.sign(x)==np.sign(y)),n_types=int(use.sum()),
                             gen_patients=p.index.get_level_values('patient').nunique(),valid_bootstraps=int(good.sum())))
            for j,c in enumerate(cats):
                el,eh=np.nanquantile(be[name][:,j],[.025,.975]);effects.append(dict(route=route,level=lv,source=name,label=c,real_effect=e['real'][j],effect=e[name][j],lo=el,hi=eh,abundant=use[j]))
        print('Cross-fitted disease effects',route,lv,flush=True)
pd.DataFrame(rows).to_csv(OUT/'results/covid_crossfit_effect_summary.csv',index=False)
pd.DataFrame(effects).to_csv(OUT/'results/covid_crossfit_effects.csv',index=False)
