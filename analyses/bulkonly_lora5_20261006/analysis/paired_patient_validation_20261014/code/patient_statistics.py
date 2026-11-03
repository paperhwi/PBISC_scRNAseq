"""Patient-weighted COVID composition effects and cluster bootstrap.

One mean composition per patient and clinical state. A patient's visits in
both states and both sources receive the same bootstrap multiplicity.
No recovery samples enter infection or severity contrasts.
"""
from common import *
import json
from scipy.stats import spearmanr

rng = np.random.default_rng(20261014); NB=2000
cross=pd.read_csv(ROOT/'paper_claim_analysis/covid_sample_patient_id_reconciliation.csv')
cross=cross.set_index('sample_id')
meta=pd.read_csv(OUT/'inputs/k1000/analysis/features/imperial_eskd/meta.csv').set_index('sample_id')
real=pd.read_csv(REPORT/'results/covid_real_labels.csv.gz',keep_default_na=False).rename(columns={'label':'fine','grp':'group'})
real['patient']=real['sample'].map(cross.patient_id)
real['wave']=real['sample'].map(cross.published_year).map({2020:'w1',2021:'w2'})
# Severity is the patient's peak among ACTIVE positive visits, never recovery.
sevmap={'mild':0,'moderate':1,'severe':2,'critical':3}
rs=cross[cross.case_control.str.lower().eq('positive')].copy()
rs['sev']=rs.WHO_severity.astype(str).str.lower().map(sevmap)
peak=rs.groupby('patient_id').sev.max()
real['severity']=real.patient.map(peak).map(lambda x:'severe' if x>=2 else 'mild' if pd.notna(x) else 'unknown')
real['lin6']=real.major.map(lin6)
sources={'real':real}
for K in ['k1000','k4000']:
    t=pd.read_csv(REPORT/f'results/covid_{K}_transfer.csv.gz',keep_default_na=False).rename(columns={'sample_id':'sample','label':'fine'})
    t['patient']=t['sample'].map(meta.donor_id);t['wave']=t['sample'].map(meta.wave)
    t['severity']=t['sample'].map(meta.severity_peak).map({'severe/critical':'severe','mild/moderate':'mild'})
    t['lin6']=t.major.map(lin6);t['sex_flag']=t['sample'].map(meta.sex_flag).astype(bool)
    t['trained']=[bool(meta.loc[s,'donor_trained_d'+str(dr).replace('draw','')]) for s,dr in zip(t['sample'],t.draw)]
    sources[K]=t
    sources[K+'_untrained_clean']=t[~t.trained & ~t.sex_flag].copy()
assert all(t.patient.notna().all() for t in sources.values())

def assemble(t,lv,cats,contrast,wave=None):
    t=t[t.group.isin(['covid_positive','covid_negative'])].copy()
    if wave:t=t[t.wave.eq(wave)]
    if contrast=='severity':
        t=t[t.group.eq('covid_positive') & t.severity.isin(['mild','severe'])]
        t['state']=t.severity
    else:t['state']=t.group.map({'covid_positive':'case','covid_negative':'control'})
    p=sample_props(t,['patient','sample','state'],lv,cats)
    # Samples, then patients: never let patients with many visits dominate.
    p=p.groupby(level=['patient','state']).mean()
    return p

rows=[];effects=[];roster=[]
for lv in ['lin6','major','fine']:
    cats=sorted(real[lv].unique())
    for contrast in ['infection','severity']:
        for wave in [None,'w1','w2']:
            case,control=('severe','mild') if contrast=='severity' else ('case','control')
            props={k:assemble(t,lv,cats,contrast,wave) for k,t in sources.items()}
            if any(s not in props['real'].index.get_level_values('state') for s in [case,control]):continue
            patients=sorted(set().union(*(set(p.index.get_level_values('patient')) for p in props.values())))
            weights=rng.multinomial(len(patients),np.full(len(patients),1/len(patients)),size=NB)
            point={};boot={};support={}
            for name,p in props.items():
                means=[];bs=[];ns=[]
                for state in [case,control]:
                    if state not in p.index.get_level_values('state'):break
                    z=p.xs(state,level='state');ns.append(len(z));lc=clr(z)
                    arr=pd.DataFrame(lc,index=z.index).reindex(patients,fill_value=0).to_numpy()
                    present=np.isin(patients,z.index).astype(float);den=weights@present
                    v=(weights@arr)/np.maximum(den[:,None],1)
                    v[den==0]=np.nan;bs.append(v);means.append(lc.mean(0))
                if len(means)!=2:continue
                point[name]=means[0]-means[1];boot[name]=bs[0]-bs[1];support[name]=ns
                for state in [case,control]:
                    for pt in p.xs(state,level='state').index:roster.append(dict(level=lv,contrast=contrast,wave=wave or 'all',source=name,state=state,patient=pt))
            realp=props['real'].groupby(level='state').mean()
            use=realp.max(axis=0).to_numpy()>=.005
            if use.sum()<3:continue
            xr=point['real'][use];br=boot['real'][:,use]
            for name,e in point.items():
                if name=='real':continue
                valid=np.isfinite(br).all(1)&np.isfinite(boot[name]).all(1)
                ys=boot[name][valid][:,use];xs=br[valid]
                slopes=((xs-xs.mean(1,keepdims=True))*(ys-ys.mean(1,keepdims=True))).sum(1)/np.maximum(((xs-xs.mean(1,keepdims=True))**2).sum(1),1e-12)
                slope=float(np.polyfit(xr,e[use],1)[0]);lo,hi=np.quantile(slopes,[.025,.975])
                rows.append(dict(level=lv,contrast=contrast,wave=wave or 'all',source=name,slope=slope,slope_lo=lo,slope_hi=hi,
                                 spearman=spearmanr(xr,e[use]).statistic,sign_agree=float(np.mean(np.sign(xr)==np.sign(e[use]))),
                                 n_types=int(use.sum()),real_case_patients=support['real'][0],real_control_patients=support['real'][1],
                                 gen_case_patients=support[name][0],gen_control_patients=support[name][1],valid_bootstraps=int(valid.sum())))
                for i,c in enumerate(cats):
                    rlo,rhi=np.nanquantile(boot['real'][:,i],[.025,.975]);elo,ehi=np.nanquantile(boot[name][:,i],[.025,.975])
                    effects.append(dict(level=lv,contrast=contrast,wave=wave or 'all',source=name,label=c,real_effect=point['real'][i],real_lo=rlo,real_hi=rhi,effect=e[i],lo=elo,hi=ehi,abundant=bool(use[i])))
            print('patient statistics',lv,contrast,wave or 'all',flush=True)
pd.DataFrame(rows).to_csv(OUT/'results/covid_patient_effect_summary.csv',index=False)
pd.DataFrame(effects).to_csv(OUT/'results/covid_patient_effects.csv',index=False)
pd.DataFrame(roster).drop_duplicates().to_csv(OUT/'results/covid_statistics_patient_roster.csv',index=False)
old=pd.read_csv(REPORT/'proportion_qualification/results/F_summary.csv')
new=pd.DataFrame(rows)
z=old[(old.cohort=='covid') & old.source.isin(['k1000','k4000'])]
comparison=z.merge(new[(new.contrast=='infection')&(new.wave=='all')],on=['level','source'],suffixes=('_old_sample','_patient'))
comparison.to_csv(OUT/'results/covid_sample_vs_patient_effects.csv',index=False)
json.dump({'bootstrap_reps':NB,'unit':'patient','shared_patient_multiplicities_across_states_and_sources':True,
           'infection_adjustment':'wave-stratified sensitivity; no post-infection severity adjustment for total infection effect',
           'severity_definition':'peak severity in active positive samples only; no recovery',
           'labels':'existing transferred labels; patient-held-out classifier assessed separately'},open(OUT/'results/patient_statistics_methods.json','w'),indent=2)
print('Patient statistics saved',flush=True)
