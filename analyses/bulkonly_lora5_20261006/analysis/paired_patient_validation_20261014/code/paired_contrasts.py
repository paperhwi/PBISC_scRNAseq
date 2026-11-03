"""Paired method differences, retaining covariance on common patients."""
from common import *
from scipy.stats import wilcoxon
df=pd.read_csv(OUT/'results/paired_sample_errors.csv.gz');rng=np.random.default_rng(14);rows=[]
for (coh,scope,lv),d in df[df.draw.eq(-1)].groupby(['cohort','scope','level']):
    p=d.pivot(index=['sample_id','donor_id'],columns='method',values='L1')
    for method in p:
        if method.startswith('PBISC_'):
            for base in p:
                if base==method or base.startswith('PBISC_'):continue
                z=p[[method,base]].dropna();delta=(z[method]-z[base]).groupby(level='donor_id').mean()
                mean,lo,hi=boot_mean(delta,rng)
                rows.append(dict(cohort=coh,scope=scope,level=lv,method=method,baseline=base,n_samples=len(z),n_donors=len(delta),
                                 method_L1=z[method].groupby(level='donor_id').mean().mean(),baseline_L1=z[base].groupby(level='donor_id').mean().mean(),
                                 paired_delta_L1=mean,lo=lo,hi=hi,wilcoxon_p=wilcoxon(delta).pvalue))
pd.DataFrame(rows).to_csv(OUT/'results/paired_method_contrasts.csv',index=False)
print('Paired contrasts saved')
