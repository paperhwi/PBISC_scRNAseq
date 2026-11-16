"""Caveat check: Real vs PBISC within the 28 patients profiled by both (severity trend and per-patient agreement)."""
import os, pandas as pd
from scipy import stats
RUN=os.environ.get('RUN','/home/claude/run/'); os.makedirs('data',exist_ok=True)
S=pd.read_csv(RUN+'out/sample_level_composition_metadata.csv')
FE={'rare_B_ASC_dividing_pct':'B ASC dividing','rare_NK CD56 bright_pct':'NK CD56bright','rare_HSPC_pct':'HSPC'}
def pt(src):
    d=S[(S.source==src)&(S.condition=='positive')&S.in_primary_scope.astype(bool)&S.qc_pass_min200cells.astype(bool)]
    g=d.groupby('patient_id'); T=g[list(FE)].mean(); T['sev']=g.severity_ord.first(); return T
R=pt('Real'); B=pt('PBISC'); sh=R.index.intersection(B.index); rows=[]
for f,n in FE.items():
    a=stats.spearmanr(R.loc[sh,'sev'],R.loc[sh,f]); b=stats.spearmanr(B.loc[sh,'sev'],B.loc[sh,f]); c=stats.spearmanr(R.loc[sh,f],B.loc[sh,f])
    rows.append(dict(population=n,n_shared=len(sh),real_rho_sev=a[0],real_p=a[1],pbisc_rho_sev=b[0],pbisc_p=b[1],per_patient_rho=c[0],per_patient_p=c[1]))
T=pd.DataFrame(rows); T.to_csv('data/shared_patient_check.csv',index=False); print(T.round(3).to_string())
