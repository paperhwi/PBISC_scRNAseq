import pandas as pd
from pathlib import Path
root=Path('paper_claim_analysis')
xlsx=root/'41590_2025_2345_MOESM9_ESM.xlsx'
df=pd.read_excel(xlsx,header=1)
df=df.rename(columns={'allmarkers.gene':'gene','avg_log2FC':'avg_log2fc','p_val':'p_value','p_val_adj':'p_adj'})
df['gene']=df['gene'].astype(str).str.strip()
df['cluster']=df['cluster'].astype(str).str.strip()
df.to_csv(root/'ye_supplementary_table6_all_markers.csv',index=False)
agg=(df.groupby('cluster',sort=False).agg(n_markers=('gene','size'),n_unique_genes=('gene','nunique'),n_up=('avg_log2fc',lambda x:int((x>0).sum())),n_down=('avg_log2fc',lambda x:int((x<0).sum())),min_adj_p=('p_adj','min'),max_pct1=('pct.1','max')).reset_index())
agg.to_csv(root/'ye_supplementary_table6_cluster_marker_counts.csv',index=False)
sets=(df.groupby('cluster',sort=False)['gene'].apply(lambda x:';'.join(x)).rename('all_markers').reset_index())
sets=agg.merge(sets,on='cluster',how='left')
sets.to_csv(root/'ye_supplementary_table6_marker_sets_by_cluster.csv',index=False)
# Exact author subset labels for the seven fixed claims, including original published cluster spelling.
claim_to_cluster={
'SEP-1':'T_CD4_c02-GPR183-NR4A2',
'SEP-2':'T_CD8_c05-GNLY-CCL4',
'SEP-3':'NKT_c05-IFNG-CCL4',
'SEP-3-NK':'NK_c03-PRF1-CCL4',
'SEP-4A':'Mono-c01-CD14-HLA-DRA',
'SEP-4B':'Mono-c03-CD14-RETN',
'SEP-5':'Mono-c06-CD14-MKI67',
}
old=pd.read_csv(root/'fixed_paper_claim_definitions.csv')
for c in ['core_markers','author_cluster_name','author_marker_count']:
 if c not in old.columns: old[c]=''
old['core_markers']=old['markers']
for i,r in old.iterrows():
 if r['study']=='Ye et al. 2026':
  cluster=claim_to_cluster.get(r['id'],'')
  old.at[i,'author_cluster_name']=cluster
  sub=df[df.cluster.eq(cluster)]
  old.at[i,'author_marker_count']=str(int(sub.gene.nunique()) if len(sub) else 0)
  old.at[i,'markers']=';'.join(sub.gene.astype(str)) if len(sub) else ''
  old.at[i,'definition_source']='full author marker roster from Supplementary Table 6; exact publisher cluster label mapped; previous core panel retained in core_markers'
old.to_csv(root/'fixed_paper_claim_definitions.csv',index=False)
print('rows',len(df),'clusters',df.cluster.nunique(),'unique genes',df.gene.nunique())
print('claim mappings:')
print(old[old.study.eq('Ye et al. 2026')][['id','author_cluster_name','author_marker_count','core_markers']].to_string(index=False))
print('summary marker sizes:',agg.n_markers.describe().to_dict())
