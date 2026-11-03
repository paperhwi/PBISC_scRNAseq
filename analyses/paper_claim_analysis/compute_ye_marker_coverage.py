from pathlib import Path
import h5py, numpy as np, pandas as pd
out=Path('paper_claim_analysis')
claims=pd.read_csv(out/'fixed_paper_claim_definitions.csv',dtype=str).fillna('')
markers=pd.read_csv(out/'ye_supplementary_table6_all_markers.csv')
base=Path('bulkonly_lora5_20261006/analysis/original_scRNA/sepsis_GSE279452_full/sepsis_GSE279452_recap.h5ad')
with h5py.File(base,'r') as f:
    x=f['var/_index'][:]
    genes={v.decode() if isinstance(v,bytes) else str(v) for v in x}
rows=[]
for cluster,g in markers.groupby('cluster',sort=False):
    panel=set(g.gene.astype(str))
    present=panel & genes
    rows.append({'author_cluster_name':cluster,'n_author_markers':len(panel),'n_present_in_original_expression':len(present),'n_missing_from_original_expression':len(panel-present),'fraction_present':len(present)/len(panel) if panel else np.nan,'missing_genes':';'.join(sorted(panel-present))})
pd.DataFrame(rows).to_csv(out/'ye_supplementary_table6_expression_coverage.csv',index=False)
print(pd.DataFrame(rows).head().to_string(index=False))
print('total marker coverage',sum(r['n_present_in_original_expression'] for r in rows),sum(r['n_author_markers'] for r in rows))
