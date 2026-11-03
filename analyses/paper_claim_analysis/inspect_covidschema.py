import h5py, numpy as np
p=r'bulkonly_lora5_20261006/analysis/original_scRNA/covid_eskd/eskd_covid19_recap.h5ad'
f=h5py.File(p,'r')
print('root',list(f.keys()))
print('obs keys',list(f['obs'].keys()))
print('attrs',dict(f['obs'].attrs))
for k in f['obs'].keys():
 o=f['obs'][k]
 if isinstance(o,h5py.Group) and 'categories' in o:
  cats=[x.decode() if isinstance(x,bytes) else str(x) for x in o['categories'][:]]
  print(k,'cats',cats[:100],len(cats))
 elif isinstance(o,h5py.Dataset) and o.ndim==1 and len(o)<1000:
  print(k,'uniq',np.unique(o[:])[:40])
f.close()
