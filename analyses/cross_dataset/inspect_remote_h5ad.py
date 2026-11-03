import fsspec,h5py
url='https://covid19.cog.sanger.ac.uk/eskd_covid19.h5ad'
f=fsspec.open(url,'rb',block_size=1024*1024,cache_type='readahead').open()
print('size',f.seek(0,2),flush=True); f.seek(0)
with h5py.File(f,'r') as h:
 print('keys',list(h.keys()),flush=True)
 for k in h.keys():
  x=h[k]
  print(k,type(x).__name__,list(x.keys())[:15] if isinstance(x,h5py.Group) else x.shape,flush=True)
f.close()
