import fsspec
import h5py

URL = "https://covid19.cog.sanger.ac.uk/eskd_covid19.h5ad"
with fsspec.open(URL, "rb", block_size=64 * 1024, cache_type="bytes") as stream:
    print("remote size", stream.seek(0, 2), flush=True)
    stream.seek(0)
    with h5py.File(stream, "r") as h5:
        obs = h5["obs"]
        print("obs columns", list(obs.keys()), flush=True)
        for key in obs.keys():
            item = obs[key]
            if isinstance(item, h5py.Dataset):
                print(key, "dataset", item.shape, item.dtype, flush=True)
            else:
                print(key, "group", list(item.keys()), flush=True)
