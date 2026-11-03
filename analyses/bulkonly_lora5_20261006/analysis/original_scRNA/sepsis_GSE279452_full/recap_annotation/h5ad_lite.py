"""Minimal AnnData (on-disk format 0.1/0.2) writer using only h5py + numpy + pandas.
Runs on the linked computer, which has no anndata/scipy. Output reads with anndata.read_h5ad / scanpy."""
import h5py, numpy as np, pandas as pd

STR = h5py.string_dtype(encoding="utf-8")


def _attrs(obj, etype, ever):
    obj.attrs["encoding-type"] = etype; obj.attrs["encoding-version"] = ever


def write_strings(g, name, values, **kw):
    d = g.create_dataset(name, data=np.asarray([str(v) for v in values], dtype=object), dtype=STR, **kw)
    _attrs(d, "string-array", "0.2.0"); return d


def write_array(g, name, values, **kw):
    d = g.create_dataset(name, data=np.asarray(values), **kw); _attrs(d, "array", "0.2.0"); return d


def write_categorical(g, name, values, categories=None):
    s = pd.Categorical(pd.Series(values).astype(str), categories=categories)
    grp = g.create_group(name); _attrs(grp, "categorical", "0.2.0"); grp.attrs["ordered"] = False
    codes = s.codes.astype(np.int8 if len(s.categories) < 127 else np.int16 if len(s.categories) < 32000 else np.int32)
    write_array(grp, "codes", codes); write_strings(grp, "categories", s.categories)


def write_column(g, name, values):
    v = pd.Series(values)
    if v.dtype == bool or str(v.dtype) == "boolean":
        write_array(g, name, v.fillna(False).to_numpy().astype(bool))
    elif pd.api.types.is_numeric_dtype(v):
        write_array(g, name, v.to_numpy())
    else:
        write_categorical(g, name, v.fillna("NA").astype(str))


def add_obs_columns(f, df):
    """append columns to an existing obs dataframe (index must already match)."""
    obs = f["obs"]; order = list(obs.attrs["column-order"])
    for c in df.columns:
        if c in obs:
            del obs[c]
        write_column(obs, c, df[c])
        if c not in order:
            order.append(c)
    obs.attrs["column-order"] = np.array(order, dtype=object)


def write_dataframe(f, key, df, index):
    g = f.create_group(key); _attrs(g, "dataframe", "0.2.0")
    g.attrs["_index"] = "_index"; write_strings(g, "_index", index)
    for c in df.columns:
        write_column(g, c, df[c])
    g.attrs["column-order"] = np.array(list(df.columns), dtype=object) if len(df.columns) else np.array([], dtype="float64")


def empty_dict(f, key):
    g = f.require_group(key); _attrs(g, "dict", "0.1.0"); return g


class CSRWriter:
    """stream rows into X (csr_matrix)."""
    def __init__(self, f, n_cols, key="X", chunk=8_000_000):
        self.g = f.create_group(key); _attrs(self.g, "csr_matrix", "0.1.0")
        self.data = self.g.create_dataset("data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(chunk,))
        self.ind = self.g.create_dataset("indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(chunk,))
        self.ip = [0]; self.n_cols = n_cols; self.nnz = 0
    def add(self, data, indices, indptr):
        n = len(data)
        self.data.resize((self.nnz + n,)); self.data[self.nnz:] = data
        self.ind.resize((self.nnz + n,)); self.ind[self.nnz:] = indices
        self.ip.extend((np.asarray(indptr[1:], dtype=np.int64) + self.nnz).tolist()); self.nnz += n
    def close(self):
        self.g.create_dataset("indptr", data=np.asarray(self.ip, dtype=np.int64))
        self.g.attrs["shape"] = np.array([len(self.ip) - 1, self.n_cols], dtype=np.int64)
