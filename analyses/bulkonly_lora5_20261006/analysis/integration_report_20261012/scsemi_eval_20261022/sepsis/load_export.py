"""Read an extracted strict export. Requires NumPy; no pickle is used."""
from pathlib import Path
import argparse
import ast
import zipfile
import numpy as np

def open_export(unit_directory):
    """Return NpzFile; use as a context manager. Accessing counts loads the dense array."""
    return np.load(Path(unit_directory) / "export" / "scsemi.npz", allow_pickle=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("unit_directory", type=Path)
    args = parser.parse_args()
    with open_export(args.unit_directory) as arrays:
        print("samples:", len(arrays["sample_ids"]))
        print("genes:", len(arrays["gene_ids"]))
        print("cells per sample:", sorted(set(arrays["cells_per_bulk"].tolist())))
    # Inspect the expression shape without loading the dense expression matrix.
    with zipfile.ZipFile(args.unit_directory / "export" / "scsemi.npz") as archive:
        with archive.open("generated_counts.npy") as handle:
            assert handle.read(6) == b"\x93NUMPY"
            version = handle.read(2)
            size = int.from_bytes(handle.read(2 if version[0] == 1 else 4), "little")
            header = ast.literal_eval(handle.read(size).decode("latin1"))
            print("expected-count matrix:", header["shape"], header["descr"])
