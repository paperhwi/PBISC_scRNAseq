import gzip
from pathlib import Path

for path in sorted(Path(__file__).parent.rglob("*.gz")):
    try:
        with gzip.open(path, "rb") as stream:
            while stream.read(8 * 1024 * 1024):
                pass
        print("OK", path.name, path.stat().st_size, flush=True)
    except Exception as error:
        print("BAD", path.name, path.stat().st_size, repr(error), flush=True)
