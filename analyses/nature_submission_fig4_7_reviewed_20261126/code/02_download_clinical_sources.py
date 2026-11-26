#!/usr/bin/env python3
"""
02_download_clinical_sources.py
-------------------------------
Optional helper script (zero hardcoded paths) to download the original public
clinical workbooks for `van der Wijst et al. (2021)` (`PMC8601717.2` Data files
S1–S9) and `Stephenson et al. (2021)` (`Nature Medicine` Supplementary Tables)
and verify their SHA-256 checksums against
`verification_20261126/clinical_sources/source_manifest.json`.

Note: Pre-extracted clinical rows (`vdw_clinical_rows.json` and
`steph2021_patient_rows.json`) are already bundled in the repository so
`01_reproduce_adjudicated_tables.py` and `03_verify_all_fig4_7.py` run offline
without downloading large external archives.
"""

from __future__ import annotations
import hashlib
import io
import json
import urllib.request
import zipfile
from pathlib import Path

import openpyxl

PKG_DIR = Path(__file__).resolve().parent.parent
CLIN_DIR = PKG_DIR / "verification_20261126" / "clinical_sources"


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def main() -> None:
    manifest = json.loads((CLIN_DIR / "source_manifest.json").read_text(encoding="utf-8"))
    expected = {entry["file"]: entry for entry in manifest}

    for entry in manifest:
        if "url" not in entry:
            continue
        name = entry["file"]
        url = entry["url"]
        out_path = CLIN_DIR / name
        if not out_path.exists():
            print(f"Downloading {name} from {url} ...")
            with urllib.request.urlopen(url, timeout=120) as resp:
                out_path.write_bytes(resp.read())
        actual_sha = sha256_bytes(out_path.read_bytes())
        assert actual_sha == entry["sha256"], f"SHA-256 mismatch for {name}: {actual_sha} != {entry['sha256']}"
        print(f"[OK] {name}: SHA-256 verified ({actual_sha[:12]}...)")

        if out_path.suffix == ".zip":
            with zipfile.ZipFile(out_path, "r") as zf:
                member = next(x for x in zf.namelist() if "Data_file_S3.xlsx" in x)
                s3_bytes = zf.read(member)
                s3_target = CLIN_DIR / "van_der_Wijst_Data_file_S3.xlsx"
                s3_target.write_bytes(s3_bytes)
                assert sha256_bytes(s3_bytes) == expected["van_der_Wijst_Data_file_S3.xlsx"]["sha256"]
                wb = openpyxl.load_workbook(io.BytesIO(s3_bytes), read_only=True, data_only=True)
                rows = list(wb["ClinicalCharacteristics"].values)
                wb.close()
                print(f"[OK] Extracted van_der_Wijst_Data_file_S3.xlsx ({len(rows) - 1} rows)")


if __name__ == "__main__":
    main()
