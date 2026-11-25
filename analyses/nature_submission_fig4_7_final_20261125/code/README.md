# Reproducible Build & Verification Scripts (`Figures 4–7`)

All scripts in `code/` use relative paths anchored to this package directory (`Path(__file__).resolve().parent.parent`) and run cleanly on Windows, macOS, and Linux (`Python >= 3.10` with `numpy`, `pandas`, `scipy`, `statsmodels`, `openpyxl`, `python-docx`, `reportlab`, `pymupdf`).

1. **`01_export_figures_and_tables.py`**
   - When run inside the primary `PBISC` workspace, exports vector/high-DPI `Figures 4–7` into `../figures/` and extracts all `22` CSV tables into `../tables/`.
   - When run in a standalone GitHub clone, validates the bundled `../figures/` and `../tables/` files and regenerates `../tables/00_table_manifest.csv` and `../PBISC_SourceData_and_Summary_Tables_Fig4_7_20261125.xlsx`.
2. **`02_build_manuscript_docx_pdf.py`**
   - Compiles `../PBISC_Manuscript_Fig4_7_Final_20261125.md` with inline embedded `../figures/Figure{4,5,6,7}_*.png` into `../PBISC_Manuscript_Fig4_7_Final_20261125.docx` and `../PBISC_Manuscript_Fig4_7_Final_20261125.pdf`.
3. **`03_verify_all_fig4_7.py`**
   - Independently recomputes and verifies every statistical, tabular, graphic, OOXML, and bibliographic assertion across Figures 4–7 and writes `../VERIFICATION_FINAL_FIG4_7_20261125.json`.
