# AI Agent Instructions for `PBISC_scRNAseq`

## 1. Active Priority Task: Figure 7 (First-Draw Anchored) Revision
If you are revising **Figure 7 (`Immune cell associations with COVID-19 severity and death in patients with ESKD`, 217 acute draws in 66 patients, first-draw mortality)** to resolve **Issues 1, 2, and 3**, start directly at the self-contained revision pack:

- **Primary Specification**: [`analyses/fig7_firstdraw_revision_pack_20261124/Figure7_FirstDraw_Revision_Guide_for_AI_Agent_20261124.md`](analyses/fig7_firstdraw_revision_pack_20261124/Figure7_FirstDraw_Revision_Guide_for_AI_Agent_20261124.md)
- **Bundled Ready-to-Plot CSVs (0 external dependencies)**: [`analyses/fig7_firstdraw_revision_pack_20261124/data/`](analyses/fig7_firstdraw_revision_pack_20261124/data/)
- **Portable Verification Script**:
  ```bash
  pip install -r requirements.txt
  python analyses/fig7_firstdraw_revision_pack_20261124/verify_and_recompute_pack.py
  ```

All paths inside `analyses/fig7_firstdraw_revision_pack_20261124/` are 100% repository-relative and OS-independent (Linux / macOS / Windows). Do **not** guess or approximate any statistical value; load directly from `analyses/fig7_firstdraw_revision_pack_20261124/data/*.csv`.

---

## 2. Scientific Figure Generation & House-Style Rules
For every publication or report figure generated or revised in this repository:
1. **Physical Dimensions (*Nature Portfolio* two-column standard)**:
   - Width = `180.0 mm` (`7.0866 in`), height `≤ 240.0 mm` (`9.4488 in`).
2. **Typography & Spines**:
   - Sans-serif font (`Arial` / `Helvetica` / `DejaVu Sans`), `pdf.fonttype = 42`, `svg.fonttype = 'none'`.
   - Panel labels: bold lowercase (`a`, `b`, `c`, ...; `8–9 pt`), panel titles `6.5–7.5 pt`, axis labels/ticks `5.5–6.5 pt` (strictly `≥ 5.0 pt` everywhere).
   - Remove top and right spines on data plots (`ax.spines[['top', 'right']].set_visible(False)`).
3. **Colorblind-Safe Source Palette**:
   - `Measured` (`Real` scRNA-seq): neutral slate grey (`#6A737D` / `#565B63`, diamond `D`).
   - `PBISC inferred` (`PBISC new`): sky/cerulean blue (`#168FE5` / `#0072B2`, circle `o`).
   - `Combined` / pooled estimate: deep navy (`#1B365D` / `#0B3C5D`, square `s`).
4. **Export & Visual Quality Control**:
   - Export vector `.pdf`, editable `.svg`, and high-resolution `.png` (`300–600 DPI`) together.
   - Visually inspect the rendered `.png` before finishing to verify zero text/marker collisions, zero clipped error bars, and exact correspondence with `verify_and_recompute_pack.py`.
