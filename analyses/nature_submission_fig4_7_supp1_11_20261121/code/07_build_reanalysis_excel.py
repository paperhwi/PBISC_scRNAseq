"""Build a clean, multi-sheet Microsoft Excel (.xlsx) workbook with zero external dependencies
(pure Python OpenXML) bundling the complete 68-patient / 217-draw inventory, 39 COVID-negative
donors, and all side-by-side Patient-Mean vs Same-Day WHO Severity & First-Draw Death reanalysis tables.
"""
import pathlib
import sys
import zipfile
from xml.sax.saxutils import escape

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")
import numpy as np
import pandas as pd

ROOT = pathlib.Path(r"F:\PBISC\deliverables\nature_submission_fig4_6_20261120\draw_sameday_sev_firstdraw_death_20261120")
TAB_DIR = ROOT / "tables"
OUT_XLSX = ROOT / "PBISC_Patient_Draw_Inventory_and_SameDaySev_FirstDrawDeath_20261120.xlsx"


def col_letter(idx: int) -> str:
    """1-based column index -> Excel column letters (A, B, ..., Z, AA, ...)."""
    s = ""
    while idx > 0:
        idx, r = divmod(idx - 1, 26)
        s = chr(65 + r) + s
    return s


def is_num(val) -> bool:
    if isinstance(val, (bool, np.bool_)):
        return False
    if isinstance(val, (int, float, np.integer, np.floating)):
        return np.isfinite(val)
    return False


def df_to_sheet_xml(df: pd.DataFrame) -> str:
    ncols = len(df.columns)
    nrows = len(df) + 1
    ref = f"A1:{col_letter(max(ncols, 1))}{max(nrows, 1)}"

    # Compute sensible column widths
    col_widths = []
    for c in df.columns:
        max_len = len(str(c))
        for v in df[c].head(250):
            if pd.isna(v):
                continue
            s = f"{v:.4g}" if isinstance(v, (float, np.floating)) and np.isfinite(v) else str(v)
            max_len = max(max_len, min(len(s), 44))
        col_widths.append(min(max(max_len + 3, 9), 46))

    cols_xml = "<cols>" + "".join(
        f'<col min="{i+1}" max="{i+1}" width="{w:.1f}" customWidth="1"/>'
        for i, w in enumerate(col_widths)
    ) + "</cols>"

    rows_xml = []
    # Header row (style 1 = bold navy fill with white text)
    h_cells = []
    for ci, cname in enumerate(df.columns, start=1):
        h_cells.append(f'<c r="{col_letter(ci)}1" t="inlineStr" s="1"><is><t>{escape(str(cname))}</t></is></c>')
    rows_xml.append(f'<row r="1" ht="22" customHeight="1">{"".join(h_cells)}</row>')

    # Data rows (zebra striping styles 2 & 3)
    for ri, row in enumerate(df.itertuples(index=False, name=None), start=2):
        st = "3" if (ri % 2 == 1) else "2"
        cells = []
        for ci, val in enumerate(row, start=1):
            coord = f"{col_letter(ci)}{ri}"
            if val is None or (isinstance(val, (float, np.floating)) and not np.isfinite(val)) or (isinstance(val, str) and val == "nan"):
                cells.append(f'<c r="{coord}" s="{st}"/>')
            elif is_num(val):
                v_str = str(int(val)) if isinstance(val, (int, np.integer)) else f"{float(val):.10g}"
                cells.append(f'<c r="{coord}" s="{st}"><v>{v_str}</v></c>')
            else:
                txt = escape(str(val)).replace("\x00", "")
                cells.append(f'<c r="{coord}" t="inlineStr" s="{st}"><is><t>{txt}</t></is></c>')
        rows_xml.append(f'<row r="{ri}">{"".join(cells)}</row>')

    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="{ref}"/>'
        '<sheetViews><sheetView workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        '</sheetView></sheetViews>'
        '<sheetFormatPr defaultRowHeight="16"/>'
        f'{cols_xml}'
        f'<sheetData>{"".join(rows_xml)}</sheetData>'
        f'<autoFilter ref="{ref}"/>'
        '</worksheet>'
    )


def write_xlsx(sheets: list[tuple[str, pd.DataFrame]], out_path: pathlib.Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(out_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        # [Content_Types].xml
        sheet_overrides = "".join(
            f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            for i in range(1, len(sheets) + 1)
        )
        zf.writestr(
            "[Content_Types].xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            f'{sheet_overrides}</Types>',
        )
        # _rels/.rels
        zf.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
            '</Relationships>',
        )
        # xl/styles.xml (0=default, 1=navy bold header, 2=plain, 3=light blue-gray stripe)
        zf.writestr(
            "xl/styles.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<fonts count="2">'
            '<font><sz val="10"/><name val="Arial"/></font>'
            '<font><b/><sz val="10"/><color rgb="FFFFFFFF"/><name val="Arial"/></font>'
            '</fonts>'
            '<fills count="4">'
            '<fill><patternFill patternType="none"/></fill>'
            '<fill><patternFill patternType="gray125"/></fill>'
            '<fill><patternFill patternType="solid"><fgColor rgb="FF0F4D92"/></patternFill></fill>'
            '<fill><patternFill patternType="solid"><fgColor rgb="FFF4F7FB"/></patternFill></fill>'
            '</fills>'
            '<borders count="2">'
            '<border><left/><right/><top/><bottom/><diagonal/></border>'
            '<border><bottom style="thin"><color rgb="FFD9E1EC"/></bottom></border>'
            '</borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="4">'
            '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
            '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1" applyAlignment="1"><alignment vertical="center"/></xf>'
            '<xf numFmtId="0" fontId="0" fillId="0" borderId="1" xfId="0" applyBorder="1"/>'
            '<xf numFmtId="0" fontId="0" fillId="3" borderId="1" xfId="0" applyFill="1" applyBorder="1"/>'
            '</cellXfs></styleSheet>',
        )
        # xl/workbook.xml + xl/_rels/workbook.xml.rels
        s_tags = []
        r_tags = []
        for i, (sname, df) in enumerate(sheets, start=1):
            clean_name = sname[:31].replace("/", "_").replace(":", "_")
            s_tags.append(f'<sheet name="{escape(clean_name)}" sheetId="{i}" r:id="rId{i}"/>')
            r_tags.append(
                f'<Relationship Id="rId{i}" '
                'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                f'Target="worksheets/sheet{i}.xml"/>'
            )
            zf.writestr(f"xl/worksheets/sheet{i}.xml", df_to_sheet_xml(df))

        zf.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f'<sheets>{"".join(s_tags)}</sheets></workbook>',
        )
        zf.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{"".join(r_tags)}'
            f'<Relationship Id="rId{len(sheets)+1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
            '</Relationships>',
        )


def main():
    readme = pd.DataFrame([
        ("Workbook Purpose", "Complete ESKD COVID-19 Real + PBISC Patient/Blood-Draw Inventory AND Full Side-by-Side Reanalysis: Previous (Patient-Average) vs NEW (All-Draws Same-Day WHO Severity + Donor First-Draw Final Death)"),
        ("Primary Cohort", "N = 66 non-CLL COVID-19+ ESKD patients (217 primary blood draws): Real-only n=10 (47 draws, 1 death), Shared-27 [Real used] n=27 (91 draws, 4 deaths), PBISC new [bulk-only] n=29 (79 draws, 6 deaths). Total verified fatal cases = 11 (TRUE_PUB_11; 9 direct COVID-19)."),
        ("Excluded & Control Donors", "2 excluded CLL patients (C116, C139) + 39 COVID-negative single-draw ESKD donors = 107 total Imperial ESKD donors. (Sepsis Ye et al. GSE279452 is strictly 1 draw/donor cross-sectional.)"),
        ("NEW Severity Framework", "ALL 217 individual blood draws paired with contemporaneous same-day WHO severity (WHO_temp_severity: 0=mild, 1=moderate, 2=severe, 3=critical; 0 missing). Patient-clustered Exchangeable GEE (unadjusted & +week time-adjusted; Huber-White & Mancl-DeRouen bias-reduced SE), Mundlak within-patient vs between-patient decomposition, and patient-block clustered bootstrap Spearman rho."),
        ("NEW Death Framework", "ONE baseline observation per patient: Donor FIRST blood draw paired with final verified death (TRUE_PUB_11; N=66, 11 deaths; Real n=37/5d, PBISC new n=29/6d). Includes exact Mann-Whitney AUC/P, same-day-WHO-adjusted first-draw death, early first-draw <=7d (N=53, 10 deaths), Direct COVID-9, WHO 2-3, and direct contrast vs late (>7d) draws."),
        ("Sheet: T1_Patients68", "All 68 COVID-19+ patients (66 primary + 2 excluded CLL): cohort arm, sampling days (from +swab & from symptom onset), draw count, same-day WHO at each draw, peak WHO, mortality flags, steroid record."),
        ("Sheet: T1b_Draws217", "All 217 primary individual blood draws with draw_order, is_first_draw, sampling day, same-day WHO, peak WHO, death flags, hospital labs, and key immune subset percentages."),
        ("Sheet: T1c_CohortSummary", "Cohort-level summary statistics of draw counts, first-draw days, same-day WHO distributions, and mortality across Real-only, Shared-27, Real ALL, PBISC new, and Combined (Survived vs Died)."),
        ("Sheet: T1d_CovidNeg39", "All 39 COVID-negative single-draw ESKD donors (Real scRNA / PBISC bulk status)."),
        ("Sheet: T1e_All362_Draws", "Complete unfiltered inventory of all 362 samples in Real scRNA (198) + PBISC (164), including Shared-27 PBISC draws, CLL, COVID-negative, and recovery samples."),
        ("Sheet: T2b_Comp43_Combined", "Master side-by-side comparison across all 43 fine-grained immune subsets (Combined cohort): Previous Patient-Mean vs NEW Same-Day WHO GEE/Mundlak/Clust-Boot Spearman & NEW First-Draw Final Death."),
        ("Sheet: T2_Comp43_AllArms", "Full 43-subset comparison stratified by Real, PBISC new, and Combined."),
        ("Sheet: T3_NewFunnel43", "5-step credibility & PBISC-rescue classification under the NEW framework."),
        ("Sheet: T4_PowerRescue", "300-iteration cohort-expansion subsampling power rescue (n=20 Real vs +29 PBISC new vs +29 shuffled) under BOTH Same-Day WHO GEE and First-Draw Final Death."),
        ("Sheet: T5_HSPC_16Mod_20Gene", "Single-cell HSPC 16 transcriptional modules (Ery5, Cycle4, HSC5, Myeloid4, ...) and 20 canonical genes reanalyzed under Same-Day WHO & First-Draw Final Death."),
        ("Sheet: T6_Tcell_ISG13", "Single-cell T-cell 19-gene ISG across 13 T-cell subsets reanalyzed under Same-Day WHO & First-Draw Final Death."),
        ("Sheet: T7_Labs_SameDay", "Real-cohort hospital blood laboratory correlations (CRP, ALC, NLR, Ferritin, D-dimer, Troponin) on same-day blood-draw pairs via patient-clustered GEE."),
    ], columns=["Section / Sheet", "Description"])

    sheet_specs = [
        ("README", readme),
        ("T1_Patients68", pd.read_csv(TAB_DIR / "T1_all68_patients_draw_schedule.csv")),
        ("T1b_Draws217", pd.read_csv(TAB_DIR / "T1b_all217_primary_draws_detail.csv")),
        ("T1c_CohortSummary", pd.read_csv(TAB_DIR / "T1c_cohort_draw_summary.csv")),
        ("T1d_CovidNeg39", pd.read_csv(TAB_DIR / "T1d_covid_negative_39_donors.csv")),
        ("T1e_All362_Draws", pd.read_csv(TAB_DIR / "T1e_all362_samples_real_and_pbisc.csv")),
        ("T2b_Comp43_Combined", pd.read_csv(TAB_DIR / "T2b_comp43_combined_full_comparison.csv")),
        ("T2_Comp43_AllArms", pd.read_csv(TAB_DIR / "T2_comp43_patientmean_vs_samedaySev_firstdrawDeath.csv")),
        ("T3_NewFunnel43", pd.read_csv(TAB_DIR / "T3_new_framework_43_funnel.csv")),
        ("T4_PowerRescue", pd.read_csv(TAB_DIR / "T4_power_rescue_new_framework.csv")),
        ("T5_HSPC_16Mod_20Gene", pd.read_csv(TAB_DIR / "T5_HSPC_16modules_20genes_new_framework.csv")),
        ("T6_Tcell_ISG13", pd.read_csv(TAB_DIR / "T6_Tcell_ISG_13subsets_new_framework.csv")),
        ("T7_Labs_SameDay", pd.read_csv(TAB_DIR / "T7_hospital_labs_sameday_pairs.csv")),
    ]
    write_xlsx(sheet_specs, OUT_XLSX)
    print(f"[OK] Wrote {OUT_XLSX} ({OUT_XLSX.stat().st_size / 1024:.1f} KB, {len(sheet_specs)} sheets)")


if __name__ == "__main__":
    main()
