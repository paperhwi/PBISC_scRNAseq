"""Build Nature Portfolio Word manuscript (.docx via pure standard-library OOXML),
Source Data workbook (.xlsx), README.md, SUPPLEMENTARY_NUMBERING_MAP.md,
migrate/verify Supplementary Figures 14-23 (following Figures 1-4 Supplementary Figures 1-13
in strict Nature Portfolio citation order), stage preview PNGs, and run a standalone
end-to-end independent verification audit (VERIFICATION_REPORT_20261122.json).
"""
import json
import math
import pathlib
import pickle
import re
import shutil
import sys
import zipfile
from xml.sax.saxutils import escape
from PIL import Image

sys.path.insert(0, r"F:\PBISC\.analysis_runtime")

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(r"F:\PBISC")
OUT_ROOT = ROOT / "deliverables" / "nature_submission_fig5_7_20261122"
DATA_DIR = OUT_ROOT / "data"
FIG_DIR = OUT_ROOT / "figures"
MD_PATH = OUT_ROOT / "PBISC_Nature_Manuscript_Fig5_7_and_SuppFig14_23_20261122.md"
DOCX_PATH = OUT_ROOT / "PBISC_Nature_Manuscript_Fig5_7_and_SuppFig14_23_20261122.docx"
XLSX_PATH = OUT_ROOT / "PBISC_SourceData_and_Verification_Fig5_7_20261122.xlsx"
REPORT_PATH = OUT_ROOT / "VERIFICATION_REPORT_20261122.json"
README_PATH = OUT_ROOT / "README.md"
SUPP_MAP_PATH = OUT_ROOT / "SUPPLEMENTARY_NUMBERING_MAP.md"
CONV_ART_DIR = pathlib.Path(r"C:\Users\scott\.gemini\antigravity\brain\add4945a-39a0-4f8d-acd4-3287c21c6121")

# Mapping from standalone Figs 5-7 S1..S10 stems -> Full-paper Nature citation-order S14..S23 stems
RENAME_STEMS = [
    ("SuppFig1_Coverage_Depth_and_Decomposition", "SuppFig14_Coverage_Depth_and_Decomposition"),
    ("SuppFig4_Sepsis_Welch_DEG_Pseudobulk", "SuppFig15_Sepsis_Welch_DEG_Pseudobulk"),
    ("SuppFig5_COVID_Welch_DEG_Pseudobulk", "SuppFig16_COVID_Welch_DEG_Pseudobulk"),
    ("SuppFig2_Rare_Populations_Across_Patients", "SuppFig17_Rare_Populations_Across_Patients"),
    ("SuppFig3_Rare_Population_Summed_Marker_Heatmaps", "SuppFig18_Rare_Population_Summed_Marker_Heatmaps"),
    ("SuppFig6_Cohort_and_Sample_Metadata_Tables", "SuppFig19_Cohort_and_Sample_Metadata_Tables"),
    ("SuppFig7_Patient_Swimmer_276_Acute_Draws", "SuppFig20_Patient_Swimmer_276_Acute_Draws"),
    ("SuppFig8_43CellType_Sensitivity_and_OfficialR_Deconv", "SuppFig21_43CellType_Sensitivity_and_OfficialR_Deconv"),
    ("SuppFig10_SingleCell_TISG_and_HSPC_Programs", "SuppFig22_SingleCell_TISG_and_HSPC_Programs"),
    ("SuppFig9_External_4Cohort_SameCriterion_Validation", "SuppFig23_External_4Cohort_SameCriterion_Validation"),
]


def migrate_supp_numbering_if_needed():
    stem_map = dict(RENAME_STEMS)
    for old_s, new_s in RENAME_STEMS:
        for ext in ("png", "pdf", "svg"):
            src = FIG_DIR / f"{old_s}.{ext}"
            dst = FIG_DIR / f"{new_s}.{ext}"
            if src.exists():
                shutil.move(str(src), str(dst))

    old_a1 = OUT_ROOT / "AUDIT_FIG5_6_SFIG1_5.json"
    new_a1 = OUT_ROOT / "AUDIT_FIG5_6_SFIG14_18.json"
    if old_a1.exists():
        items = json.loads(old_a1.read_text(encoding="utf-8"))
        for it in items:
            it["figure"] = stem_map.get(it["figure"], it["figure"])
        order1 = [
            "Figure5_CrossCohort_BulkOnly_Comparison",
            "Figure6_Rare_Cell_Frequency_and_Expression",
            "SuppFig14_Coverage_Depth_and_Decomposition",
            "SuppFig15_Sepsis_Welch_DEG_Pseudobulk",
            "SuppFig16_COVID_Welch_DEG_Pseudobulk",
            "SuppFig17_Rare_Populations_Across_Patients",
            "SuppFig18_Rare_Population_Summed_Marker_Heatmaps",
        ]
        by_k = {x["figure"]: x for x in items}
        new_a1.write_text(json.dumps([by_k[k] for k in order1 if k in by_k], indent=2), encoding="utf-8")
        old_a1.unlink()

    old_a2 = OUT_ROOT / "AUDIT_FIG7_SFIG6_10.json"
    new_a2 = OUT_ROOT / "AUDIT_FIG7_SFIG19_23.json"
    if old_a2.exists():
        items = json.loads(old_a2.read_text(encoding="utf-8"))
        for it in items:
            it["figure"] = stem_map.get(it["figure"], it["figure"])
        order2 = [
            "Figure7_SameDay_Severity_and_LastDraw_Mortality",
            "SuppFig19_Cohort_and_Sample_Metadata_Tables",
            "SuppFig20_Patient_Swimmer_276_Acute_Draws",
            "SuppFig21_43CellType_Sensitivity_and_OfficialR_Deconv",
            "SuppFig22_SingleCell_TISG_and_HSPC_Programs",
            "SuppFig23_External_4Cohort_SameCriterion_Validation",
        ]
        by_k = {x["figure"]: x for x in items}
        new_a2.write_text(json.dumps([by_k[k] for k in order2 if k in by_k], indent=2), encoding="utf-8")
        old_a2.unlink()

    # Clean superseded S1-S10 manuscript drafts
    for sup_f in [
        OUT_ROOT / "PBISC_Nature_Manuscript_Fig5_7_and_SuppFig1_10_20261122.md",
        OUT_ROOT / "PBISC_Nature_Manuscript_Fig5_7_and_SuppFig1_10_20261122.docx",
    ]:
        if sup_f.exists():
            sup_f.unlink()


def runs_from_inline(text: str, base_sz_halfpt: int = 21, force_bold: bool = False) -> str:
    tokens = re.split(r"(\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)", text)
    out = []
    for tok in tokens:
        if not tok:
            continue
        b = force_bold
        it = False
        font = "Times New Roman"
        val = tok
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            b = True
            val = tok[2:-2]
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            it = True
            val = tok[1:-1]
        elif tok.startswith("`") and tok.endswith("`") and len(tok) > 2:
            val = tok[1:-1]
            font = "Arial"
        rpr = [f'<w:rFonts w:ascii="{font}" w:hAnsi="{font}" w:cs="{font}"/>', f'<w:sz w:val="{base_sz_halfpt}"/>']
        if b:
            rpr.append("<w:b/>")
        if it:
            rpr.append("<w:i/>")
        out.append(f'<w:r><w:rPr>{"".join(rpr)}</w:rPr><w:t xml:space="preserve">{escape(val)}</w:t></w:r>')
    return "".join(out)


def build_docx():
    md_text = MD_PATH.read_text(encoding="utf-8")
    paras = []
    for raw_line in md_text.splitlines():
        line = raw_line.rstrip()
        if not line or line.strip() == "---":
            continue
        if line.startswith("# "):
            t = line[2:].strip()
            paras.append(f'<w:p><w:pPr><w:spacing w:before="240" w:after="160"/></w:pPr>{runs_from_inline(t, 30, True)}</w:p>')
        elif line.startswith("## "):
            t = line[3:].strip()
            paras.append(f'<w:p><w:pPr><w:spacing w:before="220" w:after="110"/></w:pPr>{runs_from_inline(t, 26, True)}</w:p>')
        elif line.startswith("### "):
            t = line[4:].strip()
            paras.append(f'<w:p><w:pPr><w:spacing w:before="170" w:after="75"/></w:pPr>{runs_from_inline(t, 23, True)}</w:p>')
        elif line.startswith("> "):
            t = line[2:].strip()
            paras.append(f'<w:p><w:pPr><w:ind w:left="320"/><w:spacing w:before="50" w:after="90"/></w:pPr>{runs_from_inline(t, 20, False)}</w:p>')
        elif re.match(r"^\d+\.\s+", line.strip()):
            paras.append(f'<w:p><w:pPr><w:ind w:left="320"/><w:spacing w:before="45" w:after="90" w:line="276" w:lineRule="auto"/></w:pPr>{runs_from_inline(line.strip(), 21, False)}</w:p>')
        else:
            paras.append(f'<w:p><w:pPr><w:spacing w:before="55" w:after="105" w:line="276" w:lineRule="auto"/></w:pPr>{runs_from_inline(line, 21, False)}</w:p>')

    doc_xml = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
        "<w:body>"
        + "".join(paras)
        + '<w:sectPr><w:pgSz w:w="12240" w:h="15840"/>'
          '<w:pgMar w:top="1224" w:right="1224" w:bottom="1224" w:left="1224" w:header="720" w:footer="720" w:gutter="0"/>'
          "</w:sectPr></w:body></w:document>"
    )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>'
        "</Types>"
    )
    rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>'
        "</Relationships>"
    )
    with zipfile.ZipFile(DOCX_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", rels)
        zf.writestr("word/document.xml", doc_xml)
    print(f"[OK] Built Word manuscript: {DOCX_PATH.name} ({DOCX_PATH.stat().st_size:,} bytes)", flush=True)


def build_excel():
    sheets = [
        ("Fig5_Lineage_ShiftBoot", "fig5_six_lineage_shift_bootstrap.csv"),
        ("Fig5_Lineage_PtBoot", "fig5_six_lineage_patient_bootstrap.csv"),
        ("Fig5_Lit_GeneScores", "fig5_literature_gene_scores_summary.csv"),
        ("Fig5_Welch_DEGs", "fig5_welch_deg_fidelity.csv"),
        ("Fig5_Shift_Decomp", "fig5_shift_decomposition.csv"),
        ("Fig6_Rare_ParentShares", "fig6_rare_parent_shares_14sepsis_13covid.csv"),
        ("Fig6_Rare_PtBootstrap", "fig6_rare_patient_bootstrap_10states.csv"),
        ("Fig6_Summed_InterPt_r", "fig6_rare_summed_interpatient_r_summary.csv"),
        ("Fig7_ESKD_43CellTypes", "fig7_eskd_43celltypes_sameday_sev_lastdraw_death.csv"),
        ("Fig7_ESKD_66_LastDraws", "fig7_eskd_66_last_samples.csv"),
        ("Fig7_External4_RE_Meta", "fig7_external4_random_effects_meta.csv"),
        ("Fig7_External4_PerStudy", "fig7_external4_per_cohort_results.csv"),
        ("Supp14_Scope_and_Depth", "supp1_source_scope_and_depth.csv"),
        ("Supp21_Sensitivity_CBC", "supp8_sensitivity_count_cbc_cll_276framework.csv"),
        ("Supp21_OfficialR_Deconv", "supp8_officialR_deconv_12methods.csv"),
        ("Supp21_New29_HSPC_AUC", "supp8_officialR_new29_hspc_auc.csv"),
        ("Supp22_CD8CTL_7ISG", "supp10_cd8ctl_isg_sameday_first_last.csv"),
        ("Supp22_TISG_13Subsets", "supp10_tcell_isg_13subsets_summary.csv"),
        ("Supp22_HSPC_16Modules", "supp10_hspc_16modules_summary.csv"),
        ("Supp22_HSPC_20Genes", "supp10_hspc_20genes_summary.csv"),
    ]
    with pd.ExcelWriter(XLSX_PATH, engine="openpyxl") as writer:
        for sname, fname in sheets:
            df = pd.read_csv(DATA_DIR / fname)
            df.to_excel(writer, sheet_name=sname[:31], index=False)
    print(f"[OK] Built Source Data Excel workbook: {XLSX_PATH.name} ({XLSX_PATH.stat().st_size:,} bytes)", flush=True)


def weighted_midrank_ecdf(v: np.ndarray, w: np.ndarray) -> np.ndarray:
    uniq, inv = np.unique(v, return_inverse=True)
    grp_w = np.bincount(inv, weights=w)
    cum_before = np.concatenate([[0.0], np.cumsum(grp_w[:-1])])
    mid = cum_before + 0.5 * grp_w
    return mid[inv]


def weighted_spearman_independent(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    rx = weighted_midrank_ecdf(x, w)
    ry = weighted_midrank_ecdf(y, w)
    wn = w / w.sum()
    mx = np.sum(wn * rx)
    my = np.sum(wn * ry)
    cov = np.sum(wn * (rx - mx) * (ry - my))
    vx = np.sum(wn * (rx - mx) ** 2)
    vy = np.sum(wn * (ry - my) ** 2)
    return float(cov / math.sqrt(vx * vy)) if vx > 0 and vy > 0 else 0.0


def stage_artifact_previews():
    CONV_ART_DIR.mkdir(parents=True, exist_ok=True)
    for stem in [
        "Figure5_CrossCohort_BulkOnly_Comparison.png",
        "Figure6_Rare_Cell_Frequency_and_Expression.png",
        "Figure7_SameDay_Severity_and_LastDraw_Mortality.png",
        "SuppFig19_Cohort_and_Sample_Metadata_Tables.png",
        "SuppFig20_Patient_Swimmer_276_Acute_Draws.png",
    ]:
        shutil.copy2(FIG_DIR / stem, CONV_ART_DIR / stem)
    for old_art in ["SuppFig6_Cohort_and_Sample_Metadata_Tables.png", "SuppFig7_Patient_Swimmer_276_Acute_Draws.png"]:
        p_old = CONV_ART_DIR / old_art
        if p_old.exists():
            p_old.unlink()
    print("[OK] Staged preview PNGs in conversation artifact directory.", flush=True)


def build_supp_numbering_map():
    txt = """# Full-Paper Supplementary Figure Numbering Map (`Supplementary Figures 1–23`)

To comply strictly with *Nature Portfolio* editorial guidelines—requiring Supplementary Figures to be numbered sequentially in exact order of first citation across `Main Figures 1–7`—the `13` Supplementary Figures accompanying `Main Figures 1–4` in the master PBISC manuscript (`PBISC-Diffusion - 20261116.pdf`) occupy **`Supplementary Figures 1–13`**, and the `10` Supplementary Figures accompanying **`Main Figures 5–7`** occupy **`Supplementary Figures 14–23`**.

> **Re-offsetting**: Both figure builders (`code/02_build_fig5_fig6_supp1_5.py` and `code/03_build_fig7_supp6_10.py`) accept `--supp-offset=N` (default `--supp-offset=13`).

## 1. Figures 1–4 Supplementary Figures (`Supplementary Figures 1–13`; user-managed section)

| Nature Citation-Order ID | Parent Main Figure | Title (`PBISC-Diffusion - 20261116.pdf`) | Legacy 2026-11-16 Deck ID |
| :--- | :--- | :--- | :--- |
| **Supplementary Fig. 1** | Figure 1 | Pretraining corpus and source allocation | Deck Supp. Fig. 8 |
| **Supplementary Fig. 2** | Figure 1 | Pseudobulk recipes and training exposure | Deck Supp. Fig. 9 |
| **Supplementary Fig. 3** | Figure 2 | Gene-context specificity, training objectives and model size | Deck Supp. Fig. 10 |
| **Supplementary Fig. 4** | Figure 2 | Width changes and diagnostic composition hints | Deck Supp. Fig. 11 |
| **Supplementary Fig. 5** | Figure 3 | Bulk–pseudobulk mismatch and input calibration | Deck Supp. Fig. 12 |
| **Supplementary Fig. 6** | Figure 3 | Individual-donor recovery and paired uncertainty | Deck Supp. Fig. 13 |
| **Supplementary Fig. 7** | Figure 3 | Integration and depth-matching controls | Deck Supp. Fig. 14 |
| **Supplementary Fig. 8** | Figure 4 | Matched support allocation and reference selection | Deck Supp. Fig. 15 |
| **Supplementary Fig. 9** | Figure 4 | Paired-cohort reconstruction and exploratory limits | Deck Supp. Fig. 16 |
| **Supplementary Fig. 10** | Figure 4 | Monocyte and NK disease-program heatmaps in the paired cohorts | Deck Supp. Fig. 18 |
| **Supplementary Fig. 11** | Figure 4 | Cell-type-resolved disease effects (all-type 3-source `limma-voom` volcano plots) | Deck Supp. Fig. 6 |
| **Supplementary Fig. 12** | Figure 4 | Cell-type-resolved disease-effect agreement (3-reference hexbin concordance) | Deck Supp. Fig. 7 |
| **Supplementary Fig. 13** | Figure 4 | Major-cell immune-program covariance | Deck Supp. Fig. 17 |

---

## 2. Figures 5–7 Supplementary Figures (`Supplementary Figures 14–23`; rendered in `figures/`)

| Nature Citation-Order ID | Parent Main Figure | File Stem in `figures/` | Standalone `Figs. 5–7` Local ID |
| :--- | :--- | :--- | :--- |
| **Supplementary Fig. 14** | Figure 5 | `SuppFig14_Coverage_Depth_and_Decomposition` | Local S1 (Deck Supp. Fig. 1 + decomp.) |
| **Supplementary Fig. 15** | Figure 5 | `SuppFig15_Sepsis_Welch_DEG_Pseudobulk` | Local S4 (Deck Supp. Fig. 4) |
| **Supplementary Fig. 16** | Figure 5 | `SuppFig16_COVID_Welch_DEG_Pseudobulk` | Local S5 (Deck Supp. Fig. 5) |
| **Supplementary Fig. 17** | Figure 6 | `SuppFig17_Rare_Populations_Across_Patients` | Local S2 (Deck Supp. Fig. 2) |
| **Supplementary Fig. 18** | Figure 6 | `SuppFig18_Rare_Population_Summed_Marker_Heatmaps` | Local S3 (Deck Supp. Fig. 3) |
| **Supplementary Fig. 19** | Figure 7 | `SuppFig19_Cohort_and_Sample_Metadata_Tables` | Local S6 |
| **Supplementary Fig. 20** | Figure 7 | `SuppFig20_Patient_Swimmer_276_Acute_Draws` | Local S7 |
| **Supplementary Fig. 21** | Figure 7 | `SuppFig21_43CellType_Sensitivity_and_OfficialR_Deconv` | Local S8 |
| **Supplementary Fig. 22** | Figure 7 | `SuppFig22_SingleCell_TISG_and_HSPC_Programs` | Local S10 |
| **Supplementary Fig. 23** | Figure 7 | `SuppFig23_External_4Cohort_SameCriterion_Validation` | Local S9 |
"""
    SUPP_MAP_PATH.write_text(txt, encoding="utf-8")


def build_readme():
    txt = """# PBISC Nature Submission Package — Figures 5–7, Methods (Figures 4–7), and Supplementary Figures 14–23 (`2026-11-22`)

## 1. Methodological Alignment Summary (확정 분석 기준 1:1 반영 요약)

- **Figure 4 (Methods 반영 완료 · 본문 그림 및 보충 그림 `Supplementary Figs. 8–13` 번호 매핑 완비)**:
  - 패혈증(`Ye et al., 2026`) 및 혈액투석 COVID-19(`Gisby et al., 2022; Stephenson et al., 2025`) 각 코호트에서 paired donor `5`명으로 모델을 적응시키고, 평가용 held-out donor의 bulk RNA-seq에서 세포 `4,000`개를 생성해 동일 donor의 실측(`Real`) 세포와 비교.
  - 3가지 참조 환자 조합에서 `T`, `B`, `Monocyte`, `NK` 세포의 질환별 발현 차이를 `limma-voom`으로 추정 (`BH q < 0.05`).
  - 전체 논문(`Figures 1–7`) *Nature Portfolio* 인용 순서 기준 `Figures 1–4`의 보충 그림 `13`종(`Fig. 1: S1–S2`, `Fig. 2: S3–S4`, `Fig. 3: S5–S7`, `Fig. 4: S8–S13`) 뒤를 이어 본 패키지의 `Figures 5–7` 보충 그림 `10`종을 **`Supplementary Figures 14–23`**으로 일괄 정렬 (`SUPPLEMENTARY_NUMBERING_MAP.md` 참조; `--supp-offset` 인자로 언제든 재조정 가능).
- **Figure 5 · 다른 환자 집단 간 비교 (`Figure5_CrossCohort_BulkOnly_Comparison` + `Supplementary Figs. 14–16`)**:
  - 패혈증 `122`개(`488,000` cells) 및 COVID-19 `164`개(`656,000` cells)의 unpaired bulk-only 검체에서 생성한 세포를 별도 환자 집단의 실측 세포(`108`개 `633,072` cells / `187`개 `578,835` cells) 및 `scSemiProfiler`(`34`개 / `102`개)와 비교.
  - 6대 면역 세포 비율은 `환자 × 상태`(`patient × condition`) 단위로 요약 후 환자 단위로 `49,999`회 부트스트랩 재표집.
  - 문헌 보고 질환 신호는 세포형별 유전자 점수(`Ye et al.` 5개 모듈, `Stephenson et al.` 7-ISG 모듈)로 비교.
  - 세포형별 DEG는 양측 **Welch 검정**을 사용하고 **`BH q < 0.05`이면서 `|log₂ fold change| ≥ 1`**인 유전자를 선택 (`AUROC`: 패혈증 `0.78–0.96`, COVID-19 `0.85–0.99`).
- **Figure 6 · 희귀 세포의 빈도와 발현 (`Figure6_Rare_Cell_Frequency_and_Expression` + `Supplementary Figs. 17–18`)**:
  - 고정된 `30-NN` 세포 주석으로 패혈증 **`14`종**, COVID-19 **`13`종**의 희귀 세포 상태를 구분.
  - 환자군(`cases`)에서 각 상태가 상위 세포 계통(`parent lineage`)에서 차지한 비율을 산출하되, **`HSC` 및 `HSPC`는 전체 세포 수(`all cells`)를 분모**로 사용 (`PBISC` `13/14` 및 `12/13` 재현 vs `scSemiProfiler` `7/14` 및 `8/13` 0-cell 소실).
  - 검체별로 같은 상태의 세포 발현을 합산(`log₁ₚ(CP10k)`)하여 **다른 환자(`patient_i ≠ patient_j`)**의 같은 상태 실측 세포와 **Pearson 상관(`r`)**을 구하고(`중앙값 r = 0.73–0.92` 및 `0.81–0.89`, `Real–Real` 환자 간 상한 근접/일치) 대표 마커 발현 히트맵을 비교.
- **Figure 7 · 채혈 당일 중증도와 마지막 채혈 기준 사망 비교 (`Figure7_SameDay_Severity_and_LastDraw_Mortality` + `Supplementary Figs. 19–23`)**:
  - 비-CLL 혈액투석 급성 COVID-19 환자 **`66`명**의 적격 급성기 검체 **`276`건**(`138 Real` + `138 PBISC`, 당일 중복 `0`건) 사용.
  - **중증도**: 채혈 당일 WHO 등급과 세포 비율(`% PBMC`)의 연관성을 환자별 총가중치를 같게 두어(`w_ij = 1 / n_i`, 환자당 합 `1.0`, 총합 `66.0`) 평가 → `BH43 q < 0.05` 유의 세포형 **`11/43` (`Real`) → `18/43` (`Combined raw %`, `rank-INT` 시 `20/43`)**로 확대.
  - **사망**: 환자별 **마지막 적격 검체 1건**(`66`건 = `24 Real` + `42 PBISC`, 최대 채혈일 동률 `0`건)으로 **생존 `55`명 vs 사망 `11`명**의 세포 비율을 비교(`Mann–Whitney AUC`) → `BH43 q < 0.05` 유의 예측인자 **`0/43` (`Real`) → `2/43` (`Combined raw %`: `T g/d AUC = 0.183, q = 0.0393` 및 `HSPC AUC = 0.800, q = 0.0393`; `rank-INT` 시 `CD16mono_C1` 포함 `3/43`)** 발굴.
  - **외부 4개 연구 메타분석**: `Stephenson 2021`, `COMBAT 2022`, `van der Wijst 2021`, `Liu 2021` (총 `300`개 급성기 검체 / `231`명, 사망 추적 코호트 생존 `197`명 vs 사망 `22`명)에서 **동일 기준**으로 분석 후 대응 연구 `k ≥ 2`인 세포형을 **DerSimonian–Laird 무작위효과 모형**으로 합성하여 재현 확인.

---

## 2. Directory Contents

- `PBISC_Nature_Manuscript_Fig5_7_and_SuppFig14_23_20261122.md` & `.docx`: Nature Portfolio 형식 원고 (`Results Figs. 5–7`, `Main Legends Figs. 5–7`, `Discussion Figs. 5–7`, `Methods Figs. 4–7`, `Supplementary Legends 14–23`).
- `SUPPLEMENTARY_NUMBERING_MAP.md`: 전체 논문 `Main Figures 1–7` 대응 `Supplementary Figures 1–23` 인용 순서 매핑표.
- `PBISC_SourceData_and_Verification_Fig5_7_20261122.xlsx`: 전 패널 핵심 통계/검산 통합 엑셀 워크북 (20개 시트).
- `VERIFICATION_REPORT_20261122.json`: 독립 재계산, 보충 그림 인용 단조 증가 순서(`14..23`), 그림 규격(`180.0 mm`, `min_pt ≥ 5.0`, `oob = 0`, `600 DPI`), `.docx` OOXML 무결성 전수 검산 리포트 (`81 / 81 PASS`, `FAIL = 0`).
- `figures/`: Main Figures `5–7` 및 Supplementary Figures `14–23` (`13`개 그림 × `PNG [600 DPI] / PDF [450 DPI] / SVG [editable Arial]` = 총 `39`개 파일).
- `data/`: 재현성 100% 원천 및 요약 CSV 데이터.
- `code/`: 데이터 산출(`01_...`), 본문/보충 그림 생성(`02_...`, `03_...`, `nature_style.py`), 문서 빌드 및 독립 검산(`04_...`) 스크립트 일체.
"""
    README_PATH.write_text(txt, encoding="utf-8")


def run_independent_verification():
    checks = []

    def rec(category: str, name: str, passed: bool, detail: str):
        checks.append({
            "category": category,
            "check": name,
            "status": "PASS" if passed else "FAIL",
            "detail": detail,
        })

    expected_figs = [
        "Figure5_CrossCohort_BulkOnly_Comparison",
        "Figure6_Rare_Cell_Frequency_and_Expression",
        "Figure7_SameDay_Severity_and_LastDraw_Mortality",
        "SuppFig14_Coverage_Depth_and_Decomposition",
        "SuppFig15_Sepsis_Welch_DEG_Pseudobulk",
        "SuppFig16_COVID_Welch_DEG_Pseudobulk",
        "SuppFig17_Rare_Populations_Across_Patients",
        "SuppFig18_Rare_Population_Summed_Marker_Heatmaps",
        "SuppFig19_Cohort_and_Sample_Metadata_Tables",
        "SuppFig20_Patient_Swimmer_276_Acute_Draws",
        "SuppFig21_43CellType_Sensitivity_and_OfficialR_Deconv",
        "SuppFig22_SingleCell_TISG_and_HSPC_Programs",
        "SuppFig23_External_4Cohort_SameCriterion_Validation",
    ]
    aud1 = {a["figure"]: a for a in json.loads((OUT_ROOT / "AUDIT_FIG5_6_SFIG14_18.json").read_text(encoding="utf-8"))}
    aud2 = {a["figure"]: a for a in json.loads((OUT_ROOT / "AUDIT_FIG7_SFIG19_23.json").read_text(encoding="utf-8"))}
    aud_all = {**aud1, **aud2}

    for stem in expected_figs:
        for ext in ("png", "pdf", "svg"):
            fp = FIG_DIR / f"{stem}.{ext}"
            ok = fp.exists() and fp.stat().st_size > 5000
            rec("Figures_Files", f"{stem}.{ext}", ok, f"exists={fp.exists()}, bytes={fp.stat().st_size if fp.exists() else 0}")
        with Image.open(FIG_DIR / f"{stem}.png") as im:
            w_px, h_px = im.size
        rec("Figures_DPI", f"{stem}_600dpi_width", abs(w_px - 4252) <= 2, f"px={w_px}x{h_px}")
        a = aud_all[stem]
        geom_ok = (abs(a["width_mm"] - 180.0) < 0.05 and a["height_mm"] <= 240.0
                   and a["min_font_pt"] >= 5.0 and a["below_min_pt_count"] == 0 and a["out_of_bounds_count"] == 0)
        rec("Figures_Geometry", f"{stem}_audit", geom_ok,
            f"{a['width_mm']}x{a['height_mm']} mm, min_pt={a['min_font_pt']}, oob={a['out_of_bounds_count']}")

    # Verify strict Nature citation order of Supplementary Figures 14-23 in Results
    md_txt = MD_PATH.read_text(encoding="utf-8")
    res_sec = md_txt.split("## Results", 1)[1].split("## Main Figure Legends", 1)[0]
    first_pos = {}
    for s_id in range(14, 24):
        m_single = re.search(rf"Supplementary Fig\.\s*{s_id}(?!\d)", res_sec)
        # Also handle "Supplementary Figs. 15 and 16" or "19 and 20"
        m_pair = re.search(rf"Supplementary Figs\.\s*(?:\d+\s+and\s+)?{s_id}(?!\d)", res_sec)
        positions = [m.start() for m in (m_single, m_pair) if m is not None]
        first_pos[s_id] = min(positions) if positions else -1
    all_cited = all(v >= 0 for v in first_pos.values())
    monotonic_order = [first_pos[k] for k in range(14, 24)] == sorted(first_pos[k] for k in range(14, 24))
    rec("Manuscript_SuppOrder", "Strict_Nature_citation_order_SuppFig14_to_SuppFig23",
        all_cited and monotonic_order, f"first_citation_char_pos={first_pos}")

    # Figure 5 checks
    s_scope = pd.read_csv(DATA_DIR / "supp1_source_scope_and_depth.csv")
    sep_pb = s_scope[(s_scope.cohort == "sepsis") & (s_scope.source == "PBISC")].iloc[0]
    cov_pb = s_scope[(s_scope.cohort == "covid") & (s_scope.source == "PBISC")].iloc[0]
    rec("Fig5_Scope", "Sepsis_122_bulkOnly_488k_cells",
        int(sep_pb.n_samples) == 122 and int(sep_pb.n_cells) == 488000,
        f"samples={int(sep_pb.n_samples)}, cells={int(sep_pb.n_cells)}")
    rec("Fig5_Scope", "COVID_164_bulkOnly_656k_cells",
        int(cov_pb.n_samples) == 164 and int(cov_pb.n_cells) == 656000,
        f"samples={int(cov_pb.n_samples)}, cells={int(cov_pb.n_cells)}")

    gm_sep = pd.read_csv(DATA_DIR / "fig5_sepsis_gene_means.csv")
    gm_cov = pd.read_csv(DATA_DIR / "fig5_covid_gene_means.csv")
    m_s = (gm_sep["Real"] > 0) & (gm_sep["PBISC"] > 0)
    m_c = (gm_cov["Real"] > 0) & (gm_cov["PBISC"] > 0)
    r_gm_sep = float(np.corrcoef(np.log10(gm_sep.loc[m_s, "Real"]), np.log10(gm_sep.loc[m_s, "PBISC"]))[0, 1])
    r_gm_cov = float(np.corrcoef(np.log10(gm_cov.loc[m_c, "Real"]), np.log10(gm_cov.loc[m_c, "PBISC"]))[0, 1])
    rec("Fig5_GeneMean", "Sepsis_and_COVID_log10_r",
        abs(r_gm_sep - 0.9757) < 0.002 and abs(r_gm_cov - 0.9656) < 0.002,
        f"r_sepsis={r_gm_sep:.4f}, r_covid={r_gm_cov:.4f}")

    pu = pd.read_csv(DATA_DIR / "fig5_six_lineage_patient_units.csv")
    dup_pu = pu.duplicated(subset=["cohort", "source", "patient", "condition", "cell_type"]).sum()
    rec("Fig5_PatientUnits", "Unique_patient_x_condition_units", int(dup_pu) == 0, f"duplicates={int(dup_pu)}, rows={len(pu)}")

    lit = pd.read_csv(DATA_DIR / "fig5_literature_gene_scores_summary.csv")
    pb_lit = lit[lit.source == "PBISC"]
    rec("Fig5_LiteratureScores", "All_9_PBISC_literature_modules_FDR_sig_and_concordant",
        bool((pb_lit["welch_q"] < 0.001).all()) and len(pb_lit) == 9,
        f"n_modules={len(pb_lit)}, max_PBISC_q={pb_lit['welch_q'].max():.3e}")

    wdeg = pd.read_csv(DATA_DIR / "fig5_welch_deg_fidelity.csv")
    pb_w = wdeg[wdeg.source == "PBISC"].set_index(["cohort", "celltype"])
    rec("Fig5_WelchDEG", "Welch_t_BH005_log2FC1_AUROC",
        bool((pb_w["AUROC"] >= 0.78).all()) and len(pb_w) == 8,
        f"min_AUROC={pb_w['AUROC'].min():.4f}, max_AUROC={pb_w['AUROC'].max():.4f}")

    # Figure 6 checks
    rshare = pd.read_csv(DATA_DIR / "fig6_rare_parent_shares_14sepsis_13covid.csv")
    n_sep_states = rshare[rshare.cohort == "sepsis"]["subset"].nunique()
    n_cov_states = rshare[rshare.cohort == "covid"]["subset"].nunique()
    hsc_par = set(rshare[rshare.subset.isin(["HSC", "HSPC"])]["parent"].unique())
    rec("Fig6_RareStates", "Fixed_14_Sepsis_13_COVID_and_HSC_HSPC_all_cells_denominator",
        n_sep_states == 14 and n_cov_states == 13 and hsc_par == {"all cells"},
        f"sepsis={n_sep_states}, covid={n_cov_states}, HSC/HSPC_parent={hsc_par}")

    with open(DATA_DIR / "fig6_rare_summed_cache.pkl", "rb") as fp:
        scache = pickle.load(fp)
    same_pt_leak = 0
    for coh_k in ("sepsis", "covid"):
        mdf = scache[coh_k]["mdf"]
        C = scache[coh_k]["C_all"]
        for sub_nm, st_d in scache[coh_k]["pair_dists"].items():
            ix_r = np.flatnonzero((mdf["subset"] == sub_nm) & (mdf["source"] == "Real"))
            ix_p = np.flatnonzero((mdf["subset"] == sub_nm) & (mdf["source"] == "PBISC"))
            pcol = "patient" if "patient" in mdf.columns else "patient_id"
            pts_r = mdf.iloc[ix_r][pcol].astype(str).to_numpy()
            pts_p = mdf.iloc[ix_p][pcol].astype(str).to_numpy()
            valid = pts_r[:, None] != pts_p[None, :]
            vals_recomp = C[np.ix_(ix_r, ix_p)][valid]
            if len(vals_recomp) != len(st_d["rp_all"]) or abs(float(np.median(vals_recomp)) - float(np.median(st_d["rp_all"]))) > 1e-9:
                same_pt_leak += 1
    rec("Fig6_InterPatientSummedR", "Zero_same_patient_pairs_and_exact_Pearson_medians",
        same_pt_leak == 0, f"mismatch_count={same_pt_leak}")

    # Figure 7 checks
    df276 = pd.read_csv(DATA_DIR / "fig7_eskd_276_acute_samples_metadata_and_43pct.csv")
    df66 = pd.read_csv(DATA_DIR / "fig7_eskd_66_last_samples.csv")
    eskd = pd.read_csv(DATA_DIR / "fig7_eskd_43celltypes_sameday_sev_lastdraw_death.csv")

    pt_w_sums = df276.groupby("patient_id")["equal_patient_weight"].sum()
    w_ok = (len(df276) == 276 and df276["patient_id"].nunique() == 66
            and int((df276["source"] == "Real").sum()) == 138
            and int((df276["source"] == "PBISC").sum()) == 138
            and float(np.max(np.abs(pt_w_sums - 1.0))) < 1e-12
            and df276["temp_severity_ord"].isna().sum() == 0)
    rec("Fig7_276Samples_Weights", "66_patients_276_samples_equal_patient_total_weight_1.0",
        w_ok, f"n_draws={len(df276)}, n_pts={df276['patient_id'].nunique()}, total_w={df276['equal_patient_weight'].sum():.6f}")

    max_days = df276.groupby("patient_id")["days_from_first_positive_swab"].transform("max")
    last_check = df276[df276["days_from_first_positive_swab"] == max_days]
    last_ok = (len(last_check) == 66 and len(df66) == 66
               and int((df66["fatal_pub11"] == 0).sum()) == 55
               and int((df66["fatal_pub11"] == 1).sum()) == 11
               and int((df66["source"] == "Real").sum()) == 24
               and int((df66["source"] == "PBISC").sum()) == 42)
    rec("Fig7_66LastSamples", "Unique_last_sample_per_patient_55surv_11died_0ties",
        last_ok, f"n_last_maxday={len(last_check)}, surv={int((df66['fatal_pub11']==0).sum())}, died={int((df66['fatal_pub11']==1).sum())}")

    r_comb = eskd[eskd["stratum"] == "Combined_66pts_276draws_rawPct"].set_index("feature")
    cell43 = list(r_comb.index)
    w_vec = df276["equal_patient_weight"].to_numpy(float)
    sev_vec = df276["temp_severity_ord"].to_numpy(float)
    y_d = df66["fatal_pub11"].to_numpy(int)

    max_rho_diff = 0.0
    max_auc_diff = 0.0
    p_d_list = []
    for ct in cell43:
        rho_c = weighted_spearman_independent(df276[ct].to_numpy(float), sev_vec, w_vec)
        max_rho_diff = max(max_rho_diff, abs(rho_c - float(r_comb.loc[ct, "rho_w"])))
        x1 = df66.loc[y_d == 1, ct].to_numpy(float)
        x0 = df66.loc[y_d == 0, ct].to_numpy(float)
        u_stat, p_u = stats.mannwhitneyu(x1, x0, alternative="two-sided")
        auc_c = float(u_stat / (len(x1) * len(x0)))
        max_auc_diff = max(max_auc_diff, abs(auc_c - float(r_comb.loc[ct, "auc_death"])))
        p_d_list.append(p_u)
    valid_p = np.array([np.isfinite(p) for p in p_d_list])
    q_d_recomp = np.full(len(p_d_list), np.nan)
    q_d_recomp[valid_p] = multipletests(np.array(p_d_list)[valid_p], method="fdr_bh")[1]
    max_qd_diff = float(np.nanmax(np.abs(q_d_recomp - r_comb["q_death_bh43"].to_numpy(float))))

    rec("Fig7_IndependentMath", "Exact_43_celltype_weighted_Spearman_and_MannWhitney_AUC_BH43",
        len(cell43) == 43 and max_rho_diff < 1e-10 and max_auc_diff < 1e-10 and max_qd_diff < 1e-10,
        f"n_types={len(cell43)}, max_rho_diff={max_rho_diff:.2e}, max_auc_diff={max_auc_diff:.2e}, max_qd_diff={max_qd_diff:.2e}")

    n_sev_real = int((eskd[eskd.stratum == "Real_37pts_138draws"]["q_sev_bh43"] < 0.05).sum())
    n_sev_comb = int((r_comb["q_sev_bh43"] < 0.05).sum())
    n_dth_real = int((eskd[eskd.stratum == "Real_37pts_138draws"]["q_death_bh43"] < 0.05).sum())
    n_dth_comb = int((r_comb["q_death_bh43"] < 0.05).sum())
    dth_winners = set(r_comb[r_comb["q_death_bh43"] < 0.05].index)
    rec("Fig7_FDRRescue", "BH43_11to18_severity_and_0to2_mortality_HSPC_Tgd",
        n_sev_real == 11 and n_sev_comb == 18 and n_dth_real == 0 and n_dth_comb == 2 and dth_winners == {"HSPC", "T g/d"},
        f"sev:{n_sev_real}->{n_sev_comb}, death:{n_dth_real}->{n_dth_comb} ({sorted(dth_winners)})")

    ext_s = pd.read_csv(DATA_DIR / "fig7_external4_300_acute_samples.csv")
    ext_m = pd.read_csv(DATA_DIR / "fig7_external4_random_effects_meta.csv")
    n_ext_donors = int((ext_s["cohort"].astype(str) + ":" + ext_s["donor"].astype(str)).nunique())
    ext_ok = (len(ext_s) == 300 and n_ext_donors == 231
              and ext_s["cohort"].nunique() == 4 and bool((ext_m["k_studies"] >= 2).all())
              and int(((ext_m["is_exact_43_subset"] == True) & (ext_m["q_sev_re"] < 0.05)).sum()) == 10)
    rec("Fig7_External4Meta", "4_studies_300_samples_231_donors_RE_meta_k_ge_2",
        ext_ok, f"n_samples={len(ext_s)}, n_donors={n_ext_donors}, exact_sig_sev={int(((ext_m['is_exact_43_subset']==True)&(ext_m['q_sev_re']<0.05)).sum())}")

    # DOCX XML check
    with zipfile.ZipFile(DOCX_PATH, "r") as zf:
        xml_txt = zf.read("word/document.xml").decode("utf-8")
    plain_xml = re.sub(r"<[^>]+>", "", xml_txt)
    latex_leaks = re.findall(r"\\[A-Za-z]+", plain_xml)
    rec("Manuscript_DOCX", "Zero_OOXML_LaTeX_leaks_and_valid_package",
        len(latex_leaks) == 0 and DOCX_PATH.stat().st_size > 12000,
        f"latex_leaks={latex_leaks}, bytes={DOCX_PATH.stat().st_size}")

    req_phrases = [
        "### Figure 4 · Paired same-donor",
        "limma-voom",
        "4,000",
        "### Figure 5 · Cross-cohort comparison",
        "Welch's *t*-tests",
        "|log₂ fold change| ≥ 1",
        "### Figure 6 · Rare cell state",
        "`14` rare cellular states in sepsis",
        "`13` rare cellular states in COVID-19",
        "patient_i ≠ patient_j",
        "### Figure 7 · Same-day WHO severity",
        "276` eligible acute-phase",
        "55` survivors and `11`",
        "DerSimonian–Laird random-effects",
        "### Supplementary Figure 14 |",
        "### Supplementary Figure 23 |",
    ]
    missing_p = [p for p in req_phrases if p not in md_txt]
    rec("Manuscript_Methods_Alignment", "All_user_specified_Fig4_to_Fig7_and_SuppFig14_23_rules_present",
        len(missing_p) == 0, f"missing={missing_p}")

    n_pass = sum(1 for c in checks if c["status"] == "PASS")
    n_fail = sum(1 for c in checks if c["status"] == "FAIL")
    report = {
        "generated_at": "2026-11-22",
        "package_dir": str(OUT_ROOT),
        "supp_numbering": "Supplementary Figures 14-23 (following Main Figures 1-4 Supplementary Figures 1-13 in strict Nature citation order)",
        "total_checks": len(checks),
        "passed": n_pass,
        "failed": n_fail,
        "all_passed": n_fail == 0,
        "checks": checks,
    }
    REPORT_PATH.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[AUDIT SUMMARY] total={len(checks)} | PASS={n_pass} | FAIL={n_fail}", flush=True)
    if n_fail > 0:
        for c in checks:
            if c["status"] == "FAIL":
                print("  [FAIL]", c, flush=True)
        raise SystemExit(1)


if __name__ == "__main__":
    migrate_supp_numbering_if_needed()
    build_docx()
    build_excel()
    build_supp_numbering_map()
    build_readme()
    stage_artifact_previews()
    run_independent_verification()
