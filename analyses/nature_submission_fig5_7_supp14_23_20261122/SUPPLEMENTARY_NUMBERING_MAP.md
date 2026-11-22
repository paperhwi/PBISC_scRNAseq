# Full-Paper Supplementary Figure Numbering Map (`Supplementary Figures 1–23`)

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
