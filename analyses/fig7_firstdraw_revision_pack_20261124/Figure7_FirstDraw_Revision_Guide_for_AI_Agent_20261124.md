# Figure 7 (First-Draw Anchored) Revision Specification for AI Agent

> **Repository**: `https://github.com/paperhwi/PBISC_scRNAseq`  
> **Self-Contained Pack Path**: `analyses/fig7_firstdraw_revision_pack_20261124/`  
> **Target Figure**: `Figure 7. Immune cell associations with COVID-19 severity and death in patients with ESKD` (`217 acute draws in 66 patients`, `First draw per patient: 55 survived, 11 died`)  
> **Goal**: Resolve **Consistency Issues 1, 2, and 3** using the pre-verified, zero-dependency CSV tables bundled in `analyses/fig7_firstdraw_revision_pack_20261124/data/`.

---

## 0. Quick-Start GitHub Environment & Frozen Guardrails

### 0.1 Clone, Install & Verify in 3 Commands (OS-Independent)
All data tables and verification code are self-contained inside this repository with **zero hardcoded local paths**:

```bash
git clone https://github.com/paperhwi/PBISC_scRNAseq.git
cd PBISC_scRNAseq
pip install -r requirements.txt
python analyses/fig7_firstdraw_revision_pack_20261124/verify_and_recompute_pack.py
```

Running `verify_and_recompute_pack.py` asserts `max_abs_diff < 1e-12` across every ESKD and external 4-cohort statistic and prints all ready-to-plot tables to stdout. Follow repository figure-making standards in [`AGENTS.md`](../../AGENTS.md) (`width = 180.0 mm`, `height ≤ 240.0 mm`, sans-serif `Arial`/`Helvetica` `≥ 5.0 pt`, top/right spines off, export `.pdf` + `.svg` + `300–600 DPI .png`, visual collision inspection).

### 0.2 Bundled Data Files (`analyses/fig7_firstdraw_revision_pack_20261124/data/`)
Load numbers **directly from these CSVs**—do **not** approximate or recompute with ad-hoc filters:

| Relative File Path (inside `analyses/fig7_firstdraw_revision_pack_20261124/`) | Purpose / Target Panels |
| :--- | :--- |
| `data/eskd_217_acute_draws_43subsets_rankINT.csv` | Sample-level matrix (`217` acute draws, `66` patients): draw-day WHO (`WHO_temp_severity`), raw `%`, and within-source normal scores (`*__z217`) for **Panels `a, b, e–g2`** |
| `data/eskd_66_first_draws_43subsets_rankINT.csv` | Patient-level canonical first-draw matrix (`66` patients: `55` survived, `11` died; `*__zfd`) for **Panels `a, b, h, i`** |
| `data/eskd_63_exact_first_draws_CTL_NK_ISG.csv` | Exact canonical first-draw matrix with `≥ 5` `CD8 CTL` cells (`63` patients: `52` survived, `11` died), including pre-computed `ISG_minus_CTL_minus_NK` (`(isg_zfd - ctl_zfd - nk_zfd)/3`) and `ISG_minus_CTL` (`(isg_zfd - ctl_zfd)/2`) for **Panel `j`** |
| `data/eskd_43subsets_stratified_severity_and_firstdraw_death.csv` | Stratified (`Real`, `PBISC new`, `Combined`) `GEE +wk` slopes, equal-patient-weight `Spearman ρ_w` (`sd_equal_wt_rho`), first-draw `AUC`/`q43`/`adjWHO_p`, and last-draw `AUC` across all `43` subsets for **Panels `c, d`** |
| `data/eskd_firstdraw_composite_and_cv_summary.csv` | Exact `Measured`, `PBISC inferred`, and `Combined` AUCs, 95% CIs, `P`, `q43`, `adjWHO_p`, and `50 × 5-fold` CV AUCs for `CD8 CTL`, `NK`, `ISG − CTL − NK`, `ISG − CTL`, and `CD8 CTL ISG` for **Panel `d` & Issue 2** |
| `data/eskd_leakfree_50x5_cv_benchmark.csv` | Leak-free `50 × 5-fold` stratified first-draw mortality cross-validation results (`N = 66` and `n = 63`) |
| `data/external4_300_acute_samples.csv` | Sample-level external validation matrix (`300` acute visits across `231` patients; `219` with documented outcome [`197` survivors, `22` deaths]) |
| `data/external4_per_cohort_severity_and_first_last_death.csv` | Per-cohort (`Stephenson 2021`, `COMBAT 2022`, `van der Wijst 2021`, `Liu 2021`) draw-day `ρ_w`, first-draw mortality `AUC`, and last-draw mortality `AUC` |
| `data/external4_meta_severity_and_first_last_death_23subsets.csv` | Complete DerSimonian–Laird (`DL`) random-effects meta-analysis across all `23` exact subsets (`k ≥ 2`) + `B_ASC_total`, paired side-by-side with ESKD statistics and public-atlas mapping notes |
| **`data/panel_k_ready_table.csv`** | **Ready-to-render 7-row table for revised `Panel k`** (includes `display_label` with `*`, `eskd_gee_wk_beta`, `eskd_equal_wt_rho`, `aligned_over_k`, `n_patients_ext`, `ext_rho`, `ext_ci_lo`, `ext_ci_hi`, `ext_q23`, `external_evidence`) |
| **`data/panel_l_ready_table.csv`** | **Ready-to-render table for revised `Panel l`** (includes `CD8 CTL*`, `NK`, `γδ T`, `MAIT`, `CD8 naïve T`, `CD16 mono`, `HSPC` with per-cohort 1st-draw AUCs and RE 1st- vs last-draw AUCs/`q23`) |

### 0.3 Frozen Baseline (DO NOT ALTER)
All counts and raw plotting coordinates in **Panels `a`, `b`, `e–g2`, `h`, and `i`** are 100% verified and must remain untouched:
- ESKD cohort: `37` `Measured` patients (`138` draws; `32` survived, `5` died) + `29` `PBISC` patients (`79` draws; `23` survived, `6` died) = **`66` patients (`217` acute draws; `55` survived, `11` died)**; `2` CLL patients excluded (`C116`, `C139`).
- Draw-day WHO grades (`217` draws): `Mild` `77 + 49 = 126`, `Moderate` `39 + 19 = 58`, `Severe` `15 + 6 = 21`, `Critical` `7 + 5 = 12`.
- Panel `j` sample size (`n = 63`): `52` survived (`30` Measured + `22` PBISC), `11` died (`5` Measured + `6` PBISC), excluding the `3` surviving patients with `< 5` `CD8 CTL` cells on draw 1 (`C113`, `C129`, `C140`).

---

## 1. Issue 1: Resolve Target-Switching Between ESKD First-Draw Mortality (`Panels d, h–j`) and External Mortality (`Panels k, l`)

### 1.1 Problem Diagnosis
1. **Omission of ESKD's primary mortality hits**: `Panel d` identifies **`CD8 CTL`** (`AUC = 0.213, q43 = 0.070`) and **`NK`** (`AUC = 0.217, q43 = 0.070`) as ESKD's two primary first-draw mortality predictors, and `Panel k` is titled *"External validation of severity and mortality associations in four COVID-19 cohorts"*. Yet **`Panel l` completely omits both `CD8 CTL` and `NK`**.
2. **Substitution with ESKD-null subsets**: `Panel l` displays **`CD8 naïve T`**, **`MAIT`**, and **`CD16 mono`**, none of which passed significance on ESKD first draws. In particular, **`CD8 naïve T`** (`ext 1st-draw q23 = 3.4e-5*`) is **strictly null in ESKD first draw (`Combined AUC = 0.488, P = 0.911, q43 = 0.973`)** with opposite directions across arms (`Measured 0.244` vs `PBISC 0.710`, $P_{\text{het}} = 0.017$). Showing only `CD8 naïve T`, `MAIT`, and `CD16 mono` under an *"External validation"* heading creates a target-switching inconsistency.

### 1.2 Exact Ground-Truth Table (`data/panel_l_ready_table.csv`)
External cohort `n / deaths` with documented outcomes (`219` patients: `197` survivors, `22` deaths): `Stephenson 2021 = 73/7`, `COMBAT 2022 = 75/7`, `van der Wijst 2021 = 42/5`, `Liu 2021 = 29/3`.

| Block | `display_label` | ESKD 1st AUC (`P` / `q43`) | ESKD Last AUC (`P`) | Ext. 1st-Draw Cohort AUCs (`Steph / COMBAT / vdW / Liu`) | Ext. `RE_AUC_1st` (`95% CI`) | Ext. `q (1st)` (`RE_q23_1st`) | Ext. `RE_AUC_last` (`95% CI`) | Ext. `q (last)` (`RE_q23_last`) | Interpretation |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `ESKD_1st_draw_targets` | **`CD8 CTL*`** | **`0.213`** (`0.0029` / `0.070`) | `0.236` (`0.0062`) | `0.50 / 0.46 / 0.46 / 0.49` | **`0.48`** (`0.35–0.60`) | `0.926` | `0.42` (`0.30–0.54`) | `0.482` | ESKD hit #1; coarse public `CD8 memory` labels (`*`) + ESKD uraemic CTL vulnerability |
| `ESKD_1st_draw_targets` | **`NK`** | **`0.217`** (`0.0033` / `0.070`) | `0.357` (`0.139`) | `0.67 / 0.30 / 0.46 / 0.63` | **`0.50`** (`0.31–0.69`) | `0.965` | `0.55` (`0.36–0.75`) | `0.926` | ESKD hit #2; ESKD-specific baseline cytotoxic NK attrition |
| `Public_1st_vs_last_shift` | **`γδ T`** | `0.315` (`0.055` / `0.366`) | **`0.198`** (`0.0017`) | `0.38 / 0.29 / 0.25 / 0.00` | **`0.22`** (`0.01–0.42`) | `0.051` (`P=0.0067`) | `0.26` (`0.02–0.49`) | `0.314` | **Concordant across ESKD & external** at both 1st and last draw |
| `Public_1st_vs_last_shift` | **`MAIT`** | `0.344` (`0.106` / `0.391`) | `0.331` (`0.081`) | `0.30 / 0.29 / 0.31 / 0.05` | **`0.22`** (`0.06–0.38`) | **`0.006*`** | `0.27` (`0.05–0.49`) | `0.314` | Same-direction innate-like T loss (`AUC = 0.34` in ESKD; FDR-sig. externally) |
| `Public_1st_vs_last_shift` | **`CD8 naïve T`** | `0.488` (`0.911` / `0.973`) | `0.586` (`0.376`) | `0.28 / 0.29 / 0.21 / 0.32` | **`0.27`** (`0.17–0.36`) | **`3.4e-5*`** | `0.37` (`0.22–0.52`) | `0.427` | Non-dialysis early marker; null in elderly/uraemic ESKD (`AUC = 0.49`) |
| `Public_1st_vs_last_shift` | **`CD16 mono`** | `0.476` (`0.809` / `0.916`) | **`0.269`** (`0.017`) | `0.35 / 0.55 / 0.66 / 0.42` | **`0.48`** (`0.34–0.62`) | `0.938` | **`0.32`** (`0.21–0.43`) | **`0.020*`** | **Concordant timing shift** (`1st ≈ 0.48` null → `Last ≈ 0.27–0.32` depleted) |

### 1.3 Actionable Modification Steps (`Panels k` and `l`)
1. **Update `Panel k` and `Panel l` titles to separate severity validation from cross-cohort mortality comparison**:
   - `Panel k` title → **`"External validation of draw-day severity and cross-cohort mortality comparison in four COVID-19 cohorts"`**
   - `Panel l` title → **`"External mortality at the first versus last draw: ESKD cytotoxic targets versus general COVID-19 cohorts"`**
2. **Populate `Panel l` directly from `data/panel_l_ready_table.csv` so ESKD's primary mortality targets are explicitly shown**:
   - **Recommended (Option 1A — 5 or 6 rows in `Panel l`, divided by a thin horizontal line between `ESKD_1st_draw_targets` and `Public_1st_vs_last_shift`)**:
     - Top block (`ESKD 1st-draw targets`): **`CD8 CTL* (4 cohorts)`** and **`NK (4 cohorts)`**.
     - Bottom block (`Public 1st- vs last-draw shift`): **`γδ T (4 cohorts)`**, **`MAIT (4 cohorts)`**, **`CD8 naïve T (4 cohorts)`**, and **`CD16 mono (4 cohorts)`** (if limited to 5 rows, include `CD8 CTL*`, `NK`, `γδ T` [or `MAIT`], `CD8 naïve T`, `CD16 mono`).
   - **Compact Alternative (Option 1B — if strictly keeping 4 rows in `Panel l`)**:
     - Plot **`CD8 CTL*`**, **`γδ T`** (or **`MAIT`**), **`CD8 naïve T`**, and **`CD16 mono`**, and place a direct note below `Panel l`:
       `"ESKD 1st-draw hits in public cohorts: CD8 CTL* 1st AUC = 0.48 (q23 = 0.93), NK 1st AUC = 0.50 (q23 = 0.97); *broad CD8 memory/TE labels in public atlases."`

---

## 2. Issue 2: Fix Sign Notation, Visual Separation, and CV Footnote Alignment for `CTL–NK–ISG` (`Panels d, j`)

### 2.1 Problem Diagnosis
1. **Misleading algebraic sign label**: In `data/eskd_63_exact_first_draws_CTL_NK_ISG.csv`, the plotted score is `(isg_zfd - ctl_zfd - nk_zfd) / 3` (`ISG_minus_CTL_minus_NK`), so **higher values predict death (`Combined AUC = 0.841 > 0.5`)**. Writing `CTL–NK–ISG` reads algebraically as `CTL − NK − ISG` (`+CTL` and `−ISG`, the exact opposite sign).
2. **Unseparated test status & direction flip in `Panel d`**: Rows 1–2 (`CD8 CTL`, `NK`) are single-fraction hits from the 43-subset screen (`BH q43 = 0.070 < 0.10`, **lower** in fatal, `AUC < 0.5`), whereas Row 3 is a **post-hoc 3-feature composite** (`q = —`, `n = 63`, **higher** in fatal, `AUC > 0.5`).
3. **Mismatch between plotted 3-feature composite (`CTL + NK + ISG`) and footnote 2-feature CV model (`CTL + ISG`)**:
   - `Panels d` and `j` plot the **3-feature** composite including `NK` (`in-sample AUC = 0.841`),
   - Yet the `Panel d` footnote cites `"CV AUC: clinical 0.694 → clinical + CTL + ISG 0.858"` (**excluding `NK`**).
   - From `data/eskd_leakfree_50x5_cv_benchmark.csv`, the leak-free `50 × 5-fold` CV model **including `NK`** (`Exact 1st-dr (63) + CTL + NK + ISG`) is **`0.846` (`0.8456`)**, whereas excluding `NK` (`+ CTL + ISG`) is **`0.858` (`0.8579`)**.

### 2.2 Exact Ground-Truth Table (`data/eskd_firstdraw_composite_and_cv_summary.csv`)

| `feature_label` | CSV Column (`eskd_63_exact_first_draws_CTL_NK_ISG.csv`) | Measured AUC (`95% CI`, `P`) | PBISC Inferred AUC (`95% CI`, `P`) | Combined AUC (`95% CI`, `P`) | `q43` | WHO-Adj. `P` (`adjWHO_p`) | Leak-Free `50 × 5` CV (`base → ext`) |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`CD8 CTL`** (`N = 66`) | `T CD8 CTL__zfd` (`66`-pt file) | **`0.206`** (`0.063–0.388`, `0.036`) | **`0.210`** (`0.058–0.399`, `0.031`) | **`0.213`** (`0.111–0.332`, `0.0029`) | **`0.070`** | **`0.033`** | `0.688 → 0.826` (`N=66`) |
| **`NK`** (`N = 66`) | `NK__zfd` (`66`-pt file) | **`0.256`** (`0.050–0.519`, `0.088`) | **`0.181`** (`0.029–0.377`, `0.016`) | **`0.217`** (`0.088–0.369`, `0.0033`) | **`0.070`** | **`0.030`** | `0.688 → 0.827` (`+CTL+NK`, `N=66`) |
| **`ISG − CTL − NK`** (`n = 63`) | `ISG_minus_CTL_minus_NK` | **`0.887`** (`0.727–1.000`, `0.0039`) | **`0.803`** (`0.591–0.970`, `0.024`) | **`0.841`** (`0.712–0.946`, `4.3e-4`) | `—` | **`0.013`** | **`0.694 → 0.846`** (`+CTL+NK+ISG`, `n=63`) |
| **`ISG − CTL`** (`n = 63`, alt.) | `ISG_minus_CTL` | **`0.887`** (`0.747–0.987`, `0.0039`) | **`0.742`** (`0.515–0.924`, `0.078`) | **`0.809`** (`0.675–0.921`, `0.0014`) | `—` | **`0.013`** | **`0.694 → 0.858`** (`+CTL+ISG`, `n=63`) |
| **`CD8 CTL ISG`** (`n = 63`, ref.) | `isg_zfd` | `0.760` (`0.567–0.913`, `0.069`) | `0.621` (`0.348–0.856`, `0.395`) | `0.685` (`0.521–0.834`, `0.056`) | `—` | `0.090` | `0.694 → 0.763` (`+ISG`, `n=63`) |

### 2.3 Actionable Modification Steps (`Panels d` and `j`)
1. **Fix the composite label so the signs match the formula**:
   - Change `Panel j` header from `"j CTL–NK–ISG score"` to **`"j ISG − CTL − NK score"`** (plot column `ISG_minus_CTL_minus_NK` from `data/eskd_63_exact_first_draws_CTL_NK_ISG.csv`).
   - Change `Panel d` row 3 y-tick label from `"CTL–NK–ISG"` to **`"ISG − CTL − NK"`** (with subtle sub-label `"composite, n = 63"`).
2. **Add a horizontal divider and explicit direction/FDR cues in `Panel d`**:
   - Draw a thin dashed horizontal line between Row 2 (`NK`) and Row 3 (`ISG − CTL − NK`).
   - Update `Panel d` subtitle to:
     **`"First draw · 43-type screen (q43 < 0.10) · WHO-adj. P: CTL 0.033, NK 0.030, composite 0.013."`**
   - Add subtle directional cues along the x-axis (`"← Lower in fatal"` left of `0.5`, `"Higher in fatal →"` right of `0.5`).
3. **Align the CV footnote with the plotted composite**:
   - **Recommended (Option 2A — keep 3-feature `ISG − CTL − NK` in `j` and `d`)**: Set the `Panel d` bottom note to:
     **`"95% CI · q43 across 43 fractions | 50×5 CV AUC (n = 63): clinical 0.694 → + CTL + NK + ISG 0.846 (+ CTL + ISG 0.858)"`**
   - **Alternative (Option 2B — switch `j` and `d` to 2-feature `ISG − CTL`)**: Plot column `ISG_minus_CTL` (`Combined AUC = 0.809, Measured = 0.887, PBISC = 0.742`), label it **`"ISG − CTL"`**, and set the bottom note to `"95% CI · q43 across 43 fractions | 50×5 CV AUC (n = 63): clinical 0.694 → + CTL + ISG 0.858"`.

---

## 3. Issue 3: Harmonize Severity Scales (`Panel c` vs `Panel k`), Multiplicity (`q43` vs `q23`), and Proxy Cluster Labels (`*`)

### 3.1 Problem Diagnosis
1. **Severity scale asymmetry (`GEE +wk` slope vs `Spearman ρ`)**:
   - `Panel c` (ESKD `217` draws) reports week-adjusted `GEE +wk` linear slope (`Normal-score change per WHO grade`, `BH43`), whereas `Panel k` (External `300` draws) reports unadjusted equal-patient-weight (`w_ij = 1 / n_i`) `Spearman ρ` (`BH23`).
   - *Why external `Spearman ρ` must be retained (DO NOT fit raw linear grade slopes across external cohorts)*: In `data/external4_300_acute_samples.csv`, external studies use non-identical ordinal scales (`Stephenson`: `{1, 4, 5, 6, 7}`, `COMBAT`: `{2, 4, 5, 7}`, `van der Wijst`: `{4, 5, 6, 7}`, `Liu`: `{4.5, 5.5, 7.0}`) and `days_onset` is missing in `14 / 300` visits. Weighted rank correlation `Spearman ρ` (`w_ij = 1 / n_i`) is strictly invariant to cohort-specific ordinal spacing, and across all `43` ESKD subsets `GEE +wk` $\beta$ and equal-patient-weight `ρ_w` correlate at **`r = 0.892`**. However, this rationale—and ESKD's matching `ρ_w` values—must be shown explicitly.
2. **Unmarked public-atlas proxy annotations (`Dividing ASC*`, `T CD8 CM*`)**:
   - `Stephenson 2021` has no `CD8 CM` cluster and uses **`CD8.EM` (effector memory)** as a proxy in the `2 / 2` pooled `T CD8 CM` estimate (`COMBAT`-only direct `CD8.TCM/CCL5` has `ρ_w = −0.224, P = 0.035`). `Dividing ASC` maps to `Plasmablast` (`Stephenson`) and `PB_Prolif` (`van der Wijst`).
3. **Different BH denominators (`q43` vs `q23`)**:
   - `q43` adjusts over all `43` ESKD cell types; `q23` adjusts over all `23` exact cell types mapped in `k ≥ 2` external cohorts.

### 3.2 Exact Ground-Truth Table (`data/panel_k_ready_table.csv`)

| `display_label` | `feature` | ESKD `GEE +wk` $\beta$ (`q43`) [`Panel c`] | ESKD Equal-Weight `ρ_w` (`q43`) [Same Estimator as `k`] | `aligned_over_k` | `n_patients_ext` | External `ext_rho` (`95% CI`) | `ext_q23` | `external_evidence` |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **`Dividing ASC*`** | `B_ASC_dividing` | **`+0.422`** (`1.8e-8`) | **`+0.423`** (`4.5e-6`) | `2 / 2` | `127` | **`+0.26`** (`+0.08, +0.42`) | `0.013` | `Supported` |
| **`HSPC`** | `HSPC` | **`+0.333`** (`0.00013`) | **`+0.352`** (`8.0e-4`) | `3 / 3` | `202` | **`+0.40`** (`+0.21, +0.57`) | `6.6e-4` | `Supported` |
| **`MAIT`** | `MAIT` | **`−0.164`** (`0.030`) | **`−0.200`** (`0.078`, `P=0.034`) | `4 / 4` | `231` | **`−0.29`** (`−0.45, −0.12`) | `0.004` | `Supported` |
| **`CD56-bright NK`** | `NK CD56 bright` | **`−0.276`** (`8.2e-4`) | **`−0.389`** (`1.8e-5`) | `3 / 4` | `231` | **`−0.28`** (`−0.46, −0.08`) | `0.016` | `Supported` |
| **`γδ T`** | `T g/d` | **`−0.196`** (`0.003`) | **`−0.251`** (`0.017`) | `4 / 4` | `231` | **`−0.27`** (`−0.39, −0.14`) | `6.6e-4` | `Supported` |
| **`CD16 mono C1`** | `CD16mono_C1` | **`+0.247`** (`0.004`) | **`+0.241`** (`0.017`) | `1 / 1` | `85` | **`+0.09`** (`−0.12, +0.29`) | `—` | `One study; not pooled` |
| **`T CD8 CM*`** | `T CD8 CM` | **`−0.231`** (`0.003`) | **`−0.426`** (`1.2e-6`) | `2 / 2` | `160` | **`−0.29`** (`−0.43, −0.13`) | `0.002` | `Supported` |

### 3.3 Actionable Modification Steps (`Panels c` and `k`)
1. **Render `Panel k` directly from `data/panel_k_ready_table.csv`, preserving asterisks on `Dividing ASC*` and `T CD8 CM*` (and `CD8 CTL*` in `Panel l`)**, and place the mapping note below `Panel l` (or in the legend with `*` on the figure):
   **`"* Public-atlas mapping: Dividing ASC = Plasmablast (Stephenson), PB_Prolif (van der Wijst); T CD8 CM = CD8.TCM/CCL5 (COMBAT; single-cohort ρ = −0.22, P = 0.035) + CD8.EM proxy (Stephenson); CD8 CTL = CD8.TE (Stephenson), CD8.TEMRA/TEFF/TEM (COMBAT), T8_Mem (van der Wijst), CD8 memory (Liu)."`**
2. **Bridge ESKD and External severity scales explicitly**:
   - **Recommended (Option 3A — Add compact `ESKD ρ` column to `Panel k` + clarify subtitle)**:
     - Add column **`ESKD ρ`** (`+0.42`, `+0.35`, `−0.20`, `−0.39`, `−0.25`, `+0.24`, `−0.43` from `eskd_equal_wt_rho` in `data/panel_k_ready_table.csv`) immediately left of `Aligned / k` so ESKD and external cohorts can be compared on the **exact same `ρ_w` estimator**, while `Panel c` retains swab-week-adjusted `GEE +wk` slopes.
     - Update `Panel k` table subtitle to:
       **`"300 draws in 231 patients · Equal-patient-weight (1/n_i) draw-day Spearman ρ across study-specific ordinal scales; BH q23 across all 23 pooled cell types (k ≥ 2). Of 18 ESKD findings, 11 matched ≥2 cohorts, 5 matched one and 2 matched none."`**
   - **Minimal-Space Alternative (Option 3B)**: Keep the updated `Panel k` subtitle above and state `"ESKD GEE +wk slope and equal-weight ρ correlate at r = 0.89 across 43 types"` in the subtitle/caption if horizontal width cannot accommodate an `ESKD ρ` column.

---

## 4. Final Verification Checklist Before Commit

1. Run `python analyses/fig7_firstdraw_revision_pack_20261124/verify_and_recompute_pack.py` and confirm `[PASS]`.
2. Confirm **Issue 1**: `CD8 CTL*` (`1st RE AUC = 0.48, q23 = 0.926`) and `NK` (`1st RE AUC = 0.50, q23 = 0.965`) are explicitly displayed or annotated in `Panel l`, and `Panels k, l` headers distinguish severity validation from cross-cohort mortality comparison.
3. Confirm **Issue 2**: `CTL–NK–ISG` is renamed **`ISG − CTL − NK`**, visually separated from the 43-subset screen (`q43 < 0.10`), and paired with the matching CV footnote (`clinical 0.694 → + CTL + NK + ISG 0.846 (+ CTL + ISG 0.858)`).
4. Confirm **Issue 3**: `Dividing ASC*`, `T CD8 CM*`, and `CD8 CTL*` carry `*` with public-atlas cluster definitions, and both `q43` and `q23` are explicitly defined.
5. Visually inspect the exported PNG (`300–600 DPI`) to ensure zero overlapping text, zero clipped error bars, and adherence to [`AGENTS.md`](../../AGENTS.md).
