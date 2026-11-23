# Figure 7 (First-Draw Anchored) Revision Pack (`2026-11-24`)

Self-contained, OS-independent data and instruction bundle for AI agents revising **Figure 7 (`Immune cell associations with COVID-19 severity and death in patients with ESKD`, 217 acute draws in 66 patients, first-draw mortality)** to resolve **Issues 1, 2, and 3**.

## Start Here
1. **Full AI-Agent Revision Specification**: [`Figure7_FirstDraw_Revision_Guide_for_AI_Agent_20261124.md`](Figure7_FirstDraw_Revision_Guide_for_AI_Agent_20261124.md)
2. **Zero-Dependency Verification & Table Printer**:
   ```bash
   python analyses/fig7_firstdraw_revision_pack_20261124/verify_and_recompute_pack.py
   ```
3. **Bundled Ready-to-Plot Data (`data/`)**:
   - `data/panel_k_ready_table.csv` — Ready-to-render 7-row external severity validation table (with `ESKD ρ_w`, `Aligned / k`, `ext_rho`, `ext_q23`, and proxy `*` flags)
   - `data/panel_l_ready_table.csv` — Ready-to-render external first- vs last-draw mortality table (including ESKD hits `CD8 CTL*`, `NK` + public-cohort `γδ T`, `MAIT`, `CD8 naïve T`, `CD16 mono`, `HSPC`)
   - `data/eskd_firstdraw_composite_and_cv_summary.csv` — Exact `Measured`, `PBISC inferred`, and `Combined` AUCs, 95% CIs, `adjWHO P`, and `50 × 5` CV AUCs for `CD8 CTL`, `NK`, `ISG − CTL − NK`, `ISG − CTL`, and `CD8 CTL ISG`
   - `data/eskd_63_exact_first_draws_CTL_NK_ISG.csv` — Exact 63-patient first-draw plotting coordinates (`ISG_minus_CTL_minus_NK`, `ISG_minus_CTL`)
   - `data/eskd_66_first_draws_43subsets_rankINT.csv` — Canonical 66-patient first-draw plotting coordinates (`Panels h, i`)
   - `data/eskd_217_acute_draws_43subsets_rankINT.csv` — Canonical 217-draw plotting coordinates (`Panels a, b, e–g2`)
   - `data/eskd_43subsets_stratified_severity_and_firstdraw_death.csv` — Complete 43-subset stratified ESKD statistics (`Panel c, d`)
   - `data/eskd_leakfree_50x5_cv_benchmark.csv` — Leak-free `50 × 5-fold` first-draw mortality CV benchmark
   - `data/external4_300_acute_samples.csv` — Sample-level 300 acute visits across 4 public COVID-19 scRNA-seq cohorts (`231` patients; `219` with documented outcome [`197` survivors, `22` deaths])
   - `data/external4_per_cohort_severity_and_first_last_death.csv` — Per-cohort external severity & 1st/last-draw mortality results
   - `data/external4_meta_severity_and_first_last_death_23subsets.csv` — Unified 23-subset DerSimonian–Laird random-effects meta-analysis
