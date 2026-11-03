import nbformat as nbf
nb=nbf.v4.new_notebook()
md=lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s))
code=lambda s: nb.cells.append(nbf.v4.new_code_cell(s))
md("""# Week-1 defensive analyses (2026-10-27)

A2 = is the module signal present in the matched bulk input?
C2 = are the lost modules in the low-expression regime? where do their transcripts sit?
A3 = power-corrected real-vs-real ceiling
D1 = patient-bootstrap CIs; D2 = permutation null and sepsis-vs-ICU hard negative

Inputs: `bulk_public/GSE279448` (sepsis, 207 bulk samples), `bulk_public/zenodo6497251_imperial_eskd`
(COVID, 336 bulk samples), the pseudobulk npz from `ilisi_floor_20261021/data/pb`, and the
60k balanced cell sets from `ilisi_floor_20261021/data/big`. No model change, no new generation.""")
for f,t in [("bulk_prep.py","## Bulk matrices and Ensembl↔symbol mapping"),
            ("a2_bulk.py","## A2 — module signal in the bulk vs in floor"),
            ("c2_rank.py","## C2 — expression rank, detection, and UMI concentration"),
            ("c2b_load.py","## C2b — positive-cell fraction vs transcript load (floor and prev K4000)"),
            ("a3_d1_d2.py","## A3 / D1 / D2 — power-corrected ceiling, bootstrap CIs, negative controls"),
            ("fig_week1.py","## Figures")]:
    md(t); code(open(f).read())
md("""## Outputs
`out2/a2_module_bulk.csv`, `out2/a2_bulk_background.csv`, `out2/sepsis_bulk_deg.csv`, `out2/covid_bulk_deg.csv`,
`out2/c2_module_rank.csv`, `out2/c2_concentration.csv`, `out2/c2b_module_load.csv`,
`out2/a3_power_matched_ceiling.csv`, `out2/d1_bootstrap_ci.csv`, `out2/d2_permutation_null.csv`, `out2/d2b_sepsis_vs_icu.csv`""")
nbf.write(nb,"week1_defensive_20261027.ipynb"); print("ok")
