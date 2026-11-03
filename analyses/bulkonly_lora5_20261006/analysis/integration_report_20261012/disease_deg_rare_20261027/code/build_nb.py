import nbformat as nbf, os
nb=nbf.v4.new_notebook()
md=lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s))
code=lambda s: nb.cells.append(nbf.v4.new_code_cell(s))
md("""# PBISC floor read-out: disease DEG recomputation and rare-population expression (2026-10-27)

Inputs (all local, nothing regenerated):
- `ilisi_floor_20261021/data/pb/*.npz` — per-sample x six-type pseudobulk (raw count sums, full files), built by `pb_build.py`
- `ilisi_floor_20261021/data/big/*.h5ad` — 60k-cell balanced subsets (real half A = `ref`, real half B = `realU`, `floor`, `prev4k`)
- `ilisi_floor_20261021/results/big/*_transfer.csv.gz` — kNN-transferred labels for real B and floor

Contrasts: sepsis = adult sepsis vs adult healthy control; COVID = positive vs negative, aggregated to patient level.
DEG rule: Welch t-test on log2 CPM pseudobulk, BH FDR < 0.05 and |log2FC| >= 1. Groups with < 20 cells dropped.
COVID floor excludes the 5 adapter-training-donor samples.""")
for f,title in [("deg_lib.py","## Library: pseudobulk loading, design, DEG"),
                ("run_deg.py","## 1. Disease DEG, real vs floor vs prev K4000"),
                ("run_ceiling.py","## 2. Real-vs-real ceiling (split-half of the real samples, 20 reps)"),
                ("run_genes.py","## 3. Rank-matched metrics and top DEG lists"),
                ("run_named.py","## 4. Named genes from the manuscript / source papers"),
                ("rare_lib.py","## Library: rare-population definitions (same as f_rare_parent.py / big_rare.py)"),
                ("run_rare.py","## 5. Rare-population expression: profile correlation and marker programme"),
                ("run_rare2.py","## 6. Marker retention ratios and driver genes"),
                ("run_rare3.py","## 7. Power-matched comparison with a real-A-internal noise floor")]:
    md(title); code(open(f).read())
md("""## Outputs

| file | contents |
|---|---|
| `out/deg_summary.csv` | per cohort x six-type: DEG counts, logFC r, slope, recall/precision, top-100 overlap |
| `out/real_real_ceiling.csv` | split-half real-vs-real ceiling for the same metrics |
| `out/deg_rank2000_summary.csv` | sign agreement / slope restricted to real-expression rank > 2000, generated top-100 precision |
| `out/named_genes.csv` | per-gene real vs generated logFC for the manuscript's named genes |
| `out/rare_expression.csv`, `out/rare_expression_ratios.csv` | per rare subset: cell counts, profile r, marker score / detection, identity AUROC |
| `out/rare_gene_diff.csv.gz` | per-gene mean-expression difference vs real A inside each rare subset |
| `out/rare_power_matched.csv` | equal-n comparison with the real-A-internal noise floor |
""")
nbf.write(nb,"disease_deg_rare_20261027.ipynb")
print("written")
