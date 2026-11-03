import nbformat as nbf
from nbclient import NotebookClient
OUT = "/mnt/user-data/outputs/integration_report_20261012/macrophage_markers"
src = r'''import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
OUT = '%s'
res = pd.read_csv(f'{OUT}/macrophage_marker_detection_by_group.csv.gz'); M = pd.read_csv(f'{OUT}/macrophage_marker_panel_membership.csv')
GROUPS = ['Macrophage', 'CD14 mono', 'CD16 mono', 'DC', 'T/NK', 'B/Plasma', 'Other']
show = M[M.target.isin(['macrophage', 'myeloid', 'monocyte'])].gene.tolist()
panels = [('real', 'Sepsis · real'), ('k1000', 'Sepsis · generated K = 1,000'), ('k4000', 'Sepsis · generated K = 4,000')]
plt.rcParams['font.family'] = 'DejaVu Sans'
fig, axes = plt.subplots(1, 3, figsize=(13.5, 14), sharey=True)
for ax, (src, title) in zip(axes, panels):
    d = res[(res.cohort == 'sepsis') & (res.source == src)].pivot(index='gene', columns='group', values='det').reindex(index=show, columns=GROUPS)
    im = ax.imshow(d.to_numpy(), aspect='auto', cmap='Blues', vmin=0, vmax=1)
    for i in range(d.shape[0]):
        if d.iloc[i].isna().all(): ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(GROUPS), 1, color='#dddddd'))
    ax.set_xticks(range(len(GROUPS))); ax.set_xticklabels(GROUPS, rotation=55, ha='right', fontsize=10)
    ax.set_yticks(range(len(show))); ax.set_yticklabels(show, fontsize=9)
    n = res[(res.cohort == 'sepsis') & (res.source == src)].drop_duplicates('group').set_index('group').n_cells
    ax.set_title(f'{title}\nmacrophage-labelled cells: {int(n["Macrophage"]):,}', fontsize=11)
    ax.axhline((M.target == 'macrophage').sum() - 0.5, color='k', lw=1); ax.axhline(M.target.isin(['macrophage', 'myeloid']).sum() - 0.5, color='k', lw=1)
fig.colorbar(im, ax=axes, shrink=0.35, label='fraction of cells detected (count > 0)')
fig.text(0.01, 0.005, 'Grey rows: gene not in the 6k panel. Rows: macrophage genes / shared myeloid / monocyte contrast. Generated groups: format-aware classifier (route R3); real: Ye 2026 labels, up to 3,000 cells per group.', fontsize=8.5)
fig.savefig(f'{OUT}/macrophage_marker_heatmap_sepsis_real_k1000_k4000.png', dpi=150, bbox_inches='tight'); plt.show()''' % OUT
nb = nbf.v4.new_notebook(); nb.cells = [nbf.v4.new_markdown_cell('# Macrophage marker heatmap, sepsis only: real vs generated K = 1,000 vs K = 4,000\nRe-plot of `macrophage_markers.ipynb` results (`macrophage_marker_detection_by_group.csv.gz`).'), nbf.v4.new_code_cell(src)]
nb.metadata['kernelspec'] = {'name': 'python3', 'display_name': 'Python 3', 'language': 'python'}
NotebookClient(nb, timeout=600, kernel_name='python3').execute(); nbf.write(nb, f'{OUT}/macrophage_heatmap_sepsis_3panel.ipynb'); print('done')
