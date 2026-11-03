import nbformat as nbf
nb=nbf.v4.new_notebook(); md=lambda s: nb.cells.append(nbf.v4.new_markdown_cell(s)); code=lambda s: nb.cells.append(nbf.v4.new_code_cell(s))
md("""# Why are DCs under-produced in the PBISC floor data set? (2026-10-30)
DC-1 extent by group/subtype · DC-2 does the generated pseudobulk follow its own bulk input · DC-3/4 label-free cell-level pDC assembly ·
DC-5 size of the DC signal inside the bulk. Local data only; no new generation.""")
for f,t in [("dc/dc1_extent.py","## DC-1 Extent: group and subtype"),("dc/dc2_bulk_tracking.py","## DC-2 Bulk tracking (same samples)"),
            ("dc/dc3_cell_level.py","## DC-3 Label-free marker co-expression and lineage blending"),
            ("dc/dc4_lilra4_and_composition.py","## DC-4 LILRA4 assembly and composition tracking"),
            ("dc/dc5_signal_size.py","## DC-5 DC signal size inside the bulk")]:
    md(t); code(open(f).read())
nbf.write(nb,"dc/dc_underproduction_20261030.ipynb"); print("ok")
