from pathlib import Path
p=Path('paper_claim_analysis/YE_MARKER_MANIFEST.md')
s=p.read_text(encoding='utf-8-sig')
s=s.replace('`score_claim_expression_original.py` reads the fixed definitions, so the original cohort and later K comparisons use the full, unchanged author marker lists.', '`ye_supplementary_table6_expression_coverage.csv` reports marker availability in the original sepsis expression matrix: 17,733/17,743 rows are present. All seven claim panels are present except `CYB561D2.1`, one of 307 NK-CCL4 markers (306/307 present). `score_claim_expression_original.py` reads the fixed definitions, so the original cohort and later K comparisons use the full, unchanged author marker lists, scoring only genes available in the expression matrix.')
p.write_text(s,encoding='utf-8')
