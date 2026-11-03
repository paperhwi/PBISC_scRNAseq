import sys, nbformat
from nbclient import NotebookClient
p = sys.argv[1]; nb = nbformat.read(p, 4)
try:
    NotebookClient(nb, timeout=1400, kernel_name='python3', resources={'metadata': {'path': sys.argv[2]}}).execute()
finally:
    nbformat.write(nb, p)
print('errors:', sum(1 for c in nb.cells if c.cell_type == 'code' for o in c.get('outputs', []) if o.output_type == 'error'))
