import sys, nbformat
from nbclient import NotebookClient
p = sys.argv[1]; nb = nbformat.read(p, 4); NotebookClient(nb, timeout=900, kernel_name="python3", resources={"metadata": {"path": sys.argv[2]}}).execute(); nbformat.write(nb, p); print("executed ok")
