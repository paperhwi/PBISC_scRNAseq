import sys, re, html, os, subprocess, urllib.request
from urllib.parse import urlparse
root, cohort, index_url = sys.argv[1:4]
d = os.path.join(root, cohort)
t = urllib.request.urlopen(index_url, timeout=60).read().decode()
open(os.path.join(d, "downloads_index.html"), "w", encoding="utf-8").write(t)
links = re.findall(r'<a href="([^"]+)">([^<]+)</a>', t)
procs = []
for url, name in links:
    url = html.unescape(url)
    assert urlparse(url).netloc.endswith("r2.cloudflarestorage.com")
    procs.append((name, subprocess.Popen(["curl","-sS","-f","--retry","3","-o",os.path.join(d,name),url])))
for name, p in procs:
    print(name, "rc", p.wait(), flush=True)
