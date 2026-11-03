import sys, re, html, os, subprocess
root, cohort, secs = sys.argv[1], sys.argv[2], sys.argv[3]
d = os.path.join(root, cohort)
t = open(os.path.join(d, "downloads_index.html"), encoding="utf-8").read()
for url, name in re.findall(r'<a href="([^"]+)">([^<]+)</a>', t):
    if name.endswith(".h5ad"):
        subprocess.run(["curl","-sS","-C","-","--max-time",secs,"-o",os.path.join(d,name),html.unescape(url)])
        print(cohort, os.path.getsize(os.path.join(d,name)), flush=True)
