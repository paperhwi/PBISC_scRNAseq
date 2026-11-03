import sys,re,html,os,subprocess
d=sys.argv[1]; t=open(os.path.join(d,"downloads_index.html"),encoding="utf-8").read()
ps=[]
for u,n in re.findall(r'<a href="([^"]+)">([^<]+)</a>',t):
    u=html.unescape(u); assert "r2.cloudflarestorage.com" in u.split("/")[2]
    ps.append((n,subprocess.Popen(["curl","-sS","-f","--retry","3","-o",os.path.join(d,n),u])))
for n,p in ps: print(d,n,"rc",p.wait(),flush=True)
