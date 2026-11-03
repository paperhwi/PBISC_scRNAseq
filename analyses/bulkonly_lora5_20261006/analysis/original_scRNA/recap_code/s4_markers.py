"""Cluster marker tables (FindAllMarkers-like: only.pos, min.pct 0.1, log2FC > 0.25) + a canonical panel.
Usage: python s4_markers.py <run_name> <leiden key> [panel genes comma-separated]"""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/claude/work"); import store
W = "/home/claude/work/sepsis"
name, key = sys.argv[1], sys.argv[2]
panel = sys.argv[3].split(",") if len(sys.argv) > 3 else []
P = f"{W}/run_{name}"
st = store.Store(f"{W}/store")
cells = np.load(f"{P}/embed.npz")["cells"]; lab = np.load(f"{P}/leiden.npz")[key]
groups, mean, pct, n = store.group_stats(st, cells, lab)
np.savez(f"{P}/stats_{key}.npz", groups=np.array(groups), mean=mean, pct=pct, n=n)
gn = st.var.gene_symbol.astype(str).to_numpy()
mk = store.marker_table(groups, mean, pct, n, gn, top=12)
mk.to_csv(f"{P}/markers_{key}.csv", index=False)
pd.set_option("display.width", 320); pd.set_option("display.max_columns", 80); pd.set_option("display.max_rows", 200)
if panel:
    gi = [int(np.where(gn == x)[0][0]) for x in panel]
    M = pd.DataFrame(mean[:, gi], index=groups, columns=panel).round(1); M.insert(0, "n", n.astype(int))
    M["top_markers"] = mk.groupby("cluster")["gene"].apply(lambda s: ",".join(s.head(8))).reindex(groups).values
    print(M.to_string())
