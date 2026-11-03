"""Collect BayesPrism theta (fraction of reads from each type / state), convert to cell fractions by dividing by the
type's mean UMI per cell in the run's reference, and append to sepsis_deconv_props.csv.gz as 'BayesPrism' (no calibration) and 'BayesPrism+cal' (mixture corrected with the
same 5 calibration donors as WNNLS+cal; paired folds only). The bulk-only BayesPrism run was not done (see notebook)."""
import numpy as np, pandas as pd, os
O = "/home/claude/work/prop"


def lin6(x):
    return {"Plasma cells": "Plasma", "B cells": "B", "CD4+ T cells": "T", "CD8+ T cells": "T", "NKT cells": "T", "NK cells": "NK",
            "Monocytes": "Myeloid", "Macrophages": "Myeloid", "DC": "Myeloid"}.get(x, "Other")


D = pd.read_csv(f"{O}/sepsis_deconv_props.csv.gz"); D = D[~D.method.str.startswith("BayesPrism")]
rows = []
for run, st, meth in [("fold0", "paired", "BayesPrism"), ("fold1", "paired", "BayesPrism"), ("fold0cal", "paired", "BayesPrism+cal"), ("fold1cal", "paired", "BayesPrism+cal"), ("bulkonly", "bulkonly", "BayesPrism")]:
    ft, fs = f"{O}/bp/{run}_theta_type.csv", f"{O}/bp/{run}_theta_state.csv"
    if not (os.path.exists(ft) and os.path.exists(fs)): print("missing", run); continue
    sz = pd.read_csv(f"{O}/bp/{run}_size.csv")
    for f, lv, col in [(ft, "major", "major"), (fs, "fine", "fine")]:
        th = pd.read_csv(f, index_col=0)
        size = sz.groupby(col).umi.sum() / sz.groupby(col).n.sum()
        cf = th.div(size.reindex(th.columns).to_numpy(), axis=1); cf = cf.div(cf.sum(1), axis=0)
        long = cf.stack().rename("frac").reset_index(); long.columns = ["sample", "label", "frac"]
        long["level"] = lv; long["set"] = st; long["method"] = meth; long["draw"] = "-"; rows.append(long)
        if lv == "major":
            l6 = cf.T.groupby(cf.columns.map(lin6)).sum().T.stack().rename("frac").reset_index(); l6.columns = ["sample", "label", "frac"]
            l6["level"] = "lin6"; l6["set"] = st; l6["method"] = meth; l6["draw"] = "-"; rows.append(l6)
if rows:
    B = pd.concat(rows); B["grp"] = B["sample"].map(D.drop_duplicates("sample").set_index("sample").grp)
    D = pd.concat([D, B[D.columns]], ignore_index=True)
D.to_csv(f"{O}/sepsis_deconv_props.csv.gz", index=False); print(D.method.value_counts().to_dict())
