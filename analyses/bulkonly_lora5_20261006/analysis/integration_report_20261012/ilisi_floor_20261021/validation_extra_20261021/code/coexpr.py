"""Within-lineage gene-gene co-expression (2026-10-21). Usage: python coexpr.py <coh>
Per lineage, genes = top 300 by real-A mean (high band) and real-A ranks 1,001-1,300 (low band); Pearson correlation of log1p CP10k over
5,000 cells per source. Similarity to real A = Pearson r of the upper triangles. Controls: real B (other patients), real B depth-matched to
floor / prev, real B with every gene permuted independently across cells (destroys co-expression, keeps marginals)."""
import sys, numpy as np, pandas as pd
from common import *
coh = sys.argv[1]; D = bundle(coh); rng = np.random.default_rng(1); rows = []
D["shuffle"] = D["realU"]
for lin in ("T", "NK", "B", "Myeloid"):
    Y = {}
    for s, a in D.items():
        idx = np.where(a.obs.lineage.to_numpy() == lin)[0]; idx = rng.choice(idx, min(5000, len(idx)), replace=False)
        Y[s] = L.lognorm(a.X.tocsr()[idx]).toarray().astype(np.float32)
    Y["shuffle"] = np.column_stack([rng.permutation(Y["shuffle"][:, j]) for j in range(Y["shuffle"].shape[1])])
    order = np.argsort(-Y["ref"].mean(0))
    for band, gi in (("top 300 genes", order[:300]), ("ranks 1,001-1,300", order[1000:1300])):
        def cm(M):
            X = M[:, gi]; X = X - X.mean(0); s = X.std(0); ok = s > 0; C = np.full((len(gi), len(gi)), np.nan)
            C[np.ix_(ok, ok)] = np.corrcoef(X[:, ok].T); return C
        R = cm(Y["ref"]); iu = np.triu_indices(len(gi), 1)
        for s in ("realU", "realT1f", "realT1p", "shuffle", "floor", "prev4k"):
            S = cm(Y[s]); a, b = R[iu], S[iu]; ok = np.isfinite(a) & np.isfinite(b)
            top = np.argsort(-np.abs(np.nan_to_num(a)))[:500]; okt = ok[top]
            rows.append({"cohort": coh, "lineage": lin, "genes": band, "source": s, "corr_of_corr": float(np.corrcoef(a[ok], b[ok])[0, 1]),
                         "mean_abs_corr_src": float(np.nanmean(np.abs(b))), "mean_abs_corr_real": float(np.nanmean(np.abs(a))),
                         "top500_pairs_sign_agree": float((np.sign(a[top][okt]) == np.sign(b[top][okt])).mean()),
                         "frac_pairs_defined": float(ok.mean())})
C = pd.DataFrame(rows); C.to_csv(f"{OUT}/{coh}_coexpression.csv", index=False)
pd.set_option("display.width", 250)
print(C.pivot_table(index=["genes", "lineage"], columns="source", values="corr_of_corr").round(3).to_string())
print(C.pivot_table(index=["genes", "lineage"], columns="source", values="mean_abs_corr_src").round(3).to_string())
