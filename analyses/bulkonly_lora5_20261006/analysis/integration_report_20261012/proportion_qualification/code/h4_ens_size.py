"""Ensemble-size curve from stored paired predictions: mean of k eligible draws (all subsets of size k), per-sample L1
vs the frozen scANVI six-type truth, sepsis and COVID common-test samples; equal-information WNNLS curve for sepsis."""
import itertools, numpy as np, pandas as pd
PV = "/mnt/user-data/uploads/PBISC/bulkonly_lora5_20261006/analysis/paired_patient_validation_20261014"
P = pd.read_csv(f"{PV}/results/paired_predictions.csv.gz")
ER = pd.read_csv(f"{PV}/results/paired_sample_errors.csv.gz")
rows = []
for coh in P.cohort.unique():
    common = set(ER[(ER.cohort == coh) & (ER.scope == "common_test_paired") & (ER.level == "l1") & ER.method.str.contains("ensemble")].sample_id)
    for K in ["k1000", "k4000"]:
        d = P[(P.cohort == coh) & (P.scope == "all_paired") & (P.method == f"PBISC_{K}") & P.sample_id.isin(common)]
        piv = d.pivot_table(index=["sample_id", "donor_id", "draw"], columns="cell_type", values="pred")
        tru = d.groupby(["sample_id", "cell_type"]).truth.first().unstack()
        for k in range(1, 6):
            for comb in itertools.combinations(range(5), k):
                e = []
                for (s, don), g in piv.groupby(level=[0, 1]):
                    g = g.droplevel([0, 1]); g = g.loc[[x for x in comb if x in g.index]]
                    if len(g) == 0: continue
                    e.append((don, float(np.abs(g.mean(0) - tru.loc[s]).sum())))
                e = pd.DataFrame(e, columns=["donor", "L1"]).groupby("donor").L1.mean()
                rows.append(dict(cohort=coh, K=K, k=k, combo="".join(map(str, comb)), L1=e.mean(), n_patients=len(e)))
R = pd.DataFrame(rows); R.to_csv("/home/claude/work/prop/eq/ensemble_size_curve.csv", index=False)
print(R.groupby(["cohort", "K", "k"]).L1.agg(["mean", "min", "max"]).round(3).to_string())
