"""Recreate pooled rare-subtype shares with the new LoRA kNN labels."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
OLD = ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/ilisi_floor_20261021"
sys.path.insert(0, str(ROOT / "bulkonly_lora5_20261006/analysis/integration_report_20261012/disease_deg_rare_20261027/code"))
from rare_lib import SEPSIS, PAR_S, COVID, PAR_C

rows = []
for coh, definitions, parents, case, control in (
    ("sepsis", SEPSIS, PAR_S, "adult_sepsis", "adult_HC"),
    ("covid", COVID, PAR_C, "covid_positive", "covid_negative"),
):
    labels = pd.read_csv(ROOT / "NewLora/analysis_20261102" / coh / "ceiling_cell_labels.csv.gz", keep_default_na=False)
    old = pd.read_csv(OLD / "results/big" / f"{coh}_rare_parent_share.csv")
    old = old[old.source.isin(("real B (author)", "real B (transferred)", "prev K4000 (transferred)", "floor (transferred)"))]
    rows.append(old)
    major, fine, group = (labels[col].astype(str).to_numpy() for col in ("major_tr", "fine_tr", "group"))
    for subset, parent_fn in definitions.items():
        if subset not in set(old.subset):
            continue
        parent = parent_fn(major, fine)
        target = parent & (fine == subset)
        for grp in ("all", case, control):
            gm = np.ones(len(group), dtype=bool) if grp == "all" else group == grp
            n_target, n_parent = int((target & gm).sum()), int((parent & gm).sum())
            rows.append(pd.DataFrame([dict(cohort=coh, subset=subset, parent=parents[subset], source="new LoRA ceiling (transferred)", group=grp,
                                           n_target=n_target, n_parent=n_parent, pct_of_parent=100*n_target/max(n_parent,1))]))
    print(coh, "new labels", len(labels), flush=True)

result = pd.concat(rows, ignore_index=True)
OUT.mkdir(parents=True, exist_ok=True)
result.to_csv(OUT / "rare_parent_share_with_new_lora.csv", index=False)
