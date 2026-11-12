"""Small synthetic cohorts with a known truth, for tests and the demo notebook.

A `world` fixes gene baselines, cell-type profiles, rare-state markers and the disease effect. `sample_cohort` draws patients from it, so
'real A', 'real B' (independent patients, same world) and a 'generated' cohort (same world but with attenuated effects, a depleted rare
state and a blunted low-expression tail) share one truth and can be compared with the metrics.
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp

TYPES = ("T", "B", "NK", "Mono")
# fine label -> parent
FINE = {"T_naive": "T", "T_cyc": "T", "B_naive": "B", "NK_dim": "NK", "NK_bright": "NK", "Mono_classical": "Mono", "Mono_inflam": "Mono"}
FINE_P = {"T_naive": 0.97, "T_cyc": 0.03, "B_naive": 1.0, "NK_dim": 0.93, "NK_bright": 0.07, "Mono_classical": 0.95, "Mono_inflam": 0.05}


def make_world(n_genes=600, seed=0):
    rng = np.random.default_rng(seed)
    genes = np.array([f"G{i:04d}" for i in range(n_genes)])
    base = rng.lognormal(-1.2, 1.4, n_genes)
    prof = {t: base * np.exp(rng.normal(0, 0.9, n_genes)) for t in TYPES}
    fine_prof = {}
    for f, par in FINE.items():
        mult = np.ones(n_genes)
        if FINE_P[f] < 0.1:
            mult[rng.choice(n_genes, 25, replace=False)] = 6.0
        fine_prof[f] = prof[par] * mult
    eff = {t: np.where(rng.random(n_genes) < 0.2, rng.normal(0, 0.6, n_genes), 0.0) for t in TYPES}
    p0 = np.array([0.45, 0.15, 0.15, 0.25])
    p1 = np.array([0.30, 0.10, 0.08, 0.52])
    return {"genes": genes, "fine_prof": fine_prof, "eff": eff, "p0": p0, "p1": p1, "n_genes": n_genes}


def sample_cohort(world, n_per_group=12, cells_per_sample=400, seed=1, comp_scale=1.0, within_scale=1.0, drop_fine=(), depth_mean=1500.0,
                  lowexpr_floor=0.0, prefix="p"):
    """Draw a cohort. comp_scale / within_scale shrink the disease effect (1 = real); drop_fine removes fine states; lowexpr_floor zeroes
    expected rates below that value (a crude 'cliff' for low-expression genes)."""
    rng = np.random.default_rng(seed)
    G = world["n_genes"]
    rows, blocks = [], []
    types = list(TYPES)
    for grp, pg in (("control", world["p0"]), ("case", world["p0"] + comp_scale * (world["p1"] - world["p0"]))):
        for i in range(n_per_group):
            pat = f"{prefix}_{grp}_{i:02d}"
            p = rng.dirichlet(pg * 60)
            tcount = rng.multinomial(cells_per_sample, p)
            for t, nt in zip(types, tcount):
                if nt == 0:
                    continue
                fines = [f for f, par in FINE.items() if par == t and f not in drop_fine]
                fp = np.array([FINE_P[f] for f in fines]); fp = fp / fp.sum()
                fcount = rng.multinomial(nt, fp)
                for f, nf in zip(fines, fcount):
                    if nf == 0:
                        continue
                    rate = world["fine_prof"][f] * np.exp(within_scale * (world["eff"][t] if grp == "case" else 0.0))
                    rate = rate * np.exp(rng.normal(0, 0.15, G))  # patient-level noise
                    d = rng.lognormal(np.log(depth_mean), 0.35, nf)
                    mu = d[:, None] * (rate / rate.sum())[None, :]
                    if lowexpr_floor > 0:
                        mu = np.where(mu < lowexpr_floor, 0.0, mu)
                    blocks.append(sp.csr_matrix(rng.poisson(mu).astype(np.float64)))
                    rows += [(pat, pat, grp, t, f)] * nf
    obs = pd.DataFrame(rows, columns=["sample", "patient", "group", "parent", "fine"])
    return {"X": sp.vstack(blocks).tocsr(), "genes": world["genes"], "obs": obs.reset_index(drop=True)}
