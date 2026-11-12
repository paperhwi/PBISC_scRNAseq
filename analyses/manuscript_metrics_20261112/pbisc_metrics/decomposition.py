"""Disease-effect decomposition into a composition part and a within-cell-type part, and noise-corrected comparison with real.

Model. For sample s with n_st cells of type t (p_st = n_st / n_s) and cell-mean CP10k expression mu_stg, the sample pseudobulk rate is
    R_sg = sum_t p_st * mu_stg .
For a case/control contrast define group level p_t = mean_s p_st and mu_tg = mean_s(p_st mu_stg) / p_t, so that the group mean
rate is exactly sum_t p_t mu_tg. With the midpoint (Kitagawa) identity
    p1 mu1 - p0 mu0 = (p1 - p0) (mu1 + mu0)/2  +  (p1 + p0)/2 (mu1 - mu0)
the contrast splits with no interaction residual into
    total_g = composition_g + within_g,
    composition_g = sum_t dp_t * mean(mu_t)_g      (shift of cell-type proportions),
    within_g      = sum_t mean(p_t) * dmu_tg       (change of expression inside each type).

Noise-corrected slope (patient split-half). The real effect is measured twice on disjoint halves of the patients (A and B, split
inside each group; sorted unit ids, fixed seed, `n_splits` splits). With g the generated effect vector,
    slope = cov(g, r_B) / cov(r_A, r_B)            (r_B instruments the noisy r_A; the ratio estimates g regressed on the true effect),
    r_dis = corr(g, r_B) / sqrt(corr(r_A, r_B))    (disattenuated correlation).
The same is applied to the total, composition, within and per-type (`type:<name>`) effect vectors. `split_cv` is the SD of the
slope over splits divided by |mean slope|; the composition part is the least certain.
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._utils import align_genes, cp10k, to_csr


@dataclass
class TypeAgg:
    """Sample x cell-type aggregate. S[s, t, g] = sum of CP10k over the cells of (sample, type); n[s, t] = number of cells."""
    S: np.ndarray
    n: np.ndarray
    samples: np.ndarray
    types: np.ndarray
    genes: np.ndarray
    meta: pd.DataFrame  # index = samples; columns include `group` and `patient`


def aggregate(X, genes, sample, celltype, meta, types=None):
    """Build a TypeAgg from a cells x genes count matrix and per-cell sample / cell-type labels.

    meta: DataFrame indexed by sample id with at least `group` (e.g. 'case' / 'control') and `patient`.
    """
    sample, celltype = np.asarray(sample).astype(str), np.asarray(celltype).astype(str)
    samples = np.array(sorted(set(sample)))
    types = np.array(sorted(set(celltype)) if types is None else list(types))
    si = pd.Series(np.arange(len(samples)), index=samples)[sample].to_numpy()
    ti = pd.Series(np.arange(len(types)), index=types).reindex(celltype).to_numpy()
    keep = ~np.isnan(ti)
    ti = ti[keep].astype(int)
    ind = sp.csr_matrix((np.ones(keep.sum()), (si[keep] * len(types) + ti, np.flatnonzero(keep))), shape=(len(samples) * len(types), len(sample)))
    S = (ind @ cp10k(X)).toarray().reshape(len(samples), len(types), -1)
    n = np.asarray(ind.sum(1)).reshape(len(samples), len(types)).astype(float)
    return TypeAgg(S=S, n=n, samples=samples, types=types, genes=np.asarray(genes).astype(str), meta=meta.loc[samples].copy())


def restrict_genes(agg, genes):
    """Same aggregate on a subset / reordering of genes."""
    pos = pd.Series(np.arange(len(agg.genes)), index=agg.genes)
    return TypeAgg(agg.S[:, :, pos[genes].to_numpy()], agg.n, agg.samples, agg.types, np.asarray(genes), agg.meta)


def group_terms(agg, mask):
    """Group-level (p_t, mu_tg) of the samples in `mask` (bool over agg.samples)."""
    n_s = agg.n[mask].sum(1, keepdims=True)
    p_st = agg.n[mask] / np.maximum(n_s, 1)
    p = p_st.mean(0)
    pmu = (agg.S[mask] / np.maximum(n_s, 1)[:, :, None]).mean(0)  # mean_s p_st mu_stg
    mu = pmu / np.maximum(p, 1e-12)[:, None]
    return p, mu


def decompose(agg, case, ctrl):
    """Total / composition / within effect vectors (case - control) plus per-type effects and proportion shifts.

    case, ctrl: boolean masks over agg.samples. Returns dict with
      'effects'  : DataFrame genes x [total, composition, within, type:<t> ...]  (type:<t> = unweighted mean-expression change inside t),
      'dprop_pp' : Series of proportion change per type in percentage points.
    """
    p1, m1 = group_terms(agg, case)
    p0, m0 = group_terms(agg, ctrl)
    comp = ((p1 - p0)[:, None] * (m1 + m0) / 2).sum(0)
    within = (((p1 + p0) / 2)[:, None] * (m1 - m0)).sum(0)
    total = (p1[:, None] * m1).sum(0) - (p0[:, None] * m0).sum(0)
    eff = pd.DataFrame({"total": total, "composition": comp, "within": within}, index=agg.genes)
    for t, name in enumerate(agg.types):
        eff[f"type:{name}"] = m1[t] - m0[t]
    return {"effects": eff, "dprop_pp": pd.Series((p1 - p0) * 100, index=agg.types)}


def _split_masks(agg, rng):
    """Random patient-level halves, stratified by group. Units are sorted patient ids, so a fixed seed gives a fixed split."""
    A = np.zeros(len(agg.samples), bool)
    for grp in sorted(agg.meta.group.unique()):
        pats = np.array(sorted(agg.meta.patient[agg.meta.group == grp].unique()))
        half = set(rng.permutation(pats)[: len(pats) // 2])
        A |= ((agg.meta.group == grp) & agg.meta.patient.isin(half)).to_numpy()
    return A, ~A


def split_half_slopes(real, gen_effects, case_group, ctrl_group, n_splits=200, seed=0, min_cells=0):
    """Noise-corrected slope and disattenuated r of generated effects against real effects (see module docstring).

    real: TypeAgg of real cells; gen_effects: DataFrame genes x effect names from `decompose` on the generated cells
    (the same gene order as `real.genes`; use `restrict_genes` first). Returns DataFrame indexed by effect with
    slope, slope_sd, split_cv, r_dis, r_dis_sd, reliability (corr of the two real halves), n_splits_used.
    """
    rng = np.random.default_rng(seed)
    is_case, is_ctrl = (real.meta.group == case_group).to_numpy(), (real.meta.group == ctrl_group).to_numpy()
    names = list(gen_effects.columns)
    acc = {k: {"slope": [], "r": [], "rel": []} for k in names}
    for _ in range(n_splits):
        A, B = _split_masks(real, rng)
        try:
            eA = decompose(real, A & is_case, A & is_ctrl)["effects"]
            eB = decompose(real, B & is_case, B & is_ctrl)["effects"]
        except Exception:
            continue
        for k in names:
            a, b, g = eA[k].to_numpy(), eB[k].to_numpy(), gen_effects[k].to_numpy()
            ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(g)
            if ok.sum() < 10:
                continue
            a, b, g = a[ok], b[ok], g[ok]
            cab = np.cov(a, b)[0, 1]
            rel = np.corrcoef(a, b)[0, 1]
            if cab <= 0 or rel <= 0:
                continue
            acc[k]["slope"].append(np.cov(g, b)[0, 1] / cab)
            acc[k]["r"].append(np.corrcoef(g, b)[0, 1] / np.sqrt(rel))
            acc[k]["rel"].append(rel)
    rows = []
    for k in names:
        s, r, rel = map(np.asarray, (acc[k]["slope"], acc[k]["r"], acc[k]["rel"]))
        m = float(s.mean()) if len(s) else np.nan
        rows.append({"effect": k, "slope": m, "slope_sd": float(s.std(ddof=1)) if len(s) > 1 else np.nan,
                     "split_cv": float(s.std(ddof=1) / abs(m)) if len(s) > 1 and m else np.nan,
                     "r_dis": float(r.mean()) if len(r) else np.nan, "r_dis_sd": float(r.std(ddof=1)) if len(r) > 1 else np.nan,
                     "reliability": float(rel.mean()) if len(rel) else np.nan, "n_splits_used": int(len(s))})
    return pd.DataFrame(rows).set_index("effect")


def decomposition_report(real, gen, case_group, ctrl_group, **kw):
    """Convenience wrapper: align genes, decompose generated and real (all samples), and return
    (slopes DataFrame, real effects, generated effects, dprop table with real and generated proportion shifts)."""
    ia, ib, names = align_genes(real.genes, gen.genes)
    real, gen = restrict_genes(real, names), restrict_genes(gen, names)
    mk = lambda a: ((a.meta.group == case_group).to_numpy(), (a.meta.group == ctrl_group).to_numpy())
    dr, dg = decompose(real, *mk(real)), decompose(gen, *mk(gen))
    slopes = split_half_slopes(real, dg["effects"], case_group, ctrl_group, **kw)
    dprop = pd.DataFrame({"real_pp": dr["dprop_pp"], "generated_pp": dg["dprop_pp"]})
    return slopes, dr["effects"], dg["effects"], dprop
