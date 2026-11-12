"""Rare-population statistics within parent lineages.

Definitions (PBISC notes 2026-10-12 / 2026-10-21 / 2026-10-30, f_rare_parent.py):
* share = cells of the fine label / cells of its parent lineage (so a MKI67+ T state is measured among T cells, not among all cells);
* ratio = generated share / reference (real A) share; verdict: below (< 0.5x), recovered (0.5-2x), over (> 2x),
  not_assessable when the reference share is 0 (or the real control, if given, is 0 as well: the mapping itself loses the state);
* a patient / sample bootstrap 95% CI of the ratio (units resampled with replacement in each source separately);
* low_n flag when fewer than `min_cells` generated cells carry the label; per-sample detection (samples with >= 1 cell);
* binomial power: probability that a state of a given share is not seen in N cells, (1 - share)^N.
"""
import numpy as np
import pandas as pd


def power_no_cell(share, n_cells):
    """P(no cell of a state with this share in a sample of n_cells)."""
    return float((1.0 - np.asarray(share, dtype=float)) ** np.asarray(n_cells))


def _shares(obs, fine, parent, unit, ids):
    """Per-unit counts: DataFrame units x fine labels, and parent totals units x parent labels."""
    f = pd.crosstab(obs[unit], obs[fine]).reindex(ids, fill_value=0)
    p = pd.crosstab(obs[unit], obs[parent]).reindex(ids, fill_value=0)
    return f, p


def _verdict(ratio, lo=0.5, hi=2.0):
    if not np.isfinite(ratio):
        return "not_assessable"
    return "below" if ratio < lo else ("over" if ratio > hi else "recovered")


def rare_table(obs_ref, obs_gen, fine="fine", parent="parent", unit="patient", sample="sample", obs_ctrl=None,
               rare_threshold=0.05, min_cells=20, n_boot=1000, seed=0):
    """Share, ratio, verdict and bootstrap CI of every fine label.

    obs_*: DataFrames with columns fine, parent, unit (bootstrap unit, normally patient) and sample. obs_ref = real A, obs_gen = generated,
    obs_ctrl = real B (optional ceiling: the same ratio for a second real sample).
    A label is 'is_rare' when its reference share inside its parent is below `rare_threshold`.
    """
    rng = np.random.default_rng(seed)
    labels = sorted(set(obs_ref[fine]) | set(obs_gen[fine]) | (set(obs_ctrl[fine]) if obs_ctrl is not None else set()))
    parent_of = pd.concat([o[[fine, parent]] for o in (obs_ref, obs_gen, obs_ctrl) if o is not None]).drop_duplicates(fine).set_index(fine)[parent]
    srcs = {"ref": obs_ref, "gen": obs_gen, **({"ctrl": obs_ctrl} if obs_ctrl is not None else {})}

    def counts(o):
        n_f = o[fine].value_counts().reindex(labels, fill_value=0)
        n_p = o[parent].value_counts()
        return n_f, n_f.index.map(lambda l: n_p.get(parent_of[l], 0)).to_numpy()

    tab = pd.DataFrame(index=pd.Index(labels, name=fine))
    tab["parent"] = parent_of.reindex(labels)
    for k, o in srcs.items():
        n_f, n_p = counts(o)
        tab[f"n_{k}"] = n_f.to_numpy()
        with np.errstate(divide="ignore", invalid="ignore"):
            tab[f"share_{k}"] = np.where(n_p > 0, n_f.to_numpy() / n_p, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        tab["ratio"] = np.where(tab.share_ref > 0, tab.share_gen / tab.share_ref, np.nan)
        if obs_ctrl is not None:
            tab["ratio_ctrl"] = np.where(tab.share_ref > 0, tab.share_ctrl / tab.share_ref, np.nan)
    tab["is_rare"] = tab.share_ref < rare_threshold
    tab["low_n"] = tab.n_gen < min_cells

    # patient bootstrap of the ratio
    if n_boot:
        boots = {}
        for k in ("ref", "gen"):
            o = srcs[k]
            ids = np.array(sorted(o[unit].astype(str).unique()))
            f, p = _shares(o.assign(**{unit: o[unit].astype(str)}), fine, parent, unit, ids)
            f = f.reindex(columns=labels, fill_value=0)
            pcol = parent_of.reindex(labels).to_numpy()
            P = np.column_stack([p[pp].to_numpy() if pp in p else np.zeros(len(ids)) for pp in pcol])
            boots[k] = (f.to_numpy(), P)
        B = np.empty((n_boot, len(labels)))
        for b in range(n_boot):
            sh = {}
            for k, (F, P) in boots.items():
                pick = rng.integers(0, len(F), len(F))
                num, den = F[pick].sum(0), P[pick].sum(0)
                sh[k] = np.where(den > 0, num / np.maximum(den, 1), np.nan)
            B[b] = np.where(sh["ref"] > 0, sh["gen"] / sh["ref"], np.nan)
        tab["ratio_ci_lo"], tab["ratio_ci_hi"] = np.nanpercentile(B, 2.5, axis=0), np.nanpercentile(B, 97.5, axis=0)

    # per-sample detection of the generated cells
    cs = pd.crosstab(obs_gen[sample], obs_gen[fine]).reindex(columns=labels, fill_value=0)
    tab["gen_samples_detected"] = (cs > 0).sum(axis=0).reindex(labels).to_numpy()
    tab["gen_samples_total"] = len(cs)
    n_med = float(obs_gen.groupby(sample).size().median())
    tab["p_miss_per_sample"] = [power_no_cell(x, n_med) for x in tab.n_ref / max(len(obs_ref), 1)]  # P(state absent from a median-sized sample)

    tab["verdict"] = [_verdict(r) for r in tab.ratio]
    if obs_ctrl is not None:
        tab.loc[(tab.share_ctrl == 0) & (tab.share_ref > 0), "verdict"] = "not_assessable"  # the mapping loses it for real cells too
    tab.loc[tab.share_ref == 0, "verdict"] = "not_assessable"
    return tab.reset_index()


def verdict_summary(tab, only_rare=True):
    """Number of fine labels per verdict (optionally only the rare ones)."""
    t = tab[tab.is_rare] if only_rare else tab
    return t.verdict.value_counts().reindex(["recovered", "below", "over", "not_assessable"], fill_value=0)


def per_patient_counts(obs, fine="fine", unit="patient"):
    """Patients x fine-label count matrix (for per-patient rare-population plots)."""
    return pd.crosstab(obs[unit], obs[fine])
