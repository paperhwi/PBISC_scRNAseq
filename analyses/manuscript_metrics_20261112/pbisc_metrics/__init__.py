"""PBISC manuscript metrics: gene-mean curves, iLISI, composition / within-type decomposition, rare-population statistics.

Self-contained implementations written from the definitions documented in the PBISC project notes (2026-10-20 ... 2026-11-12).
They do not read any F:/ path: every function takes counts / obs tables as arguments.
"""
from . import gene_mean, ilisi, decomposition, rare_stats, synthetic  # noqa: F401

__all__ = ["gene_mean", "ilisi", "decomposition", "rare_stats", "synthetic"]
