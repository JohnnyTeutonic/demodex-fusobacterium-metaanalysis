"""
E4 -- validation arm (PRJNA692647): does the molecular mite-burden index carry a
real Demodex signal?

This ocular-surface 16S study has two groups encoded in the sample TITLE:
"demodex infection N" vs "healthy control N" (14 vs 17). If the index means
anything, Demodex-infection samples should score higher.

We run Mann-Whitney U (two-sided + one-sided demodex>control) on the mite_index
and each Demodex-associated taxon, plus a simple ROC-AUC for the composite index.

Usage:
  python analyze_validation.py
"""

from __future__ import annotations

import os
import re
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results")
PROJECT = "PRJNA692647"

DEMODEX_TAXA = ["Snodgrassella_alvi", "Bacillus_oleronius",
                "Cutibacterium_kroppenstedtii", "Bartonella_quintana"]


def group_of(title: str) -> str:
    t = title.lower()
    if "demodex" in t:
        return "demodex"
    if "healthy" in t or "control" in t:
        return "control"
    return "unknown"


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """ROC-AUC via the Mann-Whitney U relationship (P(pos > neg))."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    u = stats.mannwhitneyu(pos, neg, alternative="greater").statistic
    return float(u / (len(pos) * len(neg)))


def main() -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    meta = pd.read_csv(os.path.join(DATA_DIR, f"meta_{PROJECT}.tsv"), sep="\t", dtype=str)
    idx = pd.read_csv(os.path.join(DATA_DIR, f"index_{PROJECT}.tsv"), sep="\t")
    df = meta.merge(idx, on="run_accession", how="inner")
    df["group"] = df["sample_title"].map(group_of)
    df = df[df["group"].isin(["demodex", "control"])]
    df.to_csv(os.path.join(RESULTS_DIR, f"merged_{PROJECT}.tsv"), sep="\t", index=False)

    n_dem = int((df["group"] == "demodex").sum())
    n_ctl = int((df["group"] == "control").sum())
    print(f"{PROJECT}: demodex={n_dem}  control={n_ctl}\n")

    rows = []
    for value in ["mite_index"] + [t + "_frac" for t in DEMODEX_TAXA]:
        pos = df.loc[df["group"] == "demodex", value].astype(float).values
        neg = df.loc[df["group"] == "control", value].astype(float).values
        u2 = stats.mannwhitneyu(pos, neg, alternative="two-sided")
        u1 = stats.mannwhitneyu(pos, neg, alternative="greater")
        rows.append({
            "metric": value,
            "median_demodex": float(np.median(pos)),
            "median_control": float(np.median(neg)),
            "log2_fold": float(np.log2((np.median(pos) + 1e-6) / (np.median(neg) + 1e-6))),
            "auc_demodex_gt_control": auc(pos, neg),
            "mwu_p_two_sided": float(u2.pvalue),
            "mwu_p_demodex_gt": float(u1.pvalue),
        })

    out = pd.DataFrame(rows)
    pd.set_option("display.width", 170)
    pd.set_option("display.max_columns", 20)
    print(out.to_string(index=False))
    out.to_csv(os.path.join(RESULTS_DIR, f"validation_tests_{PROJECT}.tsv"),
               sep="\t", index=False)
    print(f"\nWrote results to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
