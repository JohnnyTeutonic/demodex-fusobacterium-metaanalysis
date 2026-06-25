"""
E3 -- ivermectin treatment arm (PRJDB18292).

Design: 6 subjects x {lesional, non-lesional} x {pre, post} topical ivermectin
= 24 paired skin-swab 16S samples. Hypothesis from the clinical argument: a pure
antiparasitic should reduce the mite burden, so the molecular mite-burden index
(and especially Demodex-associated taxa) should fall pre -> post.

We pair within subject x site and run a Wilcoxon signed-rank test (n=6 pairs per
site; also a sign test, since n is tiny) on the mite_index and on each
Demodex-associated taxon. Honest reporting: with n=6 this is exploratory.

Usage:
  python analyze_treatment.py
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results")
PROJECT = "PRJDB18292"

DEMODEX_TAXA = ["Snodgrassella_alvi", "Bacillus_oleronius",
                "Cutibacterium_kroppenstedtii", "Bartonella_quintana"]


def load() -> pd.DataFrame:
    meta = pd.read_csv(os.path.join(DATA_DIR, f"meta_{PROJECT}.tsv"), sep="\t",
                       dtype=str)
    idx = pd.read_csv(os.path.join(DATA_DIR, f"index_{PROJECT}.tsv"), sep="\t")
    df = meta.merge(idx, on="run_accession", how="inner")
    df["host_subject_id"] = df["host_subject_id"].astype(str)
    return df


def paired_matrix(df: pd.DataFrame, site: str, value: str):
    sub = df[df["lesional"] == site]
    piv = sub.pivot_table(index="host_subject_id", columns="timepoint",
                          values=value, aggfunc="mean")
    piv = piv.dropna(subset=["pre", "post"])
    return piv


def test_pair(piv: pd.DataFrame, label: str) -> dict:
    pre, post = piv["pre"].values, piv["post"].values
    delta = post - pre
    n = len(delta)
    n_down = int(np.sum(delta < 0))
    n_up = int(np.sum(delta > 0))
    res = {"stratum": label, "n_pairs": n, "n_down": n_down, "n_up": n_up,
           "median_pre": float(np.median(pre)), "median_post": float(np.median(post)),
           "median_delta": float(np.median(delta))}
    try:
        w = stats.wilcoxon(pre, post, zero_method="wilcox", alternative="two-sided")
        res["wilcoxon_p"] = float(w.pvalue)
    except Exception:
        res["wilcoxon_p"] = float("nan")
    # exact sign test (binomial) on direction
    k = min(n_down, n_up)
    res["sign_test_p"] = float(stats.binomtest(n_down, n_down + n_up, 0.5).pvalue) \
        if (n_down + n_up) > 0 else float("nan")
    return res


def main() -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    df = load()
    df.to_csv(os.path.join(RESULTS_DIR, f"merged_{PROJECT}.tsv"), sep="\t", index=False)

    rows = []
    for site in ["lesional", "non-lesional"]:
        for value in ["mite_index"] + [t + "_frac" for t in DEMODEX_TAXA]:
            piv = paired_matrix(df, site, value)
            if piv.empty:
                continue
            r = test_pair(piv, f"{site} / {value}")
            rows.append(r)

    out = pd.DataFrame(rows)
    pd.set_option("display.width", 160)
    pd.set_option("display.max_columns", 20)
    print(out.to_string(index=False))
    out.to_csv(os.path.join(RESULTS_DIR, f"treatment_tests_{PROJECT}.tsv"),
               sep="\t", index=False)

    # also dump the per-subject paired mite_index tables for transparency
    print("\nPer-subject mite_index (pre -> post):")
    for site in ["lesional", "non-lesional"]:
        piv = paired_matrix(df, site, "mite_index")
        print(f"\n[{site}]")
        show = piv.copy()
        show["delta"] = show["post"] - show["pre"]
        print(show.round(5).to_string())
    print(f"\nWrote results to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
