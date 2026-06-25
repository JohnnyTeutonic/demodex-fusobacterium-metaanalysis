"""
Figure -- ROC curves for the composite mite-index on the Demodex-vs-control
validation cohort (PRJNA692647), comparing the NAIVE (skin-only background) and
CORRECTED (body-site-matched background) models.

Reads two index files produced by build_index.py:
  data/index_PRJNA692647_naive.tsv   (skin-only background panel)
  data/index_PRJNA692647.tsv         (body-site-matched panel; the corrected run)

Labels come from the sample title (same rule as analyze_validation.py). AUC is
computed via the Mann-Whitney U relationship so it matches the reported numbers
exactly; the curve itself is drawn by the standard threshold sweep.

Output (vector PDF for LaTeX + PNG preview):
  results/fig_roc_validation.pdf
  results/fig_roc_validation.png

Usage:
  python make_roc_figure.py
  python make_roc_figure.py --self-test
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")  # headless: no display needed
import matplotlib.pyplot as plt  # noqa: E402

from analyze_validation import group_of, PROJECT  # reuse labelling

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
RESULTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results")

NAIVE_TSV = os.path.join(DATA_DIR, f"index_{PROJECT}_naive.tsv")
CORRECTED_TSV = os.path.join(DATA_DIR, f"index_{PROJECT}.tsv")
META_TSV = os.path.join(DATA_DIR, f"meta_{PROJECT}.tsv")
VALUE = "mite_index"


def roc_curve(scores: np.ndarray, labels: np.ndarray):
    """Return (fpr, tpr) for the threshold sweep; labels are 1=positive."""
    order = np.argsort(-scores, kind="mergesort")
    y = labels[order]
    P, N = y.sum(), (1 - y).sum()
    tpr = np.concatenate([[0.0], np.cumsum(y) / P])
    fpr = np.concatenate([[0.0], np.cumsum(1 - y) / N])
    return fpr, tpr


def auc_mwu(pos: np.ndarray, neg: np.ndarray) -> float:
    u = stats.mannwhitneyu(pos, neg, alternative="greater").statistic
    return float(u / (len(pos) * len(neg)))


def load_scores(index_tsv: str):
    meta = pd.read_csv(META_TSV, sep="\t", dtype=str)
    idx = pd.read_csv(index_tsv, sep="\t")
    df = meta.merge(idx, on="run_accession", how="inner")
    df["group"] = df["sample_title"].map(group_of)
    df = df[df["group"].isin(["demodex", "control"])]
    scores = df[VALUE].astype(float).values
    labels = (df["group"] == "demodex").astype(int).values
    return scores, labels


def make_figure() -> str:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    models = [
        ("Naive (skin-only background)", NAIVE_TSV, "#9e9e9e", "--"),
        ("Corrected (body-site-matched)", CORRECTED_TSV, "#c62828", "-"),
    ]
    fig, ax = plt.subplots(figsize=(4.4, 4.2))
    for name, tsv, colour, ls in models:
        if not os.path.exists(tsv):
            raise FileNotFoundError(
                f"missing {tsv}\nGenerate it, e.g. for the naive model:\n"
                f"  python fetch_references.py --skin-only --out-dir data/kmers_naive\n"
                f"  python build_index.py --project {PROJECT} "
                f"--kmers-dir data/kmers_naive --out {tsv}")
        scores, labels = load_scores(tsv)
        fpr, tpr = roc_curve(scores, labels)
        a = auc_mwu(scores[labels == 1], scores[labels == 0])
        ax.step(fpr, tpr, where="post", color=colour, ls=ls, lw=2,
                label=f"{name}\nAUC = {a:.2f}")

    ax.plot([0, 1], [0, 1], color="0.7", lw=1, ls=":")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("Composite mite-index: Demodex+ vs control\n(PRJNA692647)",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    ax.set_aspect("equal")
    fig.tight_layout()
    pdf = os.path.join(RESULTS_DIR, "fig_roc_validation.pdf")
    png = os.path.join(RESULTS_DIR, "fig_roc_validation.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=200)
    plt.close(fig)
    print(f"Wrote {pdf}\nWrote {png}")
    return pdf


def _self_test() -> None:
    # perfectly separable -> AUC 1; fully reversed -> AUC 0
    pos = np.array([3.0, 4.0, 5.0])
    neg = np.array([0.0, 1.0, 2.0])
    assert abs(auc_mwu(pos, neg) - 1.0) < 1e-9
    scores = np.concatenate([pos, neg])
    labels = np.array([1, 1, 1, 0, 0, 0])
    fpr, tpr = roc_curve(scores, labels)
    assert fpr[0] == 0.0 and tpr[0] == 0.0
    assert abs(fpr[-1] - 1.0) < 1e-9 and abs(tpr[-1] - 1.0) < 1e-9
    # AUC by trapezoid over the swept curve should also be ~1
    area = float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2.0))
    assert abs(area - 1.0) < 1e-9
    print("self-test OK: roc_curve + auc_mwu pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    make_figure()


if __name__ == "__main__":
    main()
