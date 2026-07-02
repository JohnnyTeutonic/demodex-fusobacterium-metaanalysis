"""
Figure 2 -- closed-reference diagnostic k-mers cross-map to off-target taxa.

Parses the k-mer/ASV bridge report (E6b_kmer_bridge.md) and plots, for each
closed-reference panel, the Demodex+ vs control AUC of every ASV the panel
flagged. Each point is one ASV (size proportional to read abundance), coloured by
whether SILVA assigns it to the panel's nominal target genus (on-target) or not
(off-target). A vertical line marks chance (AUC = 0.5). The single significant
ASV (the Neisseria lead, ASV43) is highlighted.

The point: the "signal" each panel detects scatters around chance and is
overwhelmingly off-target -- the S. alvi panel flags no Snodgrassella at all --
demonstrating that closed-reference marker counting is not species-specific on
short 16S.

Output:
  results/fig_crossmap.pdf
  results/fig_crossmap.png

Usage:
  python make_crossmap_figure.py
  python make_crossmap_figure.py --self-test
"""

from __future__ import annotations

import argparse
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(__file__)
BRIDGE_MD = os.path.join(HERE, os.pardir, os.pardir, "data",
                         "dada2_PRJNA692647", "E6b_kmer_bridge.md")
RESULTS_DIR = os.path.join(HERE, os.pardir, os.pardir, "results")

# Panel -> (display label, nominal target genus)
PANELS = {
    "Snodgrassella_alvi": ("Snodgrassella alvi\npanel", "Snodgrassella"),
    "Bacillus_oleronius": ("Bacillus oleronius\npanel", "Bacillus"),
    "Bartonella_quintana": ("Bartonella quintana\npanel", "Bartonella"),
}

HEADER_RE = re.compile(r"^##\s+(\S+)\s+\(panel")
ROW_RE = re.compile(
    r"\*\*(ASV\d+)\*\*:.*?total reads=(\d+),.*?Genus=([^;]+);.*?AUC=([0-9.]+)"
)


def parse_bridge(md_path: str):
    """Return list of dicts: {panel, asv, reads, genus, auc}."""
    rows = []
    panel = None
    with open(md_path, encoding="utf-8") as fh:
        for line in fh:
            h = HEADER_RE.match(line)
            if h:
                panel = h.group(1)
                continue
            m = ROW_RE.search(line)
            if m and panel in PANELS:
                rows.append({
                    "panel": panel,
                    "asv": m.group(1),
                    "reads": int(m.group(2)),
                    "genus": m.group(3).strip(),
                    "auc": float(m.group(4)),
                })
    return rows


def make_figure(rows) -> str:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    panel_keys = list(PANELS.keys())
    y_of = {p: i for i, p in enumerate(panel_keys)}

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    rng = np.random.default_rng(0)

    def msize(reads: int) -> float:
        return 18.0 + 34.0 * np.log10(reads + 1)

    on_x, on_y, on_s = [], [], []
    off_x, off_y, off_s = [], [], []
    for r in rows:
        target = PANELS[r["panel"]][1]
        y = y_of[r["panel"]] + rng.uniform(-0.18, 0.18)
        if r["genus"].lower() == target.lower():
            on_x.append(r["auc"]); on_y.append(y); on_s.append(msize(r["reads"]))
        else:
            off_x.append(r["auc"]); off_y.append(y); off_s.append(msize(r["reads"]))

    ax.scatter(off_x, off_y, s=off_s, c="#9e9e9e", alpha=0.65,
               edgecolors="white", linewidths=0.5, label="off-target genus", zorder=2)
    ax.scatter(on_x, on_y, s=on_s, c="#2e7d32", alpha=0.85,
               edgecolors="white", linewidths=0.5, label="on-target genus", zorder=3)

    # Highlight the Neisseria lead (ASV43, the one significant separation).
    lead = next((r for r in rows if r["asv"] == "ASV43"), None)
    if lead is not None:
        ly = y_of[lead["panel"]]
        ax.scatter([lead["auc"]], [ly], s=msize(lead["reads"]) + 30,
                   facecolors="none", edgecolors="#c62828", linewidths=2.0, zorder=4)
        ax.annotate("Neisseria ASV43\nAUC 0.70 (p = 0.032)",
                    xy=(lead["auc"], ly), xytext=(lead["auc"] - 0.015, ly + 0.5),
                    fontsize=8, color="#c62828", ha="right",
                    arrowprops=dict(arrowstyle="->", color="#c62828", lw=1.2))

    ax.axvline(0.5, color="0.4", ls=":", lw=1.2, zorder=1)
    ax.text(0.5, len(panel_keys) - 0.42, "chance", rotation=90,
            va="top", ha="right", fontsize=8, color="0.4")

    ax.set_yticks(range(len(panel_keys)))
    ax.set_yticklabels([PANELS[p][0] for p in panel_keys], fontsize=9, style="italic")
    ax.set_ylim(-0.6, len(panel_keys) - 0.4)
    ax.set_xlim(0.24, 0.82)
    ax.set_xlabel("Demodex+ vs control AUC of flagged ASV")
    ax.set_title("Diagnostic k-mers cross-map to off-target taxa at ASV resolution",
                 fontsize=10)

    handles = [
        Line2D([0], [0], marker="o", ls="", mfc="#9e9e9e", mec="white",
               ms=8, label="off-target genus"),
        Line2D([0], [0], marker="o", ls="", mfc="#2e7d32", mec="white",
               ms=8, label="on-target genus"),
        Line2D([0], [0], marker="o", ls="", mfc="none", mec="#c62828",
               mew=2, ms=10, label="significant lead"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()

    pdf = os.path.join(RESULTS_DIR, "fig_crossmap.pdf")
    png = os.path.join(RESULTS_DIR, "fig_crossmap.png")
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    print(f"Wrote {pdf}\nWrote {png}")
    return pdf


def _self_test() -> None:
    sample = (
        "## Snodgrassella_alvi (panel: 891 diagnostic 31-mers)\n"
        "- **ASV43**: 36 diagnostic k-mers, total reads=4367, "
        "SILVA=[Family=Neisseriaceae; Genus=Neisseria; Species=NA], "
        "demodex-vs-control AUC=0.697 (p=0.032).\n"
        "## Bacillus_oleronius (panel: 2558 diagnostic 31-mers)\n"
        "- **ASV40**: 3 diagnostic k-mers, total reads=4712, "
        "SILVA=[Family=Bacillaceae; Genus=Bacillus; Species=NA], "
        "demodex-vs-control AUC=0.643 (p=0.092).\n"
    )
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False,
                                     encoding="utf-8") as fh:
        fh.write(sample)
        path = fh.name
    rows = parse_bridge(path)
    os.unlink(path)
    assert len(rows) == 2, rows
    a43 = next(r for r in rows if r["asv"] == "ASV43")
    assert a43["genus"] == "Neisseria" and a43["reads"] == 4367
    assert abs(a43["auc"] - 0.697) < 1e-9
    a40 = next(r for r in rows if r["asv"] == "ASV40")
    assert a40["genus"] == "Bacillus" and a40["panel"] == "Bacillus_oleronius"
    print("self-test OK: bridge parser pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    rows = parse_bridge(BRIDGE_MD)
    print(f"parsed {len(rows)} flagged ASVs across panels")
    make_figure(rows)


if __name__ == "__main__":
    main()
