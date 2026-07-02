"""Per-cohort AUC dot-plot for the three-cohort meta-analysis.
Shows which genera reproduce (all three cohorts on the same side of AUC=0.5)
versus the field's single-cohort leads that do not. Numbers are taken directly
from rosacea_meta3.md (discovery / rep2020 / shotgun2022)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# genus: (discovery, rep2020, shotgun2022) AUC for demodex>control
DATA = [
    ("Fusobacterium",  (0.332, 0.305, 0.360), "reproducible (down)"),
    ("Corynebacterium",(0.685, 0.555, 0.520), "reproducible (up)"),
    ("Streptococcus",  (0.664, 0.617, 0.376), "not reproducible"),
    ("Rothia",         (0.729, 0.590, 0.347), "not reproducible"),
    ("Neisseria",      (0.689, 0.321, 0.222), "not reproducible"),
]
COHORTS = ["discovery (16S V4)", "rep2020 (16S V3-V4)", "shotgun (WGS)"]
MARKERS = ["o", "s", "^"]
COLORS = ["#1b5e20", "#1565c0", "#c62828"]

fig, ax = plt.subplots(figsize=(7.2, 3.6))
ys = list(range(len(DATA)))[::-1]
for y, (genus, aucs, _) in zip(ys, DATA):
    for auc, mk, col in zip(aucs, MARKERS, COLORS):
        ax.scatter(auc, y, marker=mk, s=70, color=col, zorder=3,
                   edgecolor="white", linewidth=0.6)
ax.axvline(0.5, color="0.4", ls="--", lw=1, zorder=1)
ax.set_yticks(ys)
ax.set_yticklabels([g for g, _, _ in DATA], fontstyle="italic")
ax.set_xlim(0.15, 0.85)
ax.set_xlabel("per-cohort ROC-AUC (Demodex > control)")
ax.text(0.34, len(DATA) - 0.35, "depleted", ha="center", color="0.35", fontsize=9)
ax.text(0.66, len(DATA) - 0.35, "enriched", ha="center", color="0.35", fontsize=9)
handles = [plt.Line2D([0], [0], marker=mk, color="w", markerfacecolor=col,
                      markeredgecolor="white", markersize=9, label=lab)
           for mk, col, lab in zip(MARKERS, COLORS, COHORTS)]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
out = "results/fig_meta_auc.png"
import os
os.makedirs("results", exist_ok=True)
fig.savefig(out, dpi=200)
print("wrote", os.path.abspath(out))
