"""
Cross-cohort genus-level concordance for ocular Demodex disease vs control.

For each cohort (DADA2 out dir + sample sheet), compute the genus-level relative
abundance and a Mann-Whitney AUC for demodex>control. Then find genera that
separate disease from control CONCORDANTLY (same direction) across cohorts, with
a cross-cohort prevalence filter to suppress study-specific contaminants (which
are typically present in only one cohort).

A reproducible signal = concordant direction + non-trivial effect in >=2 cohorts,
present at reasonable prevalence in all of them. This is the positive question:
is there ANY microbial signature of ocular Demodex disease that replicates?

Pure stdlib. Usage:
  python cross_cohort_concordance.py \
    --cohort /home/jonat/rosacea_dada2/out:/home/jonat/rosacea_dada2/sample_sheet.tsv:discovery \
    --cohort /home/jonat/rosacea_657256/out:/home/jonat/rosacea_657256/sample_sheet.tsv:rep2020 \
    --min-prev 0.5 --out /home/jonat/rosacea_concordance.md
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Dict, List, Tuple


def _read_tsv(path: str) -> List[Dict[str, str]]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_labels(sample_sheet: str) -> Dict[str, str]:
    return {r["run_accession"]: r["group"] for r in _read_tsv(sample_sheet)
            if r.get("group") in ("demodex", "control")}


def load_cohort(out_dir: str, sample_sheet: str):
    with open(os.path.join(out_dir, "asv_table.tsv"), encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        fields = rdr.fieldnames or []
        samples = [c for c in fields if c != "ASV"]
        table = {r["ASV"]: {s: float(r[s]) for s in samples} for r in rdr}
    tax = {r["ASV"]: r for r in _read_tsv(os.path.join(out_dir, "taxonomy.tsv"))}
    labels = load_labels(sample_sheet)
    totals = {s: 0.0 for s in samples}
    for counts in table.values():
        for s in samples:
            totals[s] += counts.get(s, 0.0)
    samples = [s for s in samples if s in labels and totals[s] > 0]
    return samples, table, tax, labels, totals


def rankdata(values: List[float]) -> List[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def mann_whitney(x: List[float], y: List[float]) -> Tuple[float, float, float]:
    """Return (auc, p_two_sided, signed_z) for demodex(x) vs control(y).

    signed_z > 0 means x>y (enriched in demodex). z is the tie-corrected
    normal-approximation statistic used later for Stouffer meta-analysis.
    """
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.5, 1.0, 0.0
    allv = x + y
    ranks = rankdata(allv)
    U1 = sum(ranks[:n1]) - n1 * (n1 + 1) / 2.0
    auc = U1 / (n1 * n2)
    mu = n1 * n2 / 2.0
    counts: Dict[float, int] = {}
    for v in allv:
        counts[v] = counts.get(v, 0) + 1
    n = n1 + n2
    tie = sum(t**3 - t for t in counts.values())
    sigma = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1)))) if n > 1 else 0.0
    if sigma == 0:
        return auc, 1.0, 0.0
    signed_z = (U1 - mu) / sigma
    p_two = math.erfc(abs(signed_z) / math.sqrt(2))
    return auc, p_two, signed_z


def genus_stats(cohort) -> Dict[str, Dict[str, float]]:
    samples, table, tax, labels, totals = cohort
    demo = [s for s in samples if labels[s] == "demodex"]
    ctrl = [s for s in samples if labels[s] == "control"]
    genera: Dict[str, List[str]] = {}
    for a, t in tax.items():
        g = (t.get("Genus") or "").strip()
        if g and g.lower() not in ("", "na"):
            genera.setdefault(g, []).append(a)
    stats: Dict[str, Dict[str, float]] = {}
    for g, ids in genera.items():
        def frac(s):
            num = sum(table[a].get(s, 0.0) for a in ids)
            return num / totals[s] if totals[s] > 0 else 0.0
        x = [frac(s) for s in demo]
        y = [frac(s) for s in ctrl]
        vals = x + y
        if sum(vals) == 0:
            continue
        auc, p, z = mann_whitney(x, y)
        prev = sum(1 for v in vals if v > 0) / len(vals)
        stats[g] = {"auc": auc, "p": p, "z": z, "n": float(len(vals)),
                    "mean_demo": sum(x) / len(x),
                    "mean_ctrl": sum(y) / len(y), "prev": prev}
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", required=True,
                    help="out_dir:sample_sheet:label")
    ap.add_argument("--min-prev", type=float, default=0.5)
    ap.add_argument("--out", default="/home/jonat/rosacea_concordance.md")
    args = ap.parse_args()

    labels_order = []
    per_cohort: Dict[str, Dict[str, Dict[str, float]]] = {}
    for spec in args.cohort:
        out_dir, ss, label = spec.split(":")
        per_cohort[label] = genus_stats(load_cohort(out_dir, ss))
        labels_order.append(label)

    # Genera present at >= min-prev prevalence in ALL cohorts
    common = None
    for label in labels_order:
        keep = {g for g, s in per_cohort[label].items() if s["prev"] >= args.min_prev}
        common = keep if common is None else (common & keep)
    common = sorted(common or [])

    lines = ["# Cross-cohort genus meta-analysis (ocular Demodex vs control)\n"]
    lines.append(f"Cohorts: {', '.join(labels_order)}. "
                 f"Prevalence filter: genus present in >={args.min_prev:.0%} of samples in EVERY cohort.\n")
    lines.append(f"Genera passing prevalence filter in all cohorts: {len(common)}.\n")
    lines.append("Meta-analysis: Stouffer combined Z of per-cohort Mann-Whitney signed z, "
                 "weighted by sqrt(n). Concordant = same sign in all cohorts. "
                 "Combined p is two-sided; Bonferroni over the tested genera.\n")

    def stouffer(g, labs=None) -> Tuple[float, float, bool]:
        labs = labs if labs is not None else labels_order
        zs = [per_cohort[l][g]["z"] for l in labs]
        ws = [math.sqrt(per_cohort[l][g]["n"]) for l in labs]
        Z = sum(w * z for w, z in zip(ws, zs)) / math.sqrt(sum(w * w for w in ws))
        p = math.erfc(abs(Z) / math.sqrt(2))
        concord = all(z > 0 for z in zs) or all(z < 0 for z in zs)
        return Z, p, concord

    m = len(common)
    scored = []
    for g in common:
        Z, p, concord = stouffer(g)
        scored.append((abs(Z), Z, p, concord, g))
    scored.sort(reverse=True)

    # Benjamini-Hochberg FDR q-values over the m tested genera.
    pvals = sorted((p, g) for _, _, p, _, g in scored)
    qmap: Dict[str, float] = {}
    prev_q = 1.0
    for rank in range(len(pvals), 0, -1):
        p_i, g_i = pvals[rank - 1]
        prev_q = min(prev_q, p_i * m / rank)
        qmap[g_i] = prev_q

    hdr = "| genus | dir | combined Z | combined p | Bonferroni p | BH q | concordant | " + \
          " | ".join(f"AUC[{l}]" for l in labels_order) + " |"
    sep = "|" + "---|" * (7 + len(labels_order))
    lines.append("\n## All prevalence-filtered genera ranked by meta-analytic |Z|\n")
    lines.append(hdr)
    lines.append(sep)
    for _, Z, p, concord, g in scored:
        direction = "up" if Z > 0 else "down"
        pbonf = min(1.0, p * m)
        aucs = " | ".join(f"{per_cohort[l][g]['auc']:.3f}" for l in labels_order)
        flag = "yes" if concord else "no"
        lines.append(f"| {g} | {direction} | {Z:+.2f} | {p:.4f} | {pbonf:.3f} | {qmap[g]:.3f} | {flag} | {aucs} |")

    lines.append("\n## Concordant hits surviving multiple-testing correction\n")
    hits = [(Z, p, g) for _, Z, p, concord, g in scored
            if concord and (p * m < 0.05 or qmap[g] < 0.05)]
    if hits:
        for Z, p, g in hits:
            lines.append(f"- **{g}** ({'up' if Z>0 else 'down'} in Demodex): "
                         f"combined Z={Z:+.2f}, p={p:.4g}, Bonferroni p={min(1.0,p*m):.4g}, "
                         f"BH q={qmap[g]:.4g}")
    else:
        lines.append("_None survive Bonferroni/BH at the current cohort count._")

    # Leave-one-cohort-out robustness for the concordant leads.
    if len(labels_order) >= 3:
        lines.append("\n## Leave-one-cohort-out robustness (top concordant leads)\n")
        lines.append("Combined Z (and two-sided p) recomputed after dropping each cohort. "
                     "A robust signal keeps its sign and comparable magnitude throughout.\n")
        leads = [g for _, _, _, concord, g in scored if concord][:6]
        lines.append("| genus | full Z | " +
                     " | ".join(f"drop {l}" for l in labels_order) + " |")
        lines.append("|" + "---|" * (2 + len(labels_order)))
        for g in leads:
            Zf, _, _ = stouffer(g)
            cells = []
            for drop in labels_order:
                labs = [l for l in labels_order if l != drop]
                Zd, pd, _ = stouffer(g, labs)
                cells.append(f"{Zd:+.2f} (p={pd:.3f})")
            lines.append(f"| {g} | {Zf:+.2f} | " + " | ".join(cells) + " |")

    report = "\n".join(lines) + "\n"
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
