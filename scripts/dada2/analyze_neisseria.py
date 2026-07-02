"""
Genus/family-level replication test for the ocular Demodex *Neisseria* signal.

Reads a DADA2 ASV table + SILVA taxonomy + sample sheet (run_accession -> group in
{demodex, control}), collapses to target genera and to the Neisseriaceae family,
and tests Demodex > control separation with the SAME statistics used for the
discovery cohort (Mann-Whitney U, one-sided demodex>control, AUC = U/(n1*n2)).

Because different cohorts amplify different 16S regions (V4 515F/806R vs V3-V4
338F/806R), ASV sequences are NOT comparable across cohorts; the region-robust
unit of replication is the GENUS (and the Neisseriaceae FAMILY). This script
therefore reports genus- and family-level tests plus a full ranked table of the
strongest genus-level movers, so the honest full picture is visible.

Pure stdlib. Usage:
  python analyze_neisseria.py --out-dir ~/rosacea_657256/out --sample-sheet ~/rosacea_657256/sample_sheet.tsv --label PRJNA657256
  python analyze_neisseria.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Dict, List, Tuple

TARGET_GENERA = ["Neisseria", "Snodgrassella", "Eikenella", "Kingella"]
TARGET_FAMILY = "Neisseriaceae"


def _read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
        return rdr.fieldnames or [], rows


def load_labels(sample_sheet: str) -> Dict[str, str]:
    _, rows = _read_tsv(sample_sheet)
    return {r["run_accession"]: r["group"] for r in rows if r.get("group") in ("demodex", "control")}


def load_asv_table(path: str) -> Tuple[List[str], Dict[str, Dict[str, float]]]:
    fields, rows = _read_tsv(path)
    samples = [c for c in fields if c != "ASV"]
    table: Dict[str, Dict[str, float]] = {}
    for r in rows:
        table[r["ASV"]] = {s: float(r[s]) for s in samples}
    return samples, table


def load_taxonomy(path: str) -> Dict[str, Dict[str, str]]:
    _, rows = _read_tsv(path)
    return {r["ASV"]: r for r in rows}


def totals_per_sample(samples: List[str], table: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    totals = {s: 0.0 for s in samples}
    for counts in table.values():
        for s in samples:
            totals[s] += counts.get(s, 0.0)
    return totals


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
    """Return (U_x, auc, p_one_sided) testing x>y (x = demodex)."""
    n1, n2 = len(x), len(y)
    allv = x + y
    ranks = rankdata(allv)
    R1 = sum(ranks[:n1])
    U1 = R1 - n1 * (n1 + 1) / 2.0
    auc = U1 / (n1 * n2)
    mu = n1 * n2 / 2.0
    counts: Dict[float, int] = {}
    for v in allv:
        counts[v] = counts.get(v, 0) + 1
    n = n1 + n2
    tie = sum(t**3 - t for t in counts.values())
    sigma = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1))))
    if sigma == 0:
        return U1, auc, 1.0
    z = (U1 - mu - 0.5) / sigma
    p_one = 0.5 * math.erfc(z / math.sqrt(2))
    return U1, auc, p_one


def analyze(out_dir: str, sample_sheet: str, label: str) -> str:
    out_dir = os.path.expanduser(out_dir)
    labels = load_labels(os.path.expanduser(sample_sheet))
    samples, table = load_asv_table(os.path.join(out_dir, "asv_table.tsv"))
    tax = load_taxonomy(os.path.join(out_dir, "taxonomy.tsv"))
    totals = totals_per_sample(samples, table)

    samples = [s for s in samples if s in labels and totals[s] > 0]
    demo = [s for s in samples if labels[s] == "demodex"]
    ctrl = [s for s in samples if labels[s] == "control"]

    def frac(asv_ids: List[str], s: str) -> float:
        num = sum(table[a].get(s, 0.0) for a in asv_ids)
        return num / totals[s] if totals[s] > 0 else 0.0

    lines: List[str] = []
    lines.append(f"# Neisseria replication test -- {label} (DADA2 + SILVA 138.1)\n")
    lines.append(f"Samples: {len(demo)} demodex vs {len(ctrl)} control. Total ASVs: {len(table)}.\n")

    def feature_test(asv_ids: List[str], name: str) -> Tuple[float, float, float, float]:
        x = [frac(asv_ids, s) for s in demo]
        y = [frac(asv_ids, s) for s in ctrl]
        md = sum(x) / len(x) if x else 0.0
        mc = sum(y) / len(y) if y else 0.0
        if sum(x) + sum(y) == 0:
            lines.append(f"- **{name}**: not detected.")
            return (0.5, 1.0, md, mc)
        _, auc, p = mann_whitney(x, y)
        lines.append(f"- **{name}**: AUC={auc:.3f}, one-sided p={p:.4f}; "
                     f"mean rel.abund demodex={md:.2e} vs control={mc:.2e} "
                     f"({len(asv_ids)} ASV(s)).")
        return (auc, p, md, mc)

    # Family composite
    lines.append("\n## Neisseriaceae family composite (pre-specified)\n")
    fam_ids = [a for a, t in tax.items() if (t.get("Family") or "").strip() == TARGET_FAMILY]
    feature_test(fam_ids, f"{TARGET_FAMILY} (family)")

    # Target genera
    lines.append("\n## Target genera (pre-specified)\n")
    for genus in TARGET_GENERA:
        ids = [a for a, t in tax.items() if (t.get("Genus") or "").strip() == genus]
        feature_test(ids, f"{genus} (genus)")

    # Full genus differential ranking (context; not the pre-specified test)
    lines.append("\n## All genera ranked by demodex>control effect (context)\n")
    genera: Dict[str, List[str]] = {}
    for a, t in tax.items():
        g = (t.get("Genus") or "").strip()
        if g:
            genera.setdefault(g, []).append(a)
    scored = []
    for g, ids in genera.items():
        x = [frac(ids, s) for s in demo]
        y = [frac(ids, s) for s in ctrl]
        if sum(x) + sum(y) == 0:
            continue
        _, auc, p = mann_whitney(x, y)
        md = sum(x) / len(x) if x else 0.0
        prevalence = sum(1 for v in (x + y) if v > 0) / (len(x) + len(y))
        scored.append((auc, p, g, md, prevalence))
    scored.sort(key=lambda r: r[0], reverse=True)
    lines.append("Top 15 enriched in demodex (AUC high):")
    lines.append("| genus | AUC | one-sided p | mean RA (demodex) | prevalence |")
    lines.append("|---|---|---|---|---|")
    for auc, p, g, md, prev in scored[:15]:
        lines.append(f"| {g} | {auc:.3f} | {p:.4f} | {md:.2e} | {prev:.2f} |")

    report = "\n".join(lines) + "\n"
    rpt_path = os.path.join(out_dir, f"neisseria_replication_{label}.md")
    with open(rpt_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    print(f"Wrote {rpt_path}")
    return report


def _self_test() -> None:
    assert rankdata([10, 20, 20, 40]) == [1.0, 2.5, 2.5, 4.0]
    _, auc, p = mann_whitney([5, 6, 7], [1, 2, 3])
    assert abs(auc - 1.0) < 1e-9, auc
    assert p < 0.1, p
    _, auc2, _ = mann_whitney([1, 2, 3], [5, 6, 7])
    assert abs(auc2 - 0.0) < 1e-9, auc2
    _, auc3, _ = mann_whitney([1, 2, 3, 4], [1, 2, 3, 4])
    assert abs(auc3 - 0.5) < 1e-9, auc3
    print("self-test OK")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default="~/rosacea_657256/out")
    ap.add_argument("--sample-sheet", default="~/rosacea_657256/sample_sheet.tsv")
    ap.add_argument("--label", default="PRJNA657256")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    analyze(args.out_dir, args.sample_sheet, args.label)


if __name__ == "__main__":
    main()
