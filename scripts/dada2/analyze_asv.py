"""
E6 -- ASV-level re-analysis of the PRJNA692647 validation cohort (DADA2 output).

Reads the DADA2 ASV table + SILVA taxonomy + sample sheet, collapses to the
target taxa (Snodgrassella, Bartonella, Bacillus) at genus AND species level,
and re-tests Demodex+ vs control separation with the same statistics used for
the k-mer index (Mann-Whitney U, one-sided demodex>control, AUC = U/(n1*n2)).

This directly answers the E5 limitation: does exact-sequence (ASV) resolution
sharpen the genus-level k-mer signal -- in particular, is the Snodgrassella
signal carried by an exact *S. alvi* variant?

Pure stdlib (no scipy/numpy). Inputs default to the DADA2 out dir.

Usage:
  python analyze_asv.py --out-dir ~/rosacea_dada2/out --sample-sheet ~/rosacea_dada2/sample_sheet.tsv
  python analyze_asv.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Dict, List, Tuple

TARGET_GENERA = ["Snodgrassella", "Bartonella", "Bacillus"]


def _read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        rows = list(rdr)
        return rdr.fieldnames or [], rows


def load_labels(sample_sheet: str) -> Dict[str, str]:
    _, rows = _read_tsv(sample_sheet)
    return {r["run_accession"]: r["group"] for r in rows if r.get("group") in ("demodex", "control")}


def load_asv_table(path: str) -> Tuple[List[str], Dict[str, Dict[str, float]]]:
    """Return (sample_ids, {asv_id: {sample: count}})."""
    fields, rows = _read_tsv(path)
    samples = [c for c in fields if c != "ASV"]
    table: Dict[str, Dict[str, float]] = {}
    for r in rows:
        table[r["ASV"]] = {s: float(r[s]) for s in samples}
    return samples, table


def load_taxonomy(path: str) -> Dict[str, Dict[str, str]]:
    _, rows = _read_tsv(path)
    return {r["ASV"]: r for r in rows}


def relabundance(samples: List[str], table: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    """Total counts per sample (for normalisation)."""
    totals = {s: 0.0 for s in samples}
    for counts in table.values():
        for s in samples:
            totals[s] += counts.get(s, 0.0)
    return totals


def rankdata(values: List[float]) -> List[float]:
    """Average ranks with tie handling."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0  # 1-based average rank
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def mann_whitney(x: List[float], y: List[float]) -> Tuple[float, float, float]:
    """
    Return (U_x, auc, p_one_sided) testing x>y (x = demodex).
    Normal approximation with tie correction; AUC = U_x/(n1*n2).
    """
    n1, n2 = len(x), len(y)
    allv = x + y
    ranks = rankdata(allv)
    R1 = sum(ranks[:n1])
    U1 = R1 - n1 * (n1 + 1) / 2.0
    auc = U1 / (n1 * n2)
    mu = n1 * n2 / 2.0
    # tie correction
    counts: Dict[float, int] = {}
    for v in allv:
        counts[v] = counts.get(v, 0) + 1
    n = n1 + n2
    tie = sum(t**3 - t for t in counts.values())
    sigma = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1))))
    if sigma == 0:
        return U1, auc, 1.0
    z = (U1 - mu - 0.5) / sigma  # continuity-corrected, testing U1 large => x>y
    p_one = 0.5 * math.erfc(z / math.sqrt(2))
    return U1, auc, p_one


def species_of(tax: Dict[str, str]) -> str:
    g = (tax.get("Genus") or "").strip()
    sp = (tax.get("Species") or "").strip()
    if g and sp:
        return f"{g} {sp}"
    return g or "(unassigned)"


def analyze(out_dir: str, sample_sheet: str) -> str:
    out_dir = os.path.expanduser(out_dir)
    labels = load_labels(os.path.expanduser(sample_sheet))
    samples, table = load_asv_table(os.path.join(out_dir, "asv_table.tsv"))
    tax = load_taxonomy(os.path.join(out_dir, "taxonomy.tsv"))
    totals = relabundance(samples, table)

    samples = [s for s in samples if s in labels]
    demo = [s for s in samples if labels[s] == "demodex"]
    ctrl = [s for s in samples if labels[s] == "control"]

    lines: List[str] = []
    lines.append("# E6 -- ASV-level validation (PRJNA692647, DADA2 + SILVA 138.1)\n")
    lines.append(f"Samples: {len(demo)} demodex vs {len(ctrl)} control. "
                 f"Total ASVs: {len(table)}.\n")

    def feature_test(asv_ids: List[str], name: str) -> None:
        def frac(s: str) -> float:
            num = sum(table[a].get(s, 0.0) for a in asv_ids)
            return num / totals[s] if totals[s] > 0 else 0.0
        x = [frac(s) for s in demo]
        y = [frac(s) for s in ctrl]
        if sum(x) + sum(y) == 0:
            lines.append(f"- **{name}**: not detected.")
            return
        _, auc, p = mann_whitney(x, y)
        md = sum(x) / len(x) if x else 0
        mc = sum(y) / len(y) if y else 0
        lines.append(f"- **{name}**: AUC={auc:.3f}, one-sided p={p:.4f}; "
                     f"mean rel.abund demodex={md:.2e} vs control={mc:.2e} "
                     f"({len(asv_ids)} ASV(s)).")

    # Genus-level composites
    lines.append("\n## Genus-level composites\n")
    for genus in TARGET_GENERA:
        ids = [a for a, t in tax.items() if (t.get("Genus") or "").strip() == genus]
        feature_test(ids, f"{genus} (genus)")

    # Per-ASV for target genera (species-resolved)
    lines.append("\n## Per-ASV resolution within target genera\n")
    for genus in TARGET_GENERA:
        ids = [a for a, t in tax.items() if (t.get("Genus") or "").strip() == genus]
        if not ids:
            continue
        # rank target ASVs by total abundance
        ids.sort(key=lambda a: sum(table[a].values()), reverse=True)
        for a in ids:
            feature_test([a], f"{a} [{species_of(tax[a])}]")

    report = "\n".join(lines) + "\n"
    rpt_path = os.path.join(out_dir, "E6_asv_validation.md")
    with open(rpt_path, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    print(f"Wrote {rpt_path}")
    return report


def _self_test() -> None:
    # rankdata ties
    assert rankdata([10, 20, 20, 40]) == [1.0, 2.5, 2.5, 4.0]
    # perfect separation: x all above y -> AUC 1.0
    U, auc, p = mann_whitney([5, 6, 7], [1, 2, 3])
    assert abs(auc - 1.0) < 1e-9, auc
    assert p < 0.1, p
    # reversed -> AUC 0
    _, auc2, _ = mann_whitney([1, 2, 3], [5, 6, 7])
    assert abs(auc2 - 0.0) < 1e-9, auc2
    # equal medians -> AUC ~0.5
    _, auc3, _ = mann_whitney([1, 2, 3, 4], [1, 2, 3, 4])
    assert abs(auc3 - 0.5) < 1e-9, auc3
    assert species_of({"Genus": "Snodgrassella", "Species": "alvi"}) == "Snodgrassella alvi"
    assert species_of({"Genus": "Bacillus", "Species": ""}) == "Bacillus"
    print("self-test OK: rankdata/mann_whitney/species_of pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default="~/rosacea_dada2/out")
    ap.add_argument("--sample-sheet", default="~/rosacea_dada2/sample_sheet.tsv")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    analyze(args.out_dir, args.sample_sheet)


if __name__ == "__main__":
    main()
