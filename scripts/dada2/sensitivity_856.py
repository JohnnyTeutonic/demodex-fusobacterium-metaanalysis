"""
Sensitivity / confounding check for the shotgun cohort (PRJNA856121) Fusobacterium
depletion result. Ocular shotgun samples are ~99% host DNA, so we must verify the
Fusobacterium-down signal is not an artifact of group differences in sequencing
depth or host-DNA fraction.

Reads Kraken2 reports (~/rosacea_856/kraken/reports/*.kreport) + sample_sheet.tsv.
Per sample extracts: total reads, unclassified, Homo (9606) clade, Bacteria (2)
clade, Fusobacterium (genus) clade.

Reports (Demodex vs control):
  1. sequencing depth (total reads)                  -- should NOT differ by group
  2. host-DNA fraction (Homo / total)                -- should NOT differ by group
  3. bacterial read depth (Bacteria clade)           -- should NOT differ by group
  4. Fusobacterium depletion under 3 normalisations:
       (a) fraction of bacterial reads   (compositional -- the meta-analysis unit)
       (b) fraction of total reads        (depth/host independent numerator)
       (c) reads-per-million total        (absolute-ish)
       (d) presence/absence prevalence
  5. Spearman rho between Fusobacterium(a) and host-fraction / bacterial-depth,
     to check the depletion is not merely tracking a technical covariate.

Pure stdlib.
"""
from __future__ import annotations

import csv
import glob
import math
import os

WORK = os.path.expanduser("~/rosacea_856")
REPORTS = os.path.join(WORK, "kraken", "reports")


def parse_report(path: str) -> dict:
    """Return dict with clade read counts for key taxa."""
    unclassified = 0
    root = 0
    homo = 0
    bacteria = 0
    fuso = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            clade = int(f[1])
            rank = f[3].strip()
            taxid = f[4].strip()
            name = f[5].strip()
            if rank == "U":
                unclassified = clade
            elif taxid == "1":
                root = clade
            elif taxid == "9606":
                homo = clade
            elif taxid == "2" and rank == "D":
                bacteria = clade
            elif rank == "G" and name == "Fusobacterium":
                fuso = clade
    total = unclassified + root
    return {
        "total": total, "unclassified": unclassified, "homo": homo,
        "bacteria": bacteria, "fuso": fuso,
    }


def mannwhitney_u(a: list, b: list):
    """Two-sided Mann-Whitney U with normal approx; returns (U, z, p, auc)."""
    n1, n2 = len(a), len(b)
    combined = [(v, 0) for v in a] + [(v, 1) for v in b]
    combined.sort(key=lambda x: x[0])
    ranks = [0.0] * len(combined)
    i = 0
    while i < len(combined):
        j = i
        while j + 1 < len(combined) and combined[j + 1][0] == combined[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[k] = avg
        i = j + 1
    r1 = sum(ranks[k] for k in range(len(combined)) if combined[k][1] == 0)
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u = u1
    mu = n1 * n2 / 2.0
    # tie-corrected sigma
    from collections import Counter
    tie = Counter(v for v, _ in combined)
    tt = sum(t ** 3 - t for t in tie.values())
    n = n1 + n2
    sigma2 = (n1 * n2 / 12.0) * ((n + 1) - tt / (n * (n - 1)))
    sigma = math.sqrt(sigma2) if sigma2 > 0 else 0.0
    z = (u - mu) / sigma if sigma > 0 else 0.0
    p = math.erfc(abs(z) / math.sqrt(2))
    auc = u1 / (n1 * n2)  # AUC for group a (Demodex) being LARGER
    return u, z, p, auc


def spearman(x: list, y: list) -> float:
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else 0.0


def med(v: list) -> float:
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2.0


def main() -> None:
    groups: dict = {}
    with open(os.path.join(WORK, "sample_sheet.tsv"), encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            groups[r["run_accession"]] = r["group"]

    rows = {}
    for path in glob.glob(os.path.join(REPORTS, "*.kreport")):
        run = os.path.basename(path)[:-len(".kreport")]
        if run in groups:
            rows[run] = parse_report(path)

    db = [r for r in rows if groups[r] == "demodex"]
    ct = [r for r in rows if groups[r] == "control"]
    print(f"samples: demodex={len(db)} control={len(ct)}\n")

    def col(runs, key):
        return [rows[r][key] for r in runs]

    def report_metric(name, db_vals, ct_vals, higher_bad=False):
        u, z, p, auc = mannwhitney_u(db_vals, ct_vals)
        print(f"{name}")
        print(f"    demodex median = {med(db_vals):.4g} | control median = {med(ct_vals):.4g}")
        print(f"    MWU p = {p:.4f} | AUC(demodex>control) = {auc:.3f}")
        return p, auc

    print("=== 1-3. Technical covariates (want p>0.05 = NO group difference) ===")
    report_metric("total reads (depth)", col(db, "total"), col(ct, "total"))
    hf_db = [rows[r]["homo"] / rows[r]["total"] for r in db]
    hf_ct = [rows[r]["homo"] / rows[r]["total"] for r in ct]
    report_metric("host fraction (Homo/total)", hf_db, hf_ct)
    report_metric("bacterial reads (Bacteria clade)", col(db, "bacteria"), col(ct, "bacteria"))

    print("\n=== 4. Fusobacterium depletion under multiple normalisations ===")
    # (a) fraction of bacterial reads
    fa_db = [rows[r]["fuso"] / rows[r]["bacteria"] if rows[r]["bacteria"] else 0.0 for r in db]
    fa_ct = [rows[r]["fuso"] / rows[r]["bacteria"] if rows[r]["bacteria"] else 0.0 for r in ct]
    report_metric("(a) Fuso / bacterial reads", fa_db, fa_ct)
    # (b) fraction of total reads (independent of host/bacterial split)
    fb_db = [rows[r]["fuso"] / rows[r]["total"] for r in db]
    fb_ct = [rows[r]["fuso"] / rows[r]["total"] for r in ct]
    report_metric("(b) Fuso / total reads", fb_db, fb_ct)
    # (c) reads per million total
    fc_db = [1e6 * rows[r]["fuso"] / rows[r]["total"] for r in db]
    fc_ct = [1e6 * rows[r]["fuso"] / rows[r]["total"] for r in ct]
    report_metric("(c) Fuso reads-per-million total", fc_db, fc_ct)
    # (d) presence/absence
    pa_db = sum(1 for r in db if rows[r]["fuso"] > 0)
    pa_ct = sum(1 for r in ct if rows[r]["fuso"] > 0)
    print("(d) Fusobacterium detected (reads>0)")
    print(f"    demodex: {pa_db}/{len(db)} ({100*pa_db/len(db):.0f}%) | control: {pa_ct}/{len(ct)} ({100*pa_ct/len(ct):.0f}%)")

    print("\n=== 5. Is depletion just tracking a technical covariate? (Spearman) ===")
    all_runs = db + ct
    fuso_a = [rows[r]["fuso"] / rows[r]["bacteria"] if rows[r]["bacteria"] else 0.0 for r in all_runs]
    hostf = [rows[r]["homo"] / rows[r]["total"] for r in all_runs]
    bact = [rows[r]["bacteria"] for r in all_runs]
    print(f"    rho(Fuso_frac_bact, host_fraction)   = {spearman(fuso_a, hostf):+.3f}")
    print(f"    rho(Fuso_frac_bact, bacterial_depth) = {spearman(fuso_a, bact):+.3f}")

    print("\nInterpretation: if (1-3) show p>0.05 (no depth/host group difference) and")
    print("Fusobacterium is depleted under BOTH (a) and (b) normalisations, the signal")
    print("is not a compositional/host artifact.")


if __name__ == "__main__":
    main()
