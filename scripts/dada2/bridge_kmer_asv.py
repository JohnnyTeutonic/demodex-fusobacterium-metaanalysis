"""
E6b -- bridge the k-mer index and the DADA2 ASVs.

Applies the SAME diagnostic k-mer panel (data/kmers/<taxon>.txt) used by the
read-level index directly to the denoised ASV sequences. For each target taxon
it reports which ASVs carry >=2 diagnostic k-mers (the index's two-hit rule),
how SILVA classified those ASVs, their abundance, and the Demodex-vs-control
separation at that exact-sequence level.

Interpretation:
  - If diagnostic k-mers land on a real, abundant ASV  -> the taxon IS present;
    any "genus not detected" in SILVA is a taxonomy-resolution gap, not absence.
  - If diagnostic k-mers hit no ASV (or only ultra-rare/singleton ASVs)
    -> the read-level k-mer signal did not assemble into a real variant,
       i.e. it was error/cross-mapping noise (artefact).

Pure stdlib. Reuses the Mann-Whitney/AUC from the E6 analysis logic.

Usage:
  python bridge_kmer_asv.py
  python bridge_kmer_asv.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from typing import Dict, List, Set, Tuple

HERE = os.path.dirname(__file__)
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
KMER_DIR = os.path.join(REPO, "data", "kmers")
DADA_DIR = os.path.join(REPO, "data", "dada2_PRJNA692647")
TARGETS = ["Snodgrassella_alvi", "Bacillus_oleronius", "Bartonella_quintana"]
MIN_HITS = 2  # two-hit attribution rule, matching build_index.py

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def kmerize(seq: str, k: int) -> Set[str]:
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def load_kmers(taxon: str) -> Tuple[Set[str], int]:
    path = os.path.join(KMER_DIR, f"{taxon}.txt")
    with open(path, encoding="utf-8") as fh:
        kms = {ln.strip().upper() for ln in fh if ln.strip()}
    k = len(next(iter(kms))) if kms else 31
    return kms, k


def read_fasta(path: str) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    name = None
    buf: List[str] = []
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = ln[1:].split()[0]
                buf = []
            else:
                buf.append(ln)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def _read_tsv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        return rdr.fieldnames or [], list(rdr)


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


def mann_whitney(x: List[float], y: List[float]) -> Tuple[float, float]:
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return 0.5, 1.0
    allv = x + y
    ranks = rankdata(allv)
    U1 = sum(ranks[:n1]) - n1 * (n1 + 1) / 2.0
    auc = U1 / (n1 * n2)
    counts: Dict[float, int] = {}
    for v in allv:
        counts[v] = counts.get(v, 0) + 1
    n = n1 + n2
    tie = sum(t**3 - t for t in counts.values())
    var = n1 * n2 / 12.0 * ((n + 1) - tie / (n * (n - 1))) if n > 1 else 0.0
    if var <= 0:
        return auc, 1.0
    z = (U1 - n1 * n2 / 2.0 - 0.5) / math.sqrt(var)
    return auc, 0.5 * math.erfc(z / math.sqrt(2))


def bridge() -> str:
    asv_seqs = read_fasta(os.path.join(DADA_DIR, "asv_seqs.fasta"))
    fields, tax_rows = _read_tsv(os.path.join(DADA_DIR, "taxonomy.tsv"))
    tax = {r["ASV"]: r for r in tax_rows}
    afields, arows = _read_tsv(os.path.join(DADA_DIR, "asv_table.tsv"))
    samples = [c for c in afields if c != "ASV"]
    counts = {r["ASV"]: {s: float(r[s]) for s in samples} for r in arows}
    _, ss = _read_tsv(os.path.join(DADA_DIR, "sample_sheet.tsv"))
    labels = {r["run_accession"]: r["group"] for r in ss if r.get("group") in ("demodex", "control")}
    samples = [s for s in samples if s in labels]
    demo = [s for s in samples if labels[s] == "demodex"]
    ctrl = [s for s in samples if labels[s] == "control"]
    totals = {s: sum(counts[a].get(s, 0.0) for a in counts) for s in samples}

    # precompute ASV k-mer sets (both strands)
    asv_kmers31: Dict[str, Set[str]] = {}

    lines = ["# E6b -- k-mer panel applied directly to DADA2 ASVs\n",
             f"{len(asv_seqs)} ASVs; {len(demo)} demodex vs {len(ctrl)} control.\n"]

    for taxon in TARGETS:
        kms, k = load_kmers(taxon)
        if not asv_kmers31 or next(iter(asv_kmers31.values())) is None:
            pass
        hits = []
        for asv, seq in asv_seqs.items():
            kk = asv_kmers31.get((asv, k))
            if kk is None:
                kk = kmerize(seq, k) | kmerize(revcomp(seq), k)
                asv_kmers31[(asv, k)] = kk
            nhit = len(kk & kms)
            if nhit >= MIN_HITS:
                hits.append((asv, nhit))
        lines.append(f"\n## {taxon} (panel: {len(kms)} diagnostic {k}-mers)\n")
        if not hits:
            lines.append("- **NO ASV carries >=2 diagnostic k-mers** -> not present as a "
                         "real variant; the read-level signal did not assemble (artefact).")
            continue
        hits.sort(key=lambda t: sum(counts.get(t[0], {}).values()), reverse=True)
        for asv, nhit in hits:
            t = tax.get(asv, {})
            fam = t.get("Family", "?"); gen = t.get("Genus", "NA"); sp = t.get("Species", "NA")
            tot = sum(counts.get(asv, {}).values())
            def frac(s): return counts[asv].get(s, 0.0) / totals[s] if totals[s] > 0 else 0.0
            auc, p = mann_whitney([frac(s) for s in demo], [frac(s) for s in ctrl])
            lines.append(f"- **{asv}**: {nhit} diagnostic k-mers, total reads={int(tot)}, "
                         f"SILVA=[Family={fam}; Genus={gen}; Species={sp}], "
                         f"demodex-vs-control AUC={auc:.3f} (p={p:.3f}).")
    report = "\n".join(lines) + "\n"
    out = os.path.join(DADA_DIR, "E6b_kmer_bridge.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(report)
    print(report)
    print(f"Wrote {out}")
    return report


def _self_test() -> None:
    assert revcomp("AAAACCCG") == "CGGGTTTT"
    assert kmerize("ACGTACGT", 4) == {"ACGT", "CGTA", "GTAC", "TACG"}
    # a sequence containing 2 known kmers should be detected
    seq = "TTT" + "ACGTACGTACGT"
    ks = kmerize(seq, 4) | kmerize(revcomp(seq), 4)
    assert len(ks & {"ACGT", "CGTA"}) == 2
    auc, _ = mann_whitney([3, 4, 5], [0, 1, 2])
    assert abs(auc - 1.0) < 1e-9
    print("self-test OK: revcomp/kmerize/mann_whitney pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    bridge()


if __name__ == "__main__":
    main()
