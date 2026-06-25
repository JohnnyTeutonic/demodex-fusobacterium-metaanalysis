"""
E6c -- sequence-identity sanity check for the k-mer cross-mapping claim.

For the top ASVs that the Snodgrassella_alvi k-mer panel flags, compute their
best ungapped percent identity against our S. alvi 16S reference (both strands,
best sliding window). V4 ASVs (~253 bp) align to 16S without indels over the
region, so ungapped identity is a fair quick estimate.

A true S. alvi variant would be ~99-100% identical; a cross-mapping Neisseriaceae
relative will sit far lower (~80-90%), proving the k-mer hit is conserved-region
overlap, not genuine S. alvi.

Usage: python identity_check.py [ASV43 ASV1109 ...]
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(__file__)
REPO = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
REF = os.path.join(REPO, "data", "refs", "Snodgrassella_alvi.fasta")
ASV_FASTA = os.path.join(REPO, "data", "dada2_PRJNA692647", "asv_seqs.fasta")

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = ln[1:].split()[0], []
            else:
                buf.append(ln)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def best_ungapped_identity(query: str, ref: str) -> float:
    """Best identity of query placed ungapped against ref (slide over ref)."""
    q, best = query.upper(), 0.0
    L = len(q)
    for r in (ref.upper(), revcomp(ref.upper())):
        for off in range(0, len(r) - L + 1):
            m = sum(1 for a, b in zip(q, r[off:off + L]) if a == b)
            if m / L > best:
                best = m / L
    return best


def main() -> None:
    asv_ids = sys.argv[1:] or ["ASV43", "ASV1109", "ASV930", "ASV591", "ASV1804"]
    ref_seqs = read_fasta(REF)
    ref = max(ref_seqs.values(), key=len)  # longest S. alvi reference
    asvs = read_fasta(ASV_FASTA)
    print(f"S. alvi reference length: {len(ref)} bp; comparing {len(asv_ids)} ASVs\n")
    for a in asv_ids:
        if a not in asvs:
            print(f"{a}: not found"); continue
        ident = best_ungapped_identity(asvs[a], ref)
        verdict = ("LIKELY S. alvi" if ident >= 0.97 else
                   "NOT S. alvi (cross-mapping relative)" if ident < 0.93 else
                   "ambiguous")
        print(f"{a}: len={len(asvs[a])}  best identity to S. alvi 16S = {ident*100:.1f}%  -> {verdict}")


if __name__ == "__main__":
    main()
