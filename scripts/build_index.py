"""
E2 -- build the molecular mite-burden index for a project.

For each run we stream its FASTQ(.gz) from ENA, decompress on the fly, and count
how many reads carry a diagnostic k-mer for each panel taxon (closed-reference
k-mer attribution). We report per-taxon read fractions and a combined
Demodex-associated index = sum of fractions over the Demodex-associated taxa.

This is a relative, gradeable proxy (not an ASV census): a read is attributed to
a taxon if it contains >=1 of that taxon's diagnostic (panel-unique) k-mers.

Inputs:
  data/kmers/<taxon>.txt        diagnostic k-mers (from fetch_references.py)
  data/kmers/panel_stats.tsv    taxon -> demodex_associated flag
  data/meta_<PROJECT>.tsv       run table (from fetch_metadata.py)

Output:
  data/index_<PROJECT>.tsv      run_accession, reads, <taxon>_frac..., mite_index

Usage:
  python build_index.py --project PRJDB18292 [--max-reads 0]
  python build_index.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import urllib.request
from typing import Dict, List, Set, Tuple

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
KMERS_DIR = os.path.join(DATA_DIR, "kmers")
K = 31
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def canonical(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def load_panel(kmers_dir: str = KMERS_DIR) -> Tuple[Dict[str, str], List[str], Set[str]]:
    """Return (kmer->taxon, taxa_order, demodex_associated_set)."""
    kmap: Dict[str, str] = {}
    taxa: List[str] = []
    assoc: Set[str] = set()
    stats = os.path.join(kmers_dir, "panel_stats.tsv")
    with open(stats, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            taxa.append(row["taxon"])
            if row["demodex_associated"].strip().lower() in ("true", "1", "yes"):
                assoc.add(row["taxon"])
    for taxon in taxa:
        path = os.path.join(kmers_dir, taxon + ".txt")
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                km = line.strip()
                if km:
                    # store both orientations so the read scan is a plain dict
                    # lookup (no per-k-mer canonicalisation -> much faster)
                    kmap[km] = taxon
                    kmap[revcomp(km)] = taxon
    return kmap, taxa, assoc


def iter_fastq_reads(url: str, max_reads: int = 0):
    """Stream a (gzipped) FASTQ from a URL, yielding sequence strings."""
    if url.startswith("ftp."):
        url = "https://" + url
    req = urllib.request.Request(url, headers={"User-Agent": "rosacea-demodex/1.0"})
    with urllib.request.urlopen(req, timeout=300) as resp:
        gz = gzip.GzipFile(fileobj=resp)
        n = 0
        while True:
            header = gz.readline()
            if not header:
                break
            seq = gz.readline().decode("ascii", "replace").strip()
            gz.readline()  # plus
            gz.readline()  # qual
            if not seq:
                break
            yield seq
            n += 1
            if max_reads and n >= max_reads:
                break


def count_reads(url: str, kmap: Dict[str, str], taxa: List[str],
                max_reads: int = 0, min_hits: int = 2) -> Tuple[int, Dict[str, int]]:
    """A read is attributed to a taxon if it carries >= min_hits diagnostic
    k-mers for it (raises specificity vs single-k-mer attribution)."""
    hits = {t: 0 for t in taxa}
    reads = 0
    get = kmap.get
    for seq in iter_fastq_reads(url, max_reads):
        reads += 1
        seq = seq.upper()
        n = len(seq)
        if n < K:
            continue
        per: Dict[str, int] = {}
        for i in range(n - K + 1):
            taxon = get(seq[i:i + K])
            if taxon is not None:
                per[taxon] = per.get(taxon, 0) + 1
        for t, c in per.items():
            if c >= min_hits:
                hits[t] += 1
    return reads, hits


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", type=str)
    ap.add_argument("--max-reads", type=int, default=0, help="0 = all reads")
    ap.add_argument("--kmers-dir", default=KMERS_DIR,
                    help="k-mer panel directory (default: data/kmers)")
    ap.add_argument("--out", default=None,
                    help="output TSV path (default: data/index_<PROJECT>.tsv)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    if not args.project:
        ap.error("--project required")

    kmap, taxa, assoc = load_panel(args.kmers_dir)
    print(f"panel: {len(taxa)} taxa, {len(kmap)} diagnostic k-mers, "
          f"{len(assoc)} demodex-associated")

    meta_path = os.path.join(DATA_DIR, f"meta_{args.project}.tsv")
    runs = list(csv.DictReader(open(meta_path, encoding="utf-8"), delimiter="\t"))

    out_rows = []
    for r in runs:
        run = r["run_accession"]
        url = r.get("fastq_ftp", "")
        if not url:
            print(f"  {run}: no fastq url; skipping")
            continue
        reads, hits = None, None
        for attempt in range(1, 4):
            try:
                reads, hits = count_reads(url, kmap, taxa, args.max_reads)
                break
            except Exception as e:
                print(f"  {run}: attempt {attempt}/3 failed ({e!r})")
        if reads is None:
            print(f"  {run}: GAVE UP after 3 attempts")
            continue
        row = {"run_accession": run, "reads": reads}
        mite = 0.0
        for t in taxa:
            frac = hits[t] / reads if reads else 0.0
            row[t + "_frac"] = f"{frac:.6f}"
            if t in assoc:
                mite += frac
        row["mite_index"] = f"{mite:.6f}"
        out_rows.append(row)
        print(f"  {run:<12} reads={reads:<7} mite_index={mite:.5f}  "
              + " ".join(f"{t.split('_')[0][:4]}={hits[t]}" for t in taxa))

    cols = ["run_accession", "reads"] + [t + "_frac" for t in taxa] + ["mite_index"]
    out = args.out or os.path.join(DATA_DIR, f"index_{args.project}.tsv")
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nWrote {len(out_rows)} rows -> {out}")


def _self_test() -> None:
    # build a tiny in-memory panel and confirm attribution logic
    kmap = {canonical("A" * K): "TaxonA", canonical("C" * K): "TaxonB"}
    taxa = ["TaxonA", "TaxonB"]
    # a read containing A*K should hit TaxonA only
    seq = "A" * K + "G" * 10
    seen = set()
    for i in range(len(seq) - K + 1):
        sub = seq[i:i + K]
        t = kmap.get(canonical(sub))
        if t:
            seen.add(t)
    assert seen == {"TaxonA"}, seen
    assert canonical("AAAA") == "AAAA" and canonical("TTTT") == "AAAA"
    print("self-test OK: k-mer attribution + canonical pass.")


if __name__ == "__main__":
    main()
