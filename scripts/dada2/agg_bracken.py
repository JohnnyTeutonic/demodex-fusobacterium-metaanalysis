"""
Aggregate per-sample Bracken genus outputs (PRJNA856121 shotgun) into the
asv_table.tsv / taxonomy.tsv format consumed by cross_cohort_concordance.py.

- One "ASV" per bacterial genus (ASV id = G<taxid>); Genus column = NCBI
  scientific name (matches SILVA names for the key taxa: Fusobacterium,
  Corynebacterium, Streptococcus, Rothia, Neisseria, Pseudomonas, ...).
- Counts = Bracken new_est_reads.
- Restricted to genera whose NCBI lineage passes through Bacteria (taxid 2),
  so host (Homo) and viral reads are excluded and relative abundances are
  computed among bacteria only -- comparable to the 16S cohorts.

Pure stdlib.
Usage: python agg_bracken.py --work ~/rosacea_856 --db ~/k2db
Writes <work>/out/asv_table.tsv and <work>/out/taxonomy.tsv
"""
from __future__ import annotations

import argparse
import csv
import os
import glob


def load_taxonomy(db: str):
    """Return (parent, rank, name) dicts from NCBI nodes.dmp/names.dmp."""
    parent, rank, name = {}, {}, {}
    with open(os.path.join(db, "nodes.dmp"), encoding="utf-8") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            tid = int(f[0]); parent[tid] = int(f[1]); rank[tid] = f[2]
    with open(os.path.join(db, "names.dmp"), encoding="utf-8") as fh:
        for line in fh:
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 4 and f[3] == "scientific name":
                name[int(f[0])] = f[1]
    return parent, rank, name


def is_bacterial(tid: int, parent: dict, cache: dict) -> bool:
    seen = []
    t = tid
    while t not in (0, 1) and t in parent:
        if t in cache:
            res = cache[t]
            for s in seen:
                cache[s] = res
            return res
        if t == 2:  # Bacteria
            for s in seen + [t]:
                cache[s] = True
            return True
        seen.append(t)
        nt = parent[t]
        if nt == t:
            break
        t = nt
    for s in seen:
        cache[s] = False
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.path.expanduser("~/rosacea_856"))
    ap.add_argument("--db", default=os.path.expanduser("~/k2db"))
    args = ap.parse_args()

    parent, rank, name = load_taxonomy(args.db)
    cache: dict = {}

    files = sorted(glob.glob(os.path.join(args.work, "kraken", "bracken", "*.bracken")))
    if not files:
        raise SystemExit("no bracken outputs found")
    samples = [os.path.splitext(os.path.basename(f))[0] for f in files]

    # matrix[taxid][sample] = new_est_reads (bacterial genera only)
    matrix: dict = {}
    for f, s in zip(files, samples):
        with open(f, encoding="utf-8") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for r in rdr:
                tid = int(r["taxonomy_id"])
                if r.get("taxonomy_lvl") != "G":
                    continue
                if not is_bacterial(tid, parent, cache):
                    continue
                reads = float(r["new_est_reads"])
                matrix.setdefault(tid, {})[s] = matrix.get(tid, {}).get(s, 0.0) + reads

    out = os.path.join(args.work, "out"); os.makedirs(out, exist_ok=True)
    taxids = sorted(matrix.keys())

    with open(os.path.join(out, "asv_table.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ASV"] + samples)
        for tid in taxids:
            row = [f"G{tid}"] + [f"{matrix[tid].get(s, 0.0):.2f}" for s in samples]
            w.writerow(row)

    with open(os.path.join(out, "taxonomy.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["ASV", "Domain", "Phylum", "Class", "Order", "Family", "Genus"])
        for tid in taxids:
            w.writerow([f"G{tid}", "Bacteria", "", "", "", "", name.get(tid, f"taxid_{tid}")])

    print(f"samples: {len(samples)}; bacterial genera: {len(taxids)}")
    print(f"wrote {out}/asv_table.tsv and taxonomy.tsv")


if __name__ == "__main__":
    main()
