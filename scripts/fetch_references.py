"""
E1a -- fetch 16S rRNA reference sequences for the Demodex-associated taxon panel
and build per-taxon diagnostic k-mer sets.

The panel (from the rosacea/Demodex literature):
  Snodgrassella alvi          -- core Demodex microbiota; falls post-doxycycline
  Bacillus oleronius          -- classic Demodex symbiont
  Cutibacterium kroppenstedtii-- enriched in papulopustular rosacea
  Bartonella quintana         -- detected from Demodex in ETR rosacea
  Cutibacterium acnes         -- dominant skin/mite taxon (context / normaliser)
  Staphylococcus epidermidis  -- dominant skin taxon (context / normaliser)

Diagnostic k-mers: for each taxon we keep the k-mers (canonical, k=31) that occur
in that taxon's 16S references but in NO other panel member. Reads carrying a
diagnostic k-mer are attributed to that taxon. Context taxa are included so their
k-mers are subtracted out (they must not inflate the Demodex-associated signal).

Outputs:
  data/refs/<taxon>.fasta        raw 16S references
  data/kmers/<taxon>.txt         diagnostic canonical k-mers (one per line)
  data/kmers/panel_stats.tsv     taxon, n_refs, n_kmers_raw, n_kmers_diagnostic

Usage:
  python fetch_references.py
  python fetch_references.py --self-test
"""

from __future__ import annotations

import argparse
import os
import xml.etree.ElementTree as ET
from typing import Dict, List, Set

import sra_audit as A  # _eutil, esearch, ensure_data_dir

DATA_DIR = A.DATA_DIR
REFS_DIR = os.path.join(DATA_DIR, "refs")
KMERS_DIR = os.path.join(DATA_DIR, "kmers")
K = 31

# taxon -> (search term, is_demodex_associated)
PANEL = {
    "Snodgrassella_alvi":          ("Snodgrassella alvi[Organism] AND 16S ribosomal RNA[Title]", True),
    "Bacillus_oleronius":          ("Bacillus oleronius[Organism] AND 16S ribosomal RNA[Title]", True),
    "Cutibacterium_kroppenstedtii":("(Cutibacterium kroppenstedtii[Organism] OR Corynebacterium kroppenstedtii[Organism]) AND 16S ribosomal RNA[All Fields]", True),
    "Bartonella_quintana":         ("Bartonella quintana[Organism] AND 16S ribosomal RNA[Title]", True),
    "Cutibacterium_acnes":         ("Cutibacterium acnes[Organism] AND 16S ribosomal RNA[Title]", False),
    "Staphylococcus_epidermidis":  ("Staphylococcus epidermidis[Organism] AND 16S ribosomal RNA[Title]", False),
}
MAX_REFS = 8  # per taxon; plenty for stable diagnostic k-mers

# Broad background of common human skin / facial bacteria. Diagnostic k-mers must
# be absent from ALL of these (i.e. they must not fall in conserved 16S regions
# shared with the dominant skin flora), otherwise rare Demodex-associated taxa
# get inflated by cross-mapping. These are background only -- never scored.
# --- skin / facial flora (the original "naive" background) ---
SKIN_BACKGROUND = [
    "Staphylococcus aureus", "Staphylococcus capitis", "Staphylococcus hominis",
    "Streptococcus mitis", "Streptococcus pyogenes", "Micrococcus luteus",
    "Corynebacterium tuberculostearicum", "Corynebacterium amycolatum",
    "Lawsonella clevelandensis", "Cutibacterium granulosum",
    "Pseudomonas aeruginosa", "Acinetobacter baumannii", "Moraxella catarrhalis",
    "Enhydrobacter aerosaccus", "Finegoldia magna", "Anaerococcus prevotii",
    "Veillonella parvula", "Prevotella intermedia", "Bacillus subtilis",
    "Bacillus cereus", "Escherichia coli", "Lactobacillus crispatus",
    "Malassezia restricta", "Neisseria flava",
]
# --- ocular-surface / environmental flora + Rhizobiales relatives of Bartonella
#     and Bacillus relatives, added to remove the body-site cross-mapping that
#     inflated Bartonella/Bacillus at the ocular surface ---
OCULAR_BACKGROUND = [
    "Methylobacterium radiotolerans", "Sphingomonas paucimobilis",
    "Brevundimonas diminuta", "Bradyrhizobium japonicum", "Rhizobium leguminosarum",
    "Paracoccus denitrificans", "Ochrobactrum anthropi", "Brucella melitensis",
    "Acinetobacter johnsonii", "Bacillus pumilus", "Bacillus licheniformis",
    "Ralstonia pickettii", "Stenotrophomonas maltophilia", "Delftia acidovorans",
]
BACKGROUND_TAXA = SKIN_BACKGROUND + OCULAR_BACKGROUND
MAX_BG_REFS = 3  # per background taxon

_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def canonical(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def kmers_of(seq: str, k: int = K) -> Set[str]:
    seq = "".join(c for c in seq.upper() if c in "ACGT")
    out: Set[str] = set()
    for i in range(len(seq) - k + 1):
        out.add(canonical(seq[i:i + k]))
    return out


def parse_fasta(text: str) -> List[str]:
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def fetch_16s(term: str, max_refs: int) -> List[str]:
    ids = A.esearch("nucleotide", term, retmax=max_refs)
    if not ids:
        return []
    raw = A._eutil("efetch.fcgi",
                   {"db": "nucleotide", "id": ",".join(ids),
                    "rettype": "fasta", "retmode": "text"})
    return parse_fasta(raw.decode("utf-8", "replace"))


def build_panel(kmers_dir: str = KMERS_DIR, background_taxa=None) -> None:
    background_taxa = BACKGROUND_TAXA if background_taxa is None else background_taxa
    os.makedirs(REFS_DIR, exist_ok=True)
    os.makedirs(kmers_dir, exist_ok=True)

    raw_kmers: Dict[str, Set[str]] = {}
    n_refs: Dict[str, int] = {}
    for taxon, (term, _assoc) in PANEL.items():
        seqs = fetch_16s(term, MAX_REFS)
        n_refs[taxon] = len(seqs)
        with open(os.path.join(REFS_DIR, taxon + ".fasta"), "w", encoding="utf-8") as fh:
            for i, s in enumerate(seqs):
                fh.write(f">{taxon}_{i}\n{s}\n")
        km: Set[str] = set()
        for s in seqs:
            km |= kmers_of(s)
        raw_kmers[taxon] = km
        print(f"  {taxon:<30} refs={len(seqs):<3} raw_kmers={len(km)}")

    # broad background of common flora; diagnostic k-mers must avoid it
    print(f"Building background k-mer set ({len(background_taxa)} taxa) ...")
    background: Set[str] = set()
    bg_refs = 0
    for org in background_taxa:
        try:
            seqs = fetch_16s(f"{org}[Organism] AND 16S ribosomal RNA[Title]", MAX_BG_REFS)
        except Exception as e:
            print(f"    {org}: fetch failed {e!r}")
            continue
        bg_refs += len(seqs)
        for s in seqs:
            background |= kmers_of(s)
    print(f"  background: {bg_refs} refs, {len(background)} k-mers")
    with open(os.path.join(kmers_dir, "_background.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(sorted(background)))

    # diagnostic = unique to this taxon across the panel AND absent from background
    stats = []
    for taxon, km in raw_kmers.items():
        others: Set[str] = set()
        for other, okm in raw_kmers.items():
            if other != taxon:
                others |= okm
        diag = km - others - background
        with open(os.path.join(kmers_dir, taxon + ".txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(sorted(diag)))
        assoc = PANEL[taxon][1]
        stats.append((taxon, n_refs[taxon], len(km), len(diag), assoc))
        print(f"  {taxon:<30} diagnostic_kmers={len(diag)}  demodex_assoc={assoc}")

    with open(os.path.join(kmers_dir, "panel_stats.tsv"), "w", encoding="utf-8") as fh:
        fh.write("taxon\tn_refs\tn_kmers_raw\tn_kmers_diagnostic\tdemodex_associated\n")
        for row in stats:
            fh.write("\t".join(str(x) for x in row) + "\n")


def _self_test() -> None:
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAAA") == "TTTT"
    assert canonical("TTTT") == "AAAA"
    km = kmers_of("ACGTACGTACGT", k=4)
    assert all(len(x) == 4 for x in km)
    fa = parse_fasta(">a\nACGT\nACGT\n>b\nTTTT\n")
    assert fa == ["ACGTACGT", "TTTT"], fa
    print("self-test OK: revcomp/canonical/kmers/fasta parsers pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--out-dir", default=KMERS_DIR,
                    help="where to write the k-mer panel (default: data/kmers)")
    ap.add_argument("--skin-only", action="store_true",
                    help="use the skin-only (naive) background, no body-site extension")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    bg = SKIN_BACKGROUND if args.skin_only else BACKGROUND_TAXA
    label = "skin-only (naive)" if args.skin_only else "body-site-matched"
    print(f"Fetching 16S references + building diagnostic k-mers [{label}] ...")
    build_panel(kmers_dir=args.out_dir, background_taxa=bg)


if __name__ == "__main__":
    main()
