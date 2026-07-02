"""
WS2 -- build a Neisseriaceae-resolving diagnostic k-mer panel to identify the
Demodex-associated bacterium (Snodgrassella vs Neisseria vs Eikenella).

Reuses fetch_references.build_panel but with an EXTENDED panel in which the three
candidate Neisseriaceae genera are all scored members, so their mutually shared
(conserved) 16S k-mers cancel out and each retained diagnostic set is
genus-specific. Written to a SEPARATE dir (data/kmers_demodex) so the paper's
original panel in data/kmers is untouched.

Key correction: the default background contains "Neisseria flava"; we drop it
here because Neisseria is now a scored panel member (otherwise its diagnostic
k-mers would be subtracted away).

Usage:  python build_neisseriaceae_panel.py
"""
from __future__ import annotations
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
import fetch_references as F

OUT_DIR = os.path.join(F.DATA_DIR, "kmers_demodex")

# Extended panel: original taxa + Neisseria (genus) + Eikenella corrodens.
F.PANEL = {
    "Snodgrassella_alvi":          ("Snodgrassella alvi[Organism] AND 16S ribosomal RNA[Title]", True),
    "Neisseria":                   ("Neisseria[Organism] AND 16S ribosomal RNA[Title]", True),
    "Eikenella_corrodens":         ("Eikenella corrodens[Organism] AND 16S ribosomal RNA[Title]", True),
    "Bacillus_oleronius":          ("Bacillus oleronius[Organism] AND 16S ribosomal RNA[Title]", True),
    "Bartonella_quintana":         ("Bartonella quintana[Organism] AND 16S ribosomal RNA[Title]", True),
    "Cutibacterium_kroppenstedtii":("(Cutibacterium kroppenstedtii[Organism] OR Corynebacterium kroppenstedtii[Organism]) AND 16S ribosomal RNA[All Fields]", False),
    "Cutibacterium_acnes":         ("Cutibacterium acnes[Organism] AND 16S ribosomal RNA[Title]", False),
    "Staphylococcus_epidermidis":  ("Staphylococcus epidermidis[Organism] AND 16S ribosomal RNA[Title]", False),
}

# Background = default minus Neisseria flava (now a scored member).
bg = [x for x in F.BACKGROUND_TAXA if x != "Neisseria flava"]

if __name__ == "__main__":
    print(f"Building Neisseriaceae-resolving panel -> {OUT_DIR}")
    print(f"background taxa: {len(bg)} (Neisseria flava removed)")
    F.build_panel(kmers_dir=OUT_DIR, background_taxa=bg)
    print("PANEL_BUILD_OK")
