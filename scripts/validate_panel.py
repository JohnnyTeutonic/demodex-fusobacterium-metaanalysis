"""
E1c -- characterise the diagnostic k-mer panel (sensitivity + specificity).

This converts the closed-reference k-mer method from "approximate, untested" to
"approximate, characterised". For each panel taxon we:

  POSITIVE control (sensitivity): take HELD-OUT 16S references of the SAME taxon
  (accessions beyond those used to build the panel), slide simulated amplicon
  reads across them, and measure the fraction detected by that taxon's
  diagnostic k-mers under the same >=2-hit rule used by build_index.

  NEGATIVE control (specificity): do the same with CONGENERIC / CONFAMILIAL
  challenge organisms that must NOT be called as the target (e.g. Bartonella
  henselae vs B. quintana; a second Neisseriaceae vs Snodgrassella). High
  detection here would reveal genus-level rather than species-level resolution.

Output: results/panel_validation.tsv (taxon, sensitivity, specificity-by-challenge).

Usage:
  python validate_panel.py
  python validate_panel.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import os
from typing import Dict, List, Set

import sra_audit as A
from fetch_references import fetch_16s, kmers_of, KMERS_DIR, K
from build_index import revcomp

RESULTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results")
WIN, STEP, MIN_HITS = 250, 100, 2

# held-out positive search (we take refs[8:] as held-out, beyond the build set)
POSITIVE = {
    "Snodgrassella_alvi":           "Snodgrassella alvi[Organism] AND 16S ribosomal RNA[Title]",
    "Bacillus_oleronius":           "Bacillus oleronius[Organism] AND 16S ribosomal RNA[Title]",
    "Cutibacterium_kroppenstedtii": "(Cutibacterium kroppenstedtii[Organism] OR Corynebacterium kroppenstedtii[Organism]) AND 16S ribosomal RNA[All Fields]",
    "Bartonella_quintana":          "Bartonella quintana[Organism] AND 16S ribosomal RNA[Title]",
}
# congeneric / confamilial challengers that must NOT be called as the target
CHALLENGE = {
    "Snodgrassella_alvi":           ["Neisseria gonorrhoeae", "Kingella kingae"],
    "Bacillus_oleronius":           ["Bacillus subtilis", "Heyndrickxia coagulans"],
    "Cutibacterium_kroppenstedtii": ["Cutibacterium acnes", "Corynebacterium diphtheriae"],
    "Bartonella_quintana":          ["Bartonella henselae", "Brucella abortus"],
}


def load_diag(taxon: str) -> Set[str]:
    km: Set[str] = set()
    with open(os.path.join(KMERS_DIR, taxon + ".txt"), encoding="utf-8") as fh:
        for line in fh:
            k = line.strip()
            if k:
                km.add(k)
                km.add(revcomp(k))
    return km


def detect_fraction(seqs: List[str], diag: Set[str]) -> float:
    """Fraction of simulated WINdows carrying >= MIN_HITS diagnostic k-mers."""
    windows, hits = 0, 0
    for s in seqs:
        s = "".join(c for c in s.upper() if c in "ACGT")
        if len(s) < WIN:
            if len(s) >= K:
                cnt = sum(1 for i in range(len(s) - K + 1) if s[i:i+K] in diag)
                windows += 1
                hits += 1 if cnt >= MIN_HITS else 0
            continue
        for start in range(0, len(s) - WIN + 1, STEP):
            w = s[start:start + WIN]
            cnt = sum(1 for i in range(WIN - K + 1) if w[i:i+K] in diag)
            windows += 1
            hits += 1 if cnt >= MIN_HITS else 0
    return hits / windows if windows else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return

    os.makedirs(RESULTS_DIR, exist_ok=True)
    rows = []
    for taxon, term in POSITIVE.items():
        diag = load_diag(taxon)
        # held-out positives: refs 9..18 (beyond the 8 used to build)
        held = fetch_16s(term, 18)[8:]
        sens = detect_fraction(held, diag) if held else float("nan")
        chal_results = {}
        for org in CHALLENGE[taxon]:
            seqs = fetch_16s(f"{org}[Organism] AND 16S ribosomal RNA[Title]", 4)
            chal_results[org] = detect_fraction(seqs, diag) if seqs else float("nan")
        rows.append({"taxon": taxon, "n_diag_kmers": len(diag) // 2,
                     "sensitivity_heldout": round(sens, 3),
                     "challenges": "; ".join(f"{o}={v:.3f}" for o, v in chal_results.items())})
        print(f"  {taxon:<30} sens={sens:.3f}  " +
              "  ".join(f"{o.split()[0][:4]}.{o.split()[-1][:4]}={v:.3f}"
                       for o, v in chal_results.items()))

    cols = ["taxon", "n_diag_kmers", "sensitivity_heldout", "challenges"]
    out = os.path.join(RESULTS_DIR, "panel_validation.tsv")
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out}")


def _self_test() -> None:
    diag = {"A" * K, revcomp("A" * K), "C" * K, revcomp("C" * K)}
    # a sequence with two adjacent diagnostic k-mers -> detected
    s = "A" * (K + 1)  # contains 2 overlapping A-kmers
    assert detect_fraction([s], diag) == 1.0, detect_fraction([s], diag)
    # a non-homopolymer sequence has no diagnostic k-mer -> not detected
    assert detect_fraction(["ACGT" * 12], diag) == 0.0
    print("self-test OK: detect_fraction sensitivity/specificity logic passes.")


if __name__ == "__main__":
    main()
