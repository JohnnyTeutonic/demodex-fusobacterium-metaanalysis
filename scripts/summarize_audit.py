"""
E0b -- turn the raw SRA run audit into a ranked, human-readable project table.

NCBI's BioProject esummary endpoint is currently throwing a server-side DTD
error (verified 2026-06), so we recover study-level metadata from the SRA side:
for each distinct BioProject we efetch one full SRA record, which carries
STUDY_TITLE / STUDY_ABSTRACT and SAMPLE attributes (host, disease, etc.).

Outputs results/project_summary.tsv and prints a ranked view, flagging which
projects mention rosacea / Demodex and what marker / strategy they used.

Usage:
  python summarize_audit.py
  python summarize_audit.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import os
import collections
import xml.etree.ElementTree as ET
from typing import Dict, List

import sra_audit as A  # reuse _eutil, esearch, ensure_data_dir

DATA_DIR = A.DATA_DIR
RESULTS_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "results")


def efetch_sra_full(uid: str) -> bytes:
    return A._eutil("efetch.fcgi",
                    {"db": "sra", "id": uid, "rettype": "full", "retmode": "xml"})


def parse_study(xml_bytes: bytes) -> Dict[str, str]:
    """Pull study title/abstract + a sampling of sample attributes."""
    root = ET.fromstring(xml_bytes)
    title = (root.findtext(".//STUDY/DESCRIPTOR/STUDY_TITLE", default="") or "").strip()
    abstract = (root.findtext(".//STUDY/DESCRIPTOR/STUDY_ABSTRACT", default="") or "").strip()
    org = ""
    org_el = root.find(".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME")
    if org_el is not None and org_el.text:
        org = org_el.text.strip()
    # sample attribute keys give a hint of disease/condition design
    attrs = []
    for a in root.findall(".//SAMPLE_ATTRIBUTE"):
        tag = a.findtext("TAG", default="")
        if tag:
            attrs.append(tag.lower())
    return {"title": title, "abstract": abstract[:400], "organism": org,
            "sample_attrs": ";".join(sorted(set(attrs))[:12])}


def first_uid_for_bioproject(acc: str) -> str:
    ids = A.esearch("sra", f"{acc}[BioProject]", retmax=1)
    return ids[0] if ids else ""


def relevance(title: str, abstract: str) -> str:
    t = (title + " " + abstract).lower()
    tags = []
    if "rosacea" in t:
        tags.append("ROSACEA")
    if "demodex" in t:
        tags.append("DEMODEX")
    if "doxycyc" in t or "antibiotic" in t or "ivermectin" in t or "treatment" in t:
        tags.append("TREATMENT")
    if "16s" in t or "amplicon" in t:
        tags.append("16S")
    return ",".join(tags) if tags else "-"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return

    sra_path = os.path.join(DATA_DIR, "audit_sra.tsv")
    rows = list(csv.DictReader(open(sra_path, encoding="utf-8"), delimiter="\t"))

    by_bp: Dict[str, List[dict]] = collections.defaultdict(list)
    for r in rows:
        by_bp[r.get("bioproject", "")].append(r)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    out_rows = []
    print(f"Resolving {len(by_bp)} BioProjects ...")
    for acc, runs in sorted(by_bp.items(), key=lambda kv: -len(kv[1])):
        if not acc:
            continue
        strat = collections.Counter(x.get("library_strategy", "") for x in runs)
        spots = [int(x["spots"]) for x in runs if x.get("spots", "").isdigit()]
        meta = {"title": "", "abstract": "", "organism": "", "sample_attrs": ""}
        try:
            uid = first_uid_for_bioproject(acc)
            if uid:
                meta = parse_study(efetch_sra_full(uid))
        except Exception as e:
            meta["title"] = f"(fetch failed: {e!r})"
        tag = relevance(meta["title"], meta["abstract"])
        out_rows.append({
            "bioproject": acc,
            "n_runs": len(runs),
            "strategies": ";".join(f"{k}:{v}" for k, v in strat.most_common()),
            "median_spots": sorted(spots)[len(spots) // 2] if spots else "",
            "relevance": tag,
            "organism": meta["organism"],
            "title": meta["title"],
            "sample_attrs": meta["sample_attrs"],
        })
        print(f"  {acc:<14} n={len(runs):<4} {tag:<28} {meta['title'][:70]}")

    cols = ["bioproject", "n_runs", "strategies", "median_spots", "relevance",
            "organism", "title", "sample_attrs"]
    out_path = os.path.join(RESULTS_DIR, "project_summary.tsv")
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    hits = [r for r in out_rows if r["relevance"] not in ("-", "16S")]
    print(f"\nWrote {len(out_rows)} projects -> {out_path}")
    print(f"On-topic (rosacea/Demodex/treatment) projects: {len(hits)}")


_STUDY_SAMPLE = b"""<EXPERIMENT_PACKAGE_SET><EXPERIMENT_PACKAGE>
<STUDY><DESCRIPTOR><STUDY_TITLE>Skin microbiome in rosacea before and after doxycycline</STUDY_TITLE>
<STUDY_ABSTRACT>16S rRNA amplicon sequencing of facial skin in rosacea patients.</STUDY_ABSTRACT>
</DESCRIPTOR></STUDY>
<SAMPLE><SAMPLE_NAME><SCIENTIFIC_NAME>human skin metagenome</SCIENTIFIC_NAME></SAMPLE_NAME>
<SAMPLE_ATTRIBUTES><SAMPLE_ATTRIBUTE><TAG>disease state</TAG><VALUE>rosacea</VALUE></SAMPLE_ATTRIBUTE>
<SAMPLE_ATTRIBUTE><TAG>body site</TAG><VALUE>cheek</VALUE></SAMPLE_ATTRIBUTE></SAMPLE_ATTRIBUTES>
</SAMPLE></EXPERIMENT_PACKAGE></EXPERIMENT_PACKAGE_SET>"""


def _self_test() -> None:
    s = parse_study(_STUDY_SAMPLE)
    assert "rosacea" in s["title"].lower(), s
    assert s["organism"] == "human skin metagenome", s
    assert "disease state" in s["sample_attrs"], s
    assert relevance(s["title"], s["abstract"]) and "ROSACEA" in relevance(s["title"], s["abstract"])
    print("self-test OK: study parser + relevance tagging pass.")


if __name__ == "__main__":
    main()
