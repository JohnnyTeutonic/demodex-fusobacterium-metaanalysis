"""
E0 -- data-availability audit for the rosacea / Demodex bioinformatic project.

Hits NCBI E-utilities (pure standard library, no Biopython) to find which
rosacea / Demodex studies actually deposited reusable data, what marker they
used (16S / 18S / shotgun), the platform, and rough sample size.

Outputs:
  data/audit_bioproject.tsv   project_uid, accession, title, n_samples
  data/audit_sra.tsv          run-level: srr, experiment_title, platform,
                              library_strategy, library_source, bioproject, spots

Usage:
  python sra_audit.py                 # full audit, writes both TSVs
  python sra_audit.py --self-test     # offline: validate the XML parsers

Polite use: <=3 req/s without an API key; we sleep 0.34s between calls.
Set NCBI_API_KEY in the environment to lift to 10 req/s.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
API_KEY = os.environ.get("NCBI_API_KEY", "")
SLEEP = 0.11 if API_KEY else 0.34

# Search terms. Broad enough to catch the known studies, narrow enough to stay
# on-topic. We OR rosacea with Demodex and require a sequencing-ish context.
SRA_TERMS = [
    'rosacea[All Fields] AND microbiome[All Fields]',
    'rosacea[All Fields] AND (16S[All Fields] OR amplicon[All Fields])',
    'Demodex[Organism] OR Demodex[All Fields]',
    'rosacea[All Fields] AND skin[All Fields] AND bacteria[All Fields]',
]
BIOPROJECT_TERMS = [
    'rosacea[Title] OR rosacea[Description]',
    'Demodex[Title] OR Demodex[Description]',
]


def ensure_data_dir() -> str:
    os.makedirs(DATA_DIR, exist_ok=True)
    return DATA_DIR


def _get(url: str, timeout: float = 60.0) -> bytes:
    req = urllib.request.Request(
        url, headers={"User-Agent": "rosacea-demodex-audit/1.0 (research)"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def _eutil(tool: str, params: Dict[str, str]) -> bytes:
    if API_KEY:
        params = {**params, "api_key": API_KEY}
    url = f"{EUTILS}/{tool}?{urllib.parse.urlencode(params)}"
    time.sleep(SLEEP)
    return _get(url)


def esearch(db: str, term: str, retmax: int = 200) -> List[str]:
    """Return a list of UIDs for a search term."""
    raw = _eutil("esearch.fcgi",
                 {"db": db, "term": term, "retmax": str(retmax), "retmode": "xml"})
    root = ET.fromstring(raw)
    return [e.text for e in root.findall(".//IdList/Id") if e.text]


# ---------------------------------------------------------------------------
# Parsers (kept pure so they can be unit-tested offline)
# ---------------------------------------------------------------------------

def parse_bioproject_summary(xml_bytes: bytes) -> List[Dict[str, str]]:
    """Parse esummary (db=bioproject) docsums into rows."""
    root = ET.fromstring(xml_bytes)
    out: List[Dict[str, str]] = []
    for ds in root.findall(".//DocumentSummary"):
        uid = ds.attrib.get("uid", "")
        acc = ds.findtext("Project_Acc", default="")
        title = (ds.findtext("Project_Title", default="")
                 or ds.findtext("Project_Name", default="")).strip()
        out.append({"uid": uid, "accession": acc, "title": title})
    return out


def parse_sra_summary(xml_bytes: bytes) -> List[Dict[str, str]]:
    """Parse esummary (db=sra) docsums. The useful fields live in the
    ExpXml / Runs attribute blobs as escaped XML fragments."""
    root = ET.fromstring(xml_bytes)
    rows: List[Dict[str, str]] = []
    for ds in root.findall(".//DocSum"):
        fields = {it.attrib.get("Name", ""): (it.text or "")
                  for it in ds.findall("Item")}
        exp = fields.get("ExpXml", "")
        runs = fields.get("Runs", "")
        rows.append(_parse_exp_blob(exp, runs))
    # Newer esummary uses DocumentSummary with raw child tags too.
    for ds in root.findall(".//DocumentSummary"):
        exp = ds.findtext("ExpXml", default="")
        runs = ds.findtext("Runs", default="")
        if exp or runs:
            rows.append(_parse_exp_blob(exp, runs))
    return rows


def _wrap(fragment: str) -> ET.Element:
    """ExpXml/Runs are XML fragments without a single root; wrap them."""
    return ET.fromstring(f"<root>{fragment}</root>")


def _parse_exp_blob(exp: str, runs: str) -> Dict[str, str]:
    row = {"srr": "", "experiment_title": "", "platform": "",
           "library_strategy": "", "library_source": "", "bioproject": "",
           "organism": "", "spots": ""}
    if exp:
        try:
            e = _wrap(exp)
            row["experiment_title"] = (e.findtext(".//Title", default="") or "").strip()
            plat = e.find(".//Platform")
            if plat is not None:
                row["platform"] = (plat.text or plat.attrib.get("instrument_model", "")).strip()
            lib = e.find(".//Library_descriptor")
            if lib is not None:
                row["library_strategy"] = lib.findtext("LIBRARY_STRATEGY", default="")
                row["library_source"] = lib.findtext("LIBRARY_SOURCE", default="")
            bp = e.find(".//Bioproject")
            if bp is not None and bp.text:
                row["bioproject"] = bp.text.strip()
            org = e.find(".//Organism")
            if org is not None:
                row["organism"] = org.attrib.get("ScientificName", "")
        except ET.ParseError:
            pass
    if runs:
        try:
            r = _wrap(runs)
            run0 = r.find(".//Run")
            if run0 is not None:
                row["srr"] = run0.attrib.get("acc", "")
                row["spots"] = run0.attrib.get("total_spots", "")
        except ET.ParseError:
            pass
    return row


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def audit_bioprojects() -> List[Dict[str, str]]:
    seen: Dict[str, Dict[str, str]] = {}
    for term in BIOPROJECT_TERMS:
        ids = esearch("bioproject", term, retmax=100)
        print(f"  bioproject term '{term[:40]}...' -> {len(ids)} hits")
        for i in range(0, len(ids), 100):
            chunk = ids[i:i + 100]
            raw = _eutil("esummary.fcgi",
                         {"db": "bioproject", "id": ",".join(chunk), "retmode": "xml"})
            for row in parse_bioproject_summary(raw):
                if row["uid"]:
                    seen[row["uid"]] = row
    return list(seen.values())


def audit_sra() -> List[Dict[str, str]]:
    seen: Dict[str, Dict[str, str]] = {}
    for term in SRA_TERMS:
        ids = esearch("sra", term, retmax=300)
        print(f"  sra term '{term[:40]}...' -> {len(ids)} hits")
        for i in range(0, len(ids), 100):
            chunk = ids[i:i + 100]
            raw = _eutil("esummary.fcgi",
                         {"db": "sra", "id": ",".join(chunk), "retmode": "xml"})
            for row in parse_sra_summary(raw):
                key = row.get("srr") or row.get("experiment_title")
                if key:
                    seen[key] = row
    return list(seen.values())


def write_tsv(rows: List[Dict[str, str]], path: str, cols: List[str]) -> int:
    ensure_data_dir()
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    return len(rows)


# ---------------------------------------------------------------------------
# Offline self-test
# ---------------------------------------------------------------------------

_BP_SAMPLE = b"""<eSummaryResult><DocumentSummary uid="123456">
<Project_Acc>PRJNA000001</Project_Acc>
<Project_Title>Skin microbiome in rosacea patients</Project_Title>
</DocumentSummary></eSummaryResult>"""

_SRA_SAMPLE = (
    b"<eSummaryResult><DocSum><Id>1</Id>"
    b"<Item Name=\"ExpXml\" Type=\"String\">"
    b"&lt;Summary&gt;&lt;Title&gt;16S rRNA V3V4 rosacea cheek&lt;/Title&gt;&lt;/Summary&gt;"
    b"&lt;Platform instrument_model=\"Illumina MiSeq\"&gt;ILLUMINA&lt;/Platform&gt;"
    b"&lt;Library_descriptor&gt;&lt;LIBRARY_STRATEGY&gt;AMPLICON&lt;/LIBRARY_STRATEGY&gt;"
    b"&lt;LIBRARY_SOURCE&gt;METAGENOMIC&lt;/LIBRARY_SOURCE&gt;&lt;/Library_descriptor&gt;"
    b"&lt;Bioproject&gt;PRJNA000001&lt;/Bioproject&gt;"
    b"</Item>"
    b"<Item Name=\"Runs\" Type=\"String\">"
    b"&lt;Run acc=\"SRR0000001\" total_spots=\"123456\"/&gt;"
    b"</Item>"
    b"</DocSum></eSummaryResult>"
)


def _self_test() -> None:
    bp = parse_bioproject_summary(_BP_SAMPLE)
    assert len(bp) == 1 and bp[0]["accession"] == "PRJNA000001", bp
    assert "rosacea" in bp[0]["title"].lower(), bp
    sra = parse_sra_summary(_SRA_SAMPLE)
    assert len(sra) == 1, sra
    r = sra[0]
    assert r["srr"] == "SRR0000001", r
    assert r["library_strategy"] == "AMPLICON", r
    assert r["spots"] == "123456", r
    assert "MiSeq" in r["platform"] or r["platform"] == "ILLUMINA", r
    print("self-test OK: bioproject + sra esummary parsers pass.")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return

    print("Auditing BioProjects ...")
    bp = audit_bioprojects()
    n_bp = write_tsv(bp, os.path.join(DATA_DIR, "audit_bioproject.tsv"),
                     ["uid", "accession", "title"])
    print(f"  wrote {n_bp} bioproject rows")

    print("Auditing SRA ...")
    sra = audit_sra()
    n_sra = write_tsv(sra, os.path.join(DATA_DIR, "audit_sra.tsv"),
                      ["srr", "bioproject", "organism", "platform",
                       "library_strategy", "library_source", "spots",
                       "experiment_title"])
    print(f"  wrote {n_sra} sra rows")

    # quick marker breakdown
    from collections import Counter
    strat = Counter(r.get("library_strategy", "") for r in sra)
    print("\nLibrary-strategy breakdown:")
    for k, v in strat.most_common():
        print(f"  {k or '(blank)':<16} {v}")


if __name__ == "__main__":
    main()
