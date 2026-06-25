"""
E1b -- fetch run + sample metadata for a BioProject from ENA.

ENA gives a run table (run_accession, sample_accession, fastq_ftp, read_count)
via the portal filereport API, and per-sample attributes (DESCRIPTION,
host_subject_id, chem_administration, disease, ...) via the browser XML API.
We merge them into one tidy table the index/analysis steps consume.

Outputs:
  data/meta_<PROJECT>.tsv

Usage:
  python fetch_metadata.py --project PRJDB18292
  python fetch_metadata.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import os
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Dict, List

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml"

# attribute tags we surface as columns (lower-cased match); everything else is
# still preserved in the raw dict but these get dedicated columns.
KEY_ATTRS = ["host_subject_id", "chem_administration", "dermatology_disord",
             "disease", "disease state", "collection_date", "env_medium",
             "body site", "tissue"]


def _get(url: str, timeout: float = 90.0) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "rosacea-demodex/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def run_table(project: str) -> List[Dict[str, str]]:
    fields = "run_accession,sample_accession,sample_title,library_strategy,read_count,fastq_ftp,fastq_bytes"
    url = f"{ENA_PORTAL}?{urllib.parse.urlencode({'accession': project, 'result': 'read_run', 'fields': fields, 'format': 'tsv'})}"
    text = _get(url).decode("utf-8", "replace")
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    return rows


def parse_sample_xml(xml_bytes: bytes) -> Dict[str, str]:
    root = ET.fromstring(xml_bytes)
    out: Dict[str, str] = {}
    s = root.find(".//SAMPLE")
    if s is None:
        return out
    out["description"] = (s.findtext("DESCRIPTION", default="") or "").strip()
    out["scientific_name"] = (s.findtext(".//SCIENTIFIC_NAME", default="") or "").strip()
    for a in s.findall(".//SAMPLE_ATTRIBUTE"):
        tag = (a.findtext("TAG", default="") or "").strip().lower()
        val = (a.findtext("VALUE", default="") or "").strip()
        if tag:
            out[tag] = val
    return out


def classify_timepoint(description: str) -> str:
    d = description.lower()
    if "pre-treatment" in d or "pretreatment" in d or "baseline" in d or "before" in d:
        return "pre"
    if "post-treatment" in d or "posttreatment" in d or "after" in d or "week" in d:
        return "post"
    return ""


def classify_lesional(description: str) -> str:
    d = description.lower()
    if "non-lesional" in d or "nonlesional" in d:
        return "non-lesional"
    if "lesional" in d:
        return "lesional"
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", type=str)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        _self_test()
        return
    if not args.project:
        ap.error("--project required")

    os.makedirs(DATA_DIR, exist_ok=True)
    runs = run_table(args.project)
    print(f"{args.project}: {len(runs)} runs")

    merged: List[Dict[str, str]] = []
    sample_cache: Dict[str, Dict[str, str]] = {}
    for r in runs:
        sacc = r.get("sample_accession", "")
        if sacc and sacc not in sample_cache:
            try:
                sample_cache[sacc] = parse_sample_xml(_get(f"{ENA_XML}/{sacc}"))
            except Exception as e:
                sample_cache[sacc] = {"description": f"(fetch failed: {e!r})"}
        attrs = sample_cache.get(sacc, {})
        desc = attrs.get("description", "")
        row = {
            "run_accession": r.get("run_accession", ""),
            "sample_accession": sacc,
            "sample_title": r.get("sample_title", ""),
            "read_count": r.get("read_count", ""),
            "fastq_ftp": (r.get("fastq_ftp", "") or "").split(";")[0],
            "description": desc,
            "timepoint": classify_timepoint(desc),
            "lesional": classify_lesional(desc),
        }
        for k in KEY_ATTRS:
            row[k] = attrs.get(k, "")
        merged.append(row)
        print(f"  {row['run_accession']:<12} subj={row.get('host_subject_id',''):<3} "
              f"{row['timepoint']:<4} {row['lesional']:<12} {desc[:50]}")

    cols = ["run_accession", "sample_accession", "sample_title", "read_count",
            "timepoint", "lesional", "host_subject_id", "chem_administration",
            "dermatology_disord", "disease", "disease state", "description",
            "fastq_ftp"]
    out = os.path.join(DATA_DIR, f"meta_{args.project}.tsv")
    with open(out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(merged)
    print(f"\nWrote {len(merged)} rows -> {out}")


def _self_test() -> None:
    assert classify_timepoint("lesional skin pre-treatment") == "pre"
    assert classify_timepoint("Skin swab post-treatment week 12") == "post"
    assert classify_lesional("non-lesional skin pre-treatment") == "non-lesional"
    assert classify_lesional("lesional skin") == "lesional"
    x = parse_sample_xml(_SAMPLE_XML)
    assert x["description"].startswith("Skin swab"), x
    assert x["host_subject_id"] == "1", x
    assert x["chem_administration"] == "Ivermectin", x
    print("self-test OK: timepoint/lesional/sample-xml parsers pass.")


_SAMPLE_XML = b"""<SAMPLE_SET><SAMPLE accession="X"><TITLE>t</TITLE>
<DESCRIPTION>Skin swab sample of lesional skin pre-treatment</DESCRIPTION>
<SAMPLE_NAME><SCIENTIFIC_NAME>metagenome</SCIENTIFIC_NAME></SAMPLE_NAME>
<SAMPLE_ATTRIBUTES>
<SAMPLE_ATTRIBUTE><TAG>host_subject_id</TAG><VALUE>1</VALUE></SAMPLE_ATTRIBUTE>
<SAMPLE_ATTRIBUTE><TAG>chem_administration</TAG><VALUE>Ivermectin</VALUE></SAMPLE_ATTRIBUTE>
</SAMPLE_ATTRIBUTES></SAMPLE></SAMPLE_SET>"""


if __name__ == "__main__":
    main()
