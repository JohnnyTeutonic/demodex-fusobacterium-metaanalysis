"""
E1c -- probe a BioProject's public ENA records for usable group labels.

The facial cohort PRJEB82826 (JID 2025 multi-omics) has FASTQs on ENA, but the
labelled validation is "blocked on the authors' sample key": our merged metadata
table only surfaces a fixed KEY_ATTRS whitelist and may be dropping a label tag.

This script makes no assumptions: it fetches the full ENA run table (extended
fields) plus every sample's complete SAMPLE_ATTRIBUTE set, then
  (a) inventories *every* attribute tag seen (with coverage + example values), and
  (b) scans all free-text fields and attribute values for label-like tokens
      (demodex, mite, infestation, control, healthy, lesional, density, ...).

If a deposited label exists anywhere in the public record, it surfaces here.
If nothing surfaces, that is itself the answer: the key is not public.

Outputs:
  data/probe_<PROJECT>_attributes.tsv   (tag, n_samples, n_nonempty, examples)
  data/probe_<PROJECT>_runs.tsv         (full extended run table + label hits)

Usage:
  python probe_labels.py --project PRJEB82826
  python probe_labels.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from typing import Dict, List

DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir, "data")
ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_XML = "https://www.ebi.ac.uk/ena/browser/api/xml"

# Extended free-text run fields that sometimes carry the group label directly.
RUN_FIELDS = [
    "run_accession", "sample_accession", "experiment_accession", "study_accession",
    "sample_title", "sample_alias", "sample_description", "description",
    "library_name", "experiment_title", "study_title", "run_alias",
    "read_count", "fastq_ftp",
]

# Tokens that would indicate a group/clinical label is present.
LABEL_TOKENS = [
    "demodex", "mite", "infest", "infestation", "densit", "count",
    "control", "healthy", "case", "patient", "rosacea", "lesion",
    "papulopust", "erythematotelangiectatic", "ett", "ppr",
    "positive", "negative", "pos", "neg", "affected", "unaffected", "group",
]
TOKEN_RE = re.compile("|".join(re.escape(t) for t in LABEL_TOKENS), re.IGNORECASE)


def _get(url: str, timeout: float = 90.0) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "rosacea-demodex/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def run_table(project: str) -> List[Dict[str, str]]:
    params = {
        "accession": project, "result": "read_run",
        "fields": ",".join(RUN_FIELDS), "format": "tsv",
    }
    url = f"{ENA_PORTAL}?{urllib.parse.urlencode(params)}"
    text = _get(url).decode("utf-8", "replace")
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def parse_sample_attrs(xml_bytes: bytes) -> Dict[str, str]:
    """Return ALL sample attributes (lower-cased tags) plus description/alias/title."""
    root = ET.fromstring(xml_bytes)
    out: Dict[str, str] = {}
    s = root.find(".//SAMPLE")
    if s is None:
        return out
    out["_title"] = (s.findtext("TITLE", default="") or "").strip()
    out["_description"] = (s.findtext("DESCRIPTION", default="") or "").strip()
    out["_alias"] = (s.get("alias", "") or "").strip()
    out["_scientific_name"] = (s.findtext(".//SCIENTIFIC_NAME", default="") or "").strip()
    for a in s.findall(".//SAMPLE_ATTRIBUTE"):
        tag = (a.findtext("TAG", default="") or "").strip().lower()
        val = (a.findtext("VALUE", default="") or "").strip()
        if tag:
            out[tag] = val
    return out


def label_hits(*texts: str) -> str:
    """Distinct label tokens found across the given free-text fields."""
    found = []
    for t in texts:
        for m in TOKEN_RE.findall(t or ""):
            ml = m.lower()
            if ml not in found:
                found.append(ml)
    return ";".join(found)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", type=str)
    ap.add_argument("--max-samples", type=int, default=0,
                    help="limit distinct samples fetched (0 = all)")
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

    sample_cache: Dict[str, Dict[str, str]] = {}
    tag_nonempty: Counter = Counter()
    tag_examples: Dict[str, List[str]] = defaultdict(list)
    n_samples = 0
    for r in runs:
        sacc = r.get("sample_accession", "")
        if not sacc or sacc in sample_cache:
            continue
        if args.max_samples and n_samples >= args.max_samples:
            break
        try:
            attrs = parse_sample_attrs(_get(f"{ENA_XML}/{sacc}"))
        except Exception as e:
            attrs = {"_error": repr(e)}
        sample_cache[sacc] = attrs
        n_samples += 1
        for tag, val in attrs.items():
            if val:
                tag_nonempty[tag] += 1
                if len(tag_examples[tag]) < 5 and val not in tag_examples[tag]:
                    tag_examples[tag].append(val)

    # (a) attribute inventory
    attr_out = os.path.join(DATA_DIR, f"probe_{args.project}_attributes.tsv")
    with open(attr_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tag", "n_samples_total", "n_nonempty", "is_label_like", "examples"])
        for tag, _ in tag_nonempty.most_common():
            is_label = bool(TOKEN_RE.search(tag) or any(TOKEN_RE.search(v) for v in tag_examples[tag]))
            w.writerow([tag, n_samples, tag_nonempty[tag], "YES" if is_label else "",
                        " | ".join(tag_examples[tag])])
    print(f"\nAttribute inventory ({len(tag_nonempty)} tags) -> {attr_out}")
    print("  label-like tags:")
    any_label = False
    for tag in tag_nonempty:
        if TOKEN_RE.search(tag) or any(TOKEN_RE.search(v) for v in tag_examples[tag]):
            any_label = True
            print(f"    {tag:<28} e.g. {tag_examples[tag][:3]}")
    if not any_label:
        print("    (none) -- no attribute tag or value carries a group/clinical token.")

    # (b) per-run free-text label scan
    run_out = os.path.join(DATA_DIR, f"probe_{args.project}_runs.tsv")
    cols = RUN_FIELDS + ["sample_title_x", "sample_desc_x", "sample_alias_x", "label_hits"]
    n_hits = 0
    with open(run_out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in runs:
            attrs = sample_cache.get(r.get("sample_accession", ""), {})
            stitle, sdesc, salias = attrs.get("_title", ""), attrs.get("_description", ""), attrs.get("_alias", "")
            hits = label_hits(r.get("sample_title", ""), r.get("sample_alias", ""),
                              r.get("sample_description", ""), r.get("description", ""),
                              r.get("library_name", ""), r.get("experiment_title", ""),
                              stitle, sdesc, salias, *[v for v in attrs.values()])
            if hits:
                n_hits += 1
            row = dict(r)
            row.update({"sample_title_x": stitle, "sample_desc_x": sdesc,
                        "sample_alias_x": salias, "label_hits": hits})
            w.writerow(row)
    print(f"\nRun table + label scan -> {run_out}")
    print(f"  runs with any label token: {n_hits}/{len(runs)}")
    print("\nVERDICT:", "labels may be recoverable -- inspect label-like tags/hits above."
          if any_label or n_hits else
          "no public labels found -- sample key must be obtained from the authors.")


def _self_test() -> None:
    attrs = parse_sample_attrs(_SAMPLE_XML)
    assert attrs["_description"].startswith("Facial follicle"), attrs
    assert attrs["host_disease_status"] == "Demodex-positive rosacea", attrs
    assert attrs["mite_density"] == "Dd>5/cm2", attrs
    h = label_hits(attrs["_description"], attrs["host_disease_status"], attrs["mite_density"])
    for tok in ("demodex", "rosacea", "positive"):
        assert tok in h, (tok, h)
    # tag-name detection (handled in main) catches label-like tags by their name
    assert TOKEN_RE.search("mite_density"), "tag-name scan should flag mite_density"
    assert label_hits("Skin swab, replicate 2") == "", "should find no label tokens"
    print("self-test OK: full-attr parse + label-token scan pass.")


_SAMPLE_XML = b"""<SAMPLE_SET><SAMPLE alias="P07-cheek"><TITLE>P07</TITLE>
<DESCRIPTION>Facial follicle glue-biopsy, right cheek</DESCRIPTION>
<SAMPLE_NAME><SCIENTIFIC_NAME>human skin metagenome</SCIENTIFIC_NAME></SAMPLE_NAME>
<SAMPLE_ATTRIBUTES>
<SAMPLE_ATTRIBUTE><TAG>host_disease_status</TAG><VALUE>Demodex-positive rosacea</VALUE></SAMPLE_ATTRIBUTE>
<SAMPLE_ATTRIBUTE><TAG>mite_density</TAG><VALUE>Dd>5/cm2</VALUE></SAMPLE_ATTRIBUTE>
</SAMPLE_ATTRIBUTES></SAMPLE></SAMPLE_SET>"""


if __name__ == "__main__":
    main()
