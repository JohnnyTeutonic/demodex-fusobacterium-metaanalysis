#!/usr/bin/env python3
"""Test Neisseriaceae / Snodgrassella relative abundance across the PRJDB18292
2x2 design (site: lesional vs non-lesional; time: pre vs post ivermectin).

Inputs (under $WORK/out, default ~/db18292/out, mirrored to repo if copied):
  genus_by_sample.tsv  (sample x genus counts)
  genus_family.tsv     (genus -> family)
  sample_sheet.tsv     (run, site, time, group)   [under $WORK]
"""
import os, sys
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

WORK = os.environ.get("WORK", os.path.expanduser("~/db18292"))
OUT = os.path.join(WORK, "out")
gbs = pd.read_csv(os.path.join(OUT, "genus_by_sample.tsv"), sep="\t", index_col=0)
fam = pd.read_csv(os.path.join(OUT, "genus_family.tsv"), sep="\t", index_col=0)["family"].to_dict()
ss = pd.read_csv(os.path.join(WORK, "sample_sheet.tsv"), sep="\t", index_col=0)

# relative abundance per sample
rel = gbs.div(gbs.sum(axis=1), axis=0)

# focus taxa
neisseria = [g for g in gbs.columns if g == "Neisseria"]
neisseriaceae = [g for g in gbs.columns if fam.get(g, "") == "Neisseriaceae"]
snod = [g for g in gbs.columns if g == "Snodgrassella"]

def agg(cols):
    return rel[cols].sum(axis=1) if cols else pd.Series(0.0, index=rel.index)

focus = pd.DataFrame({
    "Neisseria": agg(neisseria),
    "Neisseriaceae": agg(neisseriaceae),
    "Snodgrassella": agg(snod),
})
focus = focus.join(ss[["site", "time", "group"]])

print("== focus taxa relative abundance (per sample) ==")
print(focus.round(5).to_string())
print("\nNeisseriaceae member genera detected:", neisseriaceae)
print("Snodgrassella detected in any sample:", bool(snod) and float(rel[snod].sum().sum()) > 0)

def test(label, taxon, mask_a, name_a, mask_b, name_b):
    a = focus.loc[mask_a, taxon].values
    b = focus.loc[mask_b, taxon].values
    try:
        U, p = mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        p = float("nan")
    print(f"  [{taxon}] {name_a} (med={np.median(a):.5f}) vs {name_b} "
          f"(med={np.median(b):.5f})  Mann-Whitney p={p:.4f}")

print("\n== contrast 1: lesional vs non-lesional (12 v 12) ==")
for t in ["Neisseria", "Neisseriaceae", "Snodgrassella"]:
    test("site", t, focus.site == "lesional", "lesional",
         focus.site == "nonlesional", "nonlesional")

print("\n== contrast 2: pre vs post ivermectin (12 v 12) ==")
for t in ["Neisseria", "Neisseriaceae", "Snodgrassella"]:
    test("time", t, focus.time == "pre", "pre",
         focus.time == "post", "post")

print("\n== contrast 3: lesional pre vs post (6 v 6) ==")
for t in ["Neisseria", "Neisseriaceae", "Snodgrassella"]:
    test("lestime", t, (focus.site == "lesional") & (focus.time == "pre"), "les-pre",
         (focus.site == "lesional") & (focus.time == "post"), "les-post")

print("\n== top 15 genera overall (mean relative abundance) ==")
print((rel.mean().sort_values(ascending=False).head(15) * 100).round(2).to_string())
