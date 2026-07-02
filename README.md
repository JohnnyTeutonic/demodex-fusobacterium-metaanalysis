# demodex-fusobacterium-metaanalysis

A cross-cohort meta-analysis of the ocular-surface bacterial microbiome in
*Demodex* disease. Three independent public cohorts are uniformly reprocessed
from raw reads and combined to ask which genus-level associations *reproduce*
across cohorts and across sequencing platforms, rather than analysing one new
dataset in isolation.

**Headline:** *Fusobacterium* is depleted in *Demodex* disease concordantly in
all three cohorts (per-cohort ROC-AUC 0.33 / 0.31 / 0.36; sqrt(*n*)-weighted
Stouffer combined *Z* = -2.93, *p* = 0.0034, Benjamini-Hochberg *q* = 0.088). It
is the top-ranked genus overall and the only one that stays significant
(*p* < 0.05) in every leave-one-cohort-out, and in the shotgun cohort the
depletion is independent of host-DNA fraction and sequencing depth.
*Corynebacterium* is concordantly (but weakly) enriched. Genera prominent in
individual cohorts (*Neisseria*, *Rothia*, *Streptococcus*) do **not** replicate,
reversing direction between studies. Manuscript:
`fusobacterium_demodex_metaanalysis.tex`.

## Cohorts (111 subjects: 69 Demodex, 42 control)
| label | accession | platform / target | Demodex | control |
|---|---|---|---|---|
| discovery   | PRJNA692647 | Illumina 16S V4 (515F/806R)      | 14 | 17 |
| rep2020     | PRJNA657256 | Illumina MiSeq 16S V3-V4 (338F/806R) | 30 | 14 |
| shotgun2022 | PRJNA856121 | Illumina shotgun WGS (75 bp)     | 25 | 11 |

## Layout
- `fusobacterium_demodex_metaanalysis.tex` -- the manuscript
- `scripts/dada2/` -- full pipeline (16S DADA2, shotgun Kraken2/Bracken, meta-analysis)
- `results/`       -- figures and per-cohort / meta-analysis tables
- `data/`          -- cached genus tables and intermediates

## Pipeline (run order)
**16S cohorts (discovery, rep2020):**
```
bash    scripts/dada2/fetch_fastqs.sh PRJNA692647   # discovery FASTQs + SILVA + sheet
bash    scripts/dada2/fetch_657256.sh               # rep2020 FASTQs + sheet
bash    scripts/dada2/cut_primers.sh                # cutadapt 515F/806R (V4)
bash    scripts/dada2/cut_primers_v3v4.sh           # cutadapt 338F/806R (V3-V4)
Rscript scripts/dada2/run_dada2.R                   # DADA2 -> ASV genus tables
```
**Shotgun cohort (shotgun2022):**
```
bash   scripts/dada2/fetch_856_conj.sh   # conjunctival WGS runs
bash   scripts/dada2/fetch_k2db.sh       # Kraken2 standard database
bash   scripts/dada2/build_tools.sh      # build Kraken2 + Bracken
bash   scripts/dada2/classify_856.sh     # Kraken2 classify + Bracken (genus)
python scripts/dada2/agg_bracken.py      # aggregate per-sample genus table
python scripts/dada2/sensitivity_856.py  # host-DNA / depth confounder tests
```
**Meta-analysis + figure:**
```
python scripts/dada2/cross_cohort_concordance.py \
    --cohort <discovery_out>:<discovery_sheet>:discovery \
    --cohort <rep_out>:<rep_sheet>:rep2020 \
    --cohort <shotgun_out>:<shotgun_sheet>:shotgun2022 \
    --min-prev 0.5 --out results/meta3.md
python scripts/dada2/make_forest.py      # -> results/fig_meta_auc.png
```

## Methods
- **16S:** cutadapt primer removal; DADA2 denoising to amplicon sequence
  variants (ASVs); SILVA v138.1 taxonomy; collapsed to genus.
- **Shotgun:** Kraken2 (standard DB) + Bracken genus re-estimation, restricted
  to the Bacteria lineage so abundances are comparable to the 16S cohorts.
- **Meta-analysis:** per-cohort tie-corrected Mann-Whitney signed *z*;
  sqrt(*n*)-weighted Stouffer combined *Z* over the genera present in >= 50% of
  samples in *every* cohort; Benjamini-Hochberg FDR; leave-one-cohort-out
  robustness. Concordant = same direction in all cohorts.

## References
Manuscript bibliography is fetched programmatically from Crossref via
`scripts/dada2/refs_fetch.py` (no hand-entered author/year/DOI metadata).

## Notes
- Networked steps are read-only downloads from ENA/NCBI.
- All statistics are pure-Python (Mann-Whitney *U*, Stouffer *Z*,
  Benjamini-Hochberg, Spearman) in `cross_cohort_concordance.py` and
  `sensitivity_856.py`.
- `scripts/` also retains the earlier single-cohort exploratory code; the
  three-cohort meta-analysis above supersedes it.
