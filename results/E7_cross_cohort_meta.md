# E7 — Cross-cohort meta-analysis of the ocular *Demodex* microbiome

## Motivation
The single-cohort *Neisseria* lead (discovery PRJNA692647: genus AUC 0.689, one-sided
p=0.039) did **not** replicate in an independent ocular cohort (PRJNA657256:
AUC 0.321, opposite direction). Rather than chase one taxon, we reframed the
question positively: **is there ANY reproducible microbial signature of ocular
*Demodex* disease across independent cohorts?**

## Method
- Both cohorts reprocessed with the **same** DADA2 + SILVA 138.1 pipeline
  (discovery = V4 515F/806R; rep2020 = V3–V4 338F/806R). Because amplicon regions
  differ, the region-robust unit of comparison is the **genus**.
- `cross_cohort_concordance.py`: per-genus Mann–Whitney (demodex vs control) in each
  cohort; genera filtered to ≥50% prevalence in **every** cohort (this suppresses
  study-specific contaminants, which are the dominant per-cohort "signal" otherwise);
  Stouffer combined Z (√n weighting); Bonferroni over tested genera.

## Cohorts
| cohort | accession | region | n (DB / control) | site | status |
|---|---|---|---|---|---|
| discovery | PRJNA692647 | V4 | 14 / 17 | conjunctiva | processed |
| rep2020 (Yan) | PRJNA657256 | V3–V4 | 30 / 14 | conjunctiva | processed |
| P&V 2024 (Zou) | GSA CRA014100 | V3–V4 | 25 / 21 | eyelash | **not fetched — GSA/NGDC unreliable from AU** |
| fcimb 2022 | PRJNA856121 | shotgun | 25 / 11 | conjunctiva+meibum | not processed (needs classifier) |

## Result (2 cohorts)
Prevalence-filtered genera tested: 68. **Nothing survives Bonferroni at n=2 cohorts**
(expected in this underpowered regime). Concordant nominal trends (same direction in
both cohorts):

| genus | direction in Demodex | combined Z | combined p |
|---|---|---|---|
| **Fusobacterium** | down | −2.64 | 0.008 |
| **Rothia** | up | +2.15 | 0.031 |
| **Streptococcus** | up | +1.94 | 0.052 |
| Lachnospiraceae NK4A136 | up | +1.91 | 0.056 |
| Burkholderia group | down | −1.61 | 0.107 |
| **Corynebacterium** | up | +1.57 | 0.117 |

## External corroboration (published biomarkers, other cohorts)
- **Corynebacterium up** in Demodex: our discovery + rep2020 + **P&V 2024** (Corynebacterium_1
  biomarker) + Lee 2012 (pyrosequencing). Mechanistically grounded — the *Demodex*
  endobacterium is *Corynebacterium kroppenstedtii* subsp. *demodicis*.
- **Burkholderia down** in Demodex: P&V 2024 report Burkholderia 2.36× lower in DB (their
  top control biomarker, AUC 0.74) — matches our concordant down trend.
- **Lachnospiraceae NK4A136 up**: P&V 2024 DB biomarker — matches.

## Honest interpretation
- Corynebacterium enrichment is the most cross-validated, mechanistically coherent
  reproducible signal, but is only *nominally* trending in our 2-cohort reanalysis.
- Fusobacterium depletion is our strongest *statistical* signal (p=0.008) and concordant,
  but has weaker prior literature — it would be a more novel claim needing more cohorts.
- To reach publishable strength, add ≥1–2 independent cohorts. Stouffer Z grows with
  each concordant cohort, so a genuine signal (e.g. Corynebacterium up, Fusobacterium
  down) should clear multiple-testing correction at 3–4 cohorts.

## UPDATE — 3-cohort result (added shotgun PRJNA856121)

Shotgun cohort processed: 36 conjunctival samples (25 DB + 11 control), Kraken2
Standard-8 + Bracken genus abundances (75-mer distrib), restricted to bacterial
genera. Genera passing the >=50%-in-every-cohort prevalence filter: 38.
Output: `~/rosacea_meta3.md`.

| genus | dir | Z(3) | p(3) | Bonf p | AUC disc | AUC rep | AUC shot | concordant |
|---|---|---|---|---|---|---|---|---|
| **Fusobacterium** | down | **-2.93** | **0.0034** | 0.131 | 0.332 | 0.305 | 0.360 | **3/3** |
| Bradyrhizobium | up | +2.83 | 0.0046 | 0.176 | 0.668 | 0.520 | 0.844 | 3/3 (likely contaminant) |
| Brevundimonas | down | -2.11 | 0.035 | 1.0 | 0.366 | 0.452 | 0.291 | 3/3 (possible contaminant) |
| Corynebacterium | up | +1.40 | 0.163 | 1.0 | 0.685 | 0.555 | 0.520 | 3/3 (endobacterium) |
| Rothia | up | +0.94 | 0.347 | 1.0 | 0.729 | 0.590 | 0.347 | NO (shotgun disagrees) |
| Streptococcus | up | +0.93 | 0.352 | 1.0 | 0.664 | 0.617 | 0.376 | NO (shotgun disagrees) |
| Neisseria | mixed | -1.79 | 0.073 | 1.0 | 0.689 | 0.321 | 0.222 | NO (original hypothesis refuted) |

### Key takeaways
- **Fusobacterium depletion is the robust, reproducible signal**: same direction in
  all 3 independent cohorts across 2 methodologies (16S x2 + shotgun x1). The
  combined Z GREW from -2.64 (2 cohorts) to -2.93 (3 cohorts) -- the signature of a
  real effect gaining power with data. Not a canonical reagent contaminant.
- **Rothia and Streptococcus did not replicate** in shotgun -> earlier 2-cohort leads
  were likely amplicon-specific / noise. (Correctly falsified.)
- **Corynebacterium** stays concordant-up in all 3 but weak; mechanistically tied to
  the Demodex endobacterium C. kroppenstedtii.
- **Nothing survives Bonferroni at n=3** (Fusobacterium Bonf p=0.131). A 4th cohort
  (GSA P&V 2024 eyelash) is the clear next step to test whether Fusobacterium clears
  multiple-testing correction.
- Caution: Bradyrhizobium/Brevundimonas rank high but are common environmental/reagent
  contaminants; do not over-interpret.

## Sensitivity check (shotgun cohort, `sensitivity_856.py`) -- PASSED

Ocular shotgun is ~99% host DNA, so we tested whether Fusobacterium depletion is a
depth/host artifact (Kraken2 report clade counts; demodex n=25, control n=11):

| covariate | demodex med | control med | MWU p | verdict |
|---|---|---|---|---|
| host fraction (Homo/total) | 0.066 | 0.127 | 0.99 | no group difference |
| bacterial reads (Bacteria clade) | 3.09M | 2.47M | 0.82 | no group difference |
| total reads (depth) | 22.6M | 16.5M | 0.022 | DB HIGHER -> biases AGAINST depletion |

Fusobacterium depletion is robust across normalisations:
- (a) /bacterial reads: AUC 0.360 ; (b) /total reads: AUC 0.331 ; (c) RPM: AUC 0.331
- (d) prevalence: **detected in 21/25 (84%) DB vs 11/11 (100%) control**
- Spearman rho(Fuso, host-fraction)=+0.20, rho(Fuso, depth)=-0.15 (both weak) ->
  not tracking a technical covariate.

Conclusion: the depletion is NOT a compositional/host-DNA artifact; DB even has higher
total depth (which would inflate rare-taxon detection), making the result conservative.

## Cohort 4 (Zhang/Zou et al. 2024, Parasites & Vectors, CRA014100) -- BLOCKED

Eyelash V3-V4, 25 DB + 21 control. Attempted to add as cohort 4:
- GSA/NGDC hosts unreachable from this network: `ngdc.cncb.ac.cn` returns http=000
  (timeout); `download.cncb.ac.cn/gsa*/CRA014100/` returns 404 (autoindex off, needs
  exact CRR accessions); FTP `download.big.ac.cn` returns 550. Probe: `gsa_probe.sh`.
- No per-sample data published: data-availability = "deposited in GSA CRA014100" only;
  sole supplementary is a Burkholderia/survival .docx; the 26 differential genera are in
  Fig. 6a (image, not a table). **Fusobacterium is never reported in the paper.**
=> Cannot be folded into the pooled Fusobacterium meta-analysis without the raw fastqs.

External (non-pooled) corroboration extracted from the paper:
- Corynebacterium_1 UP in DB (6.71%; t=-1.77, p=0.084) -> matches our concordant
  Corynebacterium-up (now 4/4 cohorts agree on DIRECTION).
- Burkholderia-Caballeronia-Paraburkholderia DOWN in DB (2.36x, p=0.002; AUC 0.745)
  -> matches our concordant Burkholderia-down trend.

To actually add this cohort, one of: (a) download CRA014100 on a network with NGDC
access and run it through the DADA2 16S pipeline (like discovery/rep2020); (b) request
the per-sample genus table / raw fastqs from the corresponding author (Ting Wang).

Note: our "discovery" cohort PRJNA692647 == Liang et al. 2021 (14 DI + 17 healthy),
so no additional easy SRA-hosted ocular Demodex 16S cohort remains beyond the three used.

## Next steps (pending)
1. Obtain CRA014100 raw data via alternate network/author to test Fusobacterium at n=4.
2. Literature/mechanism check on ocular Fusobacterium depletion in Demodex disease.

## Files
- `scripts/dada2/fetch_657256.sh`, `cut_primers_v3v4.sh` — rep2020 processing
- `scripts/dada2/analyze_neisseria.py` — per-cohort target-taxa test
- `scripts/dada2/cross_cohort_concordance.py` — meta-analysis
