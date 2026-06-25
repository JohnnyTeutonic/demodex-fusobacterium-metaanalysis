# Rosacea / Demodex bioinformatic deep dive -- findings

## Question
Can public sequencing data give bioinformatic support to the clinical protocol
*"quantify the underlying Demodex mite burden, then grade antimicrobial treatment
to that load"*? Specifically, can a molecular index built from skin 16S data
serve as a non-invasive proxy for mite burden?

## What was done (all pure-Python, stdlib networking; reused the genetic_code
## fetch/parse pattern)
- **E0 Audit** -- NCBI E-utilities: 656 SRA runs / 25 BioProjects, 17 on-topic.
  Found a topical-**ivermectin** before/after 16S set (PRJDB18292), a
  rosacea+Demodex multi-omics set (PRJEB82826), and a **Demodex+/- ocular** 16S
  set (PRJNA692647), among others.
- **E1 Panel** -- fetched 16S references for a Demodex-associated taxon panel
  (*Snodgrassella alvi, Bacillus oleronius, Cutibacterium kroppenstedtii,
  Bartonella quintana*; + *C. acnes / S. epidermidis* as context) and built
  diagnostic k-mer sets, **background-subtracted** against 24 common skin taxa to
  kill conserved-region cross-mapping.
- **E2 Index** -- streamed FASTQs from ENA, attributed reads by diagnostic k-mers
  (>=2 per read), produced a per-sample `mite_index`.
- **E3 Treatment** -- paired pre/post ivermectin (6 subjects x 2 sites).
- **E4 Validation** -- Demodex-infection vs healthy control (14 vs 17).
- **E5 Method** -- held-out sensitivity/specificity characterisation of the panel.
- **E6 ASV (DADA2)** -- faithful DADA2 1.30 reprocessing of PRJNA692647
  (cutadapt primer removal + SILVA v138.1), a k-mer/ASV bridge, and an
  ASV-vs-reference sequence-identity check. (In WSL2; not pure-Python.)

## Results (honest)
1. **Treatment (E3): NULL.** No significant pre->post change in `mite_index` or
   any taxon (all p >= 0.25); direction mixed, one large outlier. A *surface*
   bacterial proxy poorly reflects a *follicular* mite that swabs under-sample,
   and mite lysis can transiently release bacteria onto the surface.
2. **Validation (E4), k-mer index: borderline composite signal.** With a
   body-site-matched background the composite separated Demodex+ from control
   (AUC 0.66; one-sided p ~ 0.074), *apparently* carried by **Snodgrassella
   alvi** (AUC 0.65) and *Bacillus* (AUC 0.58). **This attribution did NOT
   survive ASV resolution -- see E6.**
3. **Method characterised (E5).** Sensitivity 0.74-0.98 on held-out references;
   clean against distant relatives but genus-level cross-reactive against close
   congeners. The held-out specificity control used too few challengers and
   overstated real-world specificity (E6 shows why).
4. **ASV re-analysis (E6): the S. alvi lead is a CROSS-MAPPING ARTEFACT.**
   DADA2 gave 8,551 ASVs (97.4% reads kept). **Zero** are *Snodgrassella*
   (SILVA v138.1 contains the genus -- genuine absence). Our own *S. alvi*
   diagnostic k-mers map only onto commensal **Neisseriaceae** (*Neisseria*,
   *Eikenella*, ...), 92.5-94.5% identical to the *S. alvi* 16S over V4 --
   inter-genus, not the ~99% a true *S. alvi* would show. The genuine
   exact-sequence signal is a single **Neisseria ASV** enriched in Demodex+
   (**AUC 0.70, one-sided p = 0.032**) -- the strongest, species-resolved
   separation in the study. *Bacillus* panel also cross-maps across Bacillota;
   *Bartonella* onto Rhizobiaceae (confirmed).

## Bottom line
- **Methodological headline:** closed-reference k-mer / marker counting on short
  16S is **not species-specific** and can *manufacture* a plausible genus-level
  signal (here, a fake *S. alvi*). Background subtraction + a held-out
  specificity control were **not enough**; only ASV/DADA2 resolution exposed it.
  Demodex-microbiome studies using closed-reference assignment without ASV
  resolution risk false species-level findings.
- **Surface proxy fails on its own terms:** no treatment tracking; the group
  signal was borderline (and the marker carrying it was an artefact). Grading
  treatment should rest on **direct mite quantification**, not a surface readout.
- **New, defensible lead:** a single **Neisseria ASV** enriched in Demodex+
  ocular surfaces (AUC 0.70, p = 0.032), species-resolved and (to our knowledge)
  previously unreported. A curated ASV-based *Neisseria* assay with
  follicle-targeted sampling is the recommended next experiment.
- Does **not** corroborate the facial *S. alvi* observation of Homey et al.
  (JID 2025): different body site (ocular here), and at exact-sequence
  resolution the ocular signal is *Neisseria*, not *S. alvi*.

## Limitations
- The *Neisseria* lead rests on one ocular cohort (14 vs 17) and one dominant
  ASV; p = 0.032 is nominal (uncorrected). May be ocular-specific.
- Treatment arm small (6 pairs). Facial adjudication cohort PRJEB82826 (n=96,
  cyanoacrylate follicle; Homey/JID 2025) has **no public infestation labels** --
  labelled facial ASV validation blocked on the authors' sample key.
- Taxonomy rests on SILVA v138.1; a different reference could rename, not
  resurrect, the absent *Snodgrassella*.
- See `data/dada2_PRJNA692647/E6_SUMMARY.md` and `scripts/dada2/` for the full
  ASV workflow and bridge.

## Reproduce
See `README.md`. Scripts are independently `--self-test`-able; networked steps
are cached under `data/`.
