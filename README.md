# rosacea_demodex

Bioinformatic support for the Demodex-burden / graded-treatment hypothesis in
rosacea. The discovery pipeline is pure-Python (stdlib networking + stats);
the confirmatory ASV step (E6) uses a faithful DADA2/R workflow under WSL2.
Reuses the fetch/parse pattern from `../genetic_code`.

**Headline (after the E6 ASV upgrade):** the closed-reference $k$-mer index's
apparent *S. alvi* signal is a conserved-region **cross-mapping artefact** --
DADA2/SILVA finds zero *Snodgrassella* ASVs and the *S. alvi* k-mers land on
commensal Neisseriaceae (92.5-94.5% identity). The real, species-resolved
Demodex+ signal is a single **Neisseria ASV** (AUC 0.70, p=0.032). Closed-
reference marker counting on short 16S is not species-specific; ASV resolution
is required. See `FINDINGS.md` and `data/dada2_PRJNA692647/E6_SUMMARY.md`.

## Layout
- `PLAN.md`      -- project plan and panel rationale
- `FINDINGS.md`  -- headline results and conclusions (read this)
- `scripts/`     -- pipeline (each has `--self-test`)
- `data/`        -- cached fetches, k-mer panels, per-sample indices
- `results/`     -- E0/E3/E4 write-ups + test tables

## Pipeline (run order)
```
# E0  data-availability audit (NCBI E-utilities)
python scripts/sra_audit.py
python scripts/summarize_audit.py            # -> results/project_summary.tsv

# E1  reference 16S + diagnostic, body-site-matched background-subtracted k-mers
python scripts/fetch_references.py           # -> data/refs, data/kmers
python scripts/validate_panel.py             # E5: sensitivity/specificity control

# E2  per-project metadata + molecular mite-burden index (streams ENA FASTQ)
python scripts/fetch_metadata.py --project PRJDB18292
python scripts/build_index.py   --project PRJDB18292
python scripts/fetch_metadata.py --project PRJNA692647
python scripts/build_index.py   --project PRJNA692647

# E3/E4  analyses
python scripts/analyze_treatment.py          # ivermectin before/after
python scripts/analyze_validation.py         # Demodex+ vs control

# Figure  ROC: naive (skin-only) vs corrected (body-site-matched)
#   first regenerate the naive ocular index into a separate file:
python scripts/fetch_references.py --skin-only --out-dir data/kmers_naive
python scripts/build_index.py --project PRJNA692647 \
       --kmers-dir data/kmers_naive --out data/index_PRJNA692647_naive.tsv
python scripts/make_roc_figure.py            # -> results/fig_roc_validation.pdf

# E6  exact-sequence (ASV) re-analysis -- DADA2/R under WSL2 (see scripts/dada2/)
Rscript scripts/dada2/install_dada2.R        # one-time: dada2 into a user lib
bash    scripts/dada2/fetch_fastqs.sh PRJNA692647   # FASTQs + SILVA refs + sheet
bash    scripts/dada2/cut_primers.sh         # cutadapt 515F/806R removal
Rscript scripts/dada2/run_dada2.R            # -> data/dada2_PRJNA692647/*
python  scripts/dada2/analyze_asv.py         # ASV group test (E6)
python  scripts/dada2/bridge_kmer_asv.py     # k-mer panel applied to ASVs (E6b)
python  scripts/dada2/identity_check.py      # ASV-vs-S.alvi identity (E6c)
```

Key result: the $k$-mer composite gave a borderline Demodex+ vs control signal
(AUC 0.66) apparently carried by *S. alvi* -- but the **E6 DADA2/ASV upgrade
overturned it** (cross-mapping artefact; the real lead is a *Neisseria* ASV,
AUC 0.70, p=0.032). The ivermectin arm is null. See `FINDINGS.md`.

## Self-tests (offline, no network)
```
python scripts/sra_audit.py        --self-test
python scripts/summarize_audit.py  --self-test
python scripts/fetch_references.py --self-test
python scripts/validate_panel.py   --self-test
python scripts/fetch_metadata.py   --self-test
python scripts/build_index.py      --self-test
python scripts/make_roc_figure.py  --self-test
python scripts/probe_labels.py     --self-test   # PRJEB82826 label probe
python scripts/dada2/analyze_asv.py     --self-test
python scripts/dada2/bridge_kmer_asv.py --self-test
```

## Notes
- Networked steps are read-only GETs to NCBI E-utilities and EBI ENA.
- `build_index.py` retries each FASTQ up to 3x on transient timeouts.
- Index = closed-reference diagnostic-k-mer attribution (>=2 k-mers/read),
  a relative proxy, NOT an ASV census. See `FINDINGS.md` for caveats.
