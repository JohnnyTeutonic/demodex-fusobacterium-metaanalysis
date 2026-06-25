# E0 -- Data-availability audit (NCBI SRA / BioProject)

Run: `scripts/sra_audit.py` then `scripts/summarize_audit.py`.
656 SRA runs across 25 BioProjects; 17 on-topic (rosacea / Demodex / treatment).
NCBI BioProject esummary is server-side broken (DTD loader error) in both XML
and JSON; study titles recovered via SRA `efetch ... rettype=full`.

## Library-strategy breakdown (656 runs)
- AMPLICON (16S etc.): 532
- WGS (shotgun):        73
- WXS:                  24
- RNA-Seq:              17
- WGA:                  7
- OTHER:                3

## Highest-value projects
| BioProject   | n  | tags              | what it is |
|--------------|----|-------------------|------------|
| PRJDB18292   | 24 | ROSACEA,TREATMENT | Topical **ivermectin** treatment of rosacea changes the bacterial microbiome (antiparasitic before/after) |
| PRJEB82826   | 96 | ROSACEA,DEMODEX   | Multi-omics microbe-host interactions in rosacea (infestation + microbiome + transcriptome; JID 2025) |
| PRJNA1189573 | 93 | ROSACEA,16S       | Skin, blood, stool microbiome in rosacea |
| PRJNA856121  | 72 | DEMODEX           | Metagenomic ocular surface microbiome in Demodex |
| PRJNA1288008 | 66 | ROSACEA           | Skin/blood/stool mycobiome in rosacea |
| PRJNA861245  | 55 | DEMODEX           | Innate type-2 immunity controls follicle commensalism by Demodex (mouse) |
| PRJNA879945  | 42 | DEMODEX,TREATMENT | Ocular microbiota in Demodex blepharitis |
| PRJNA692647  | 31 | DEMODEX,16S       | Demodex infection changes ocular surface microbial communities |
| PRJDB / etc. | .. | DEMODEX           | reference genome (PRJEB13411), transcriptomes (PRJNA510717/630526) |

## Decisions
- **Treatment arm**: PRJDB18292 (ivermectin, n=24) -- small, paired, directly on-thesis.
- **Index build/validate**: PRJEB82826 (rosacea+Demodex multi-omics, n=96) and
  PRJNA692647 (Demodex status vs microbiome, n=31).
- **References**: pull 16S for the Demodex-associated taxon panel + PRJEB13411
  Demodex genome for completeness.

## The processing constraint (honest)
SRA stores raw FASTQ. Standard 16S -> ASV -> taxonomy needs DADA2/QIIME2
(R/conda, heavy, not installed; WSL2 available but a big lift). Tractable
alternative without that stack:
- pull FASTQ via the **ENA** per-run download API (direct links, no prefetch);
- build a **closed-reference marker index**: diagnostic k-mers from the
  Demodex-associated taxa 16S references, counted in reads, normalised per
  library size -> a continuous molecular mite-burden proxy.
This is approximate (k-mer/closed-reference, not full ASV) but reproducible in
pure Python and sufficient for a relative, gradeable index.
