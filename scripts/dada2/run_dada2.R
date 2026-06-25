#!/usr/bin/env Rscript
# DADA2 paired-end pipeline for PRJNA692647 V4 (515F/806R), primers already
# removed by cutadapt (scripts/dada2/cut_primers.sh).
#
# filterAndTrim -> learnErrors -> dada -> mergePairs -> makeSequenceTable
#   -> removeBimeraDenovo -> assignTaxonomy(SILVA 138.1) -> addSpecies
#
# Env vars:
#   WORK     working dir (default ~/rosacea_dada2)
#   TRUNC_F  forward truncation length (default 170)
#   TRUNC_R  reverse truncation length (default 150)
#   THREADS  multithread (default TRUE = all cores)
#
# Outputs under $WORK/out:
#   asv_seqs.fasta, asv_table.tsv, taxonomy.tsv, track.tsv

suppressWarnings(.libPaths(c(path.expand("~/R/library"), .libPaths())))
suppressMessages(library(dada2))

work    <- Sys.getenv("WORK", path.expand("~/rosacea_dada2"))
truncF  <- as.integer(Sys.getenv("TRUNC_F", "170"))
truncR  <- as.integer(Sys.getenv("TRUNC_R", "150"))
threads <- Sys.getenv("THREADS", "TRUE")
threads <- if (threads %in% c("TRUE","FALSE")) as.logical(threads) else as.integer(threads)

trim <- file.path(work, "trimmed")
filt <- file.path(work, "filtered")
out  <- file.path(work, "out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
ref_train <- file.path(work, "ref", "silva_nr99_v138.1_train_set.fa.gz")
ref_spec  <- file.path(work, "ref", "silva_species_assignment_v138.1.fa.gz")

fnFs <- sort(list.files(trim, pattern = "_F\\.fastq\\.gz$", full.names = TRUE))
fnRs <- sort(list.files(trim, pattern = "_R\\.fastq\\.gz$", full.names = TRUE))
stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))
sample.names <- sub("_F\\.fastq\\.gz$", "", basename(fnFs))
cat("samples:", length(sample.names), "| truncLen =", truncF, truncR, "\n")

filtFs <- file.path(filt, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names; names(filtRs) <- sample.names

cat("== filterAndTrim ==\n")
flt <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(truncF, truncR),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = threads)
print(head(flt))

# Some samples may lose all reads; keep only those with surviving filtered files.
keep <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]
sample.names <- sample.names[keep]

cat("== learnErrors (forward) ==\n");  errF <- learnErrors(filtFs, multithread = threads)
cat("== learnErrors (reverse) ==\n");  errR <- learnErrors(filtRs, multithread = threads)

cat("== dada ==\n")
ddF <- dada(filtFs, err = errF, multithread = threads)
ddR <- dada(filtRs, err = errR, multithread = threads)

cat("== mergePairs ==\n")
mergers <- mergePairs(ddF, filtFs, ddR, filtRs, verbose = TRUE)

seqtab <- makeSequenceTable(mergers)
cat("ASVs before chimera removal:", ncol(seqtab), "\n")
cat("merged length distribution:\n"); print(table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = threads, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab.nochim),
    "| frac reads kept:", round(sum(seqtab.nochim)/sum(seqtab), 4), "\n")

# ---- track reads through the pipeline ----
getN <- function(x) sum(getUniques(x))
single <- function(d) if (length(sample.names) == 1) getN(d) else sapply(d, getN)
track <- cbind(flt[keep, , drop = FALSE],
               denoisedF = single(ddF), denoisedR = single(ddR),
               merged = single(mergers), nonchim = rowSums(seqtab.nochim))
rownames(track) <- sample.names
write.table(data.frame(sample = rownames(track), track),
            file.path(out, "track.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- taxonomy ----
taxa <- NULL
if (file.exists(ref_train)) {
  cat("== assignTaxonomy (SILVA 138.1) ==\n")
  taxa <- assignTaxonomy(seqtab.nochim, ref_train, multithread = threads, tryRC = TRUE)
  if (file.exists(ref_spec)) {
    cat("== addSpecies ==\n")
    taxa <- addSpecies(taxa, ref_spec, tryRC = TRUE)
  }
} else {
  cat("!! SILVA training set missing; skipping taxonomy.\n")
}

# ---- write ASV outputs with stable IDs ----
seqs <- getSequences(seqtab.nochim)
asv.ids <- paste0("ASV", seq_along(seqs))
fa <- file(file.path(out, "asv_seqs.fasta"), "w")
for (i in seq_along(seqs)) cat(sprintf(">%s\n%s\n", asv.ids[i], seqs[i]), file = fa)
close(fa)

tab <- t(seqtab.nochim); rownames(tab) <- asv.ids
write.table(data.frame(ASV = asv.ids, tab, check.names = FALSE),
            file.path(out, "asv_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(taxa)) {
  tx <- as.data.frame(taxa); rownames(tx) <- NULL
  tx <- cbind(ASV = asv.ids, sequence = seqs, tx)
  write.table(tx, file.path(out, "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  hits <- grep("Snodgrassella|Bartonella|Bacillus", apply(as.data.frame(taxa), 1, paste, collapse=";"))
  cat("\n== ASVs hitting target genera (Snodgrassella/Bartonella/Bacillus):", length(hits), "==\n")
  if (length(hits)) print(cbind(asv.ids[hits], as.data.frame(taxa)[hits, c("Genus"), drop=FALSE]))
}

cat("DADA2_OK -> outputs in", out, "\n")
