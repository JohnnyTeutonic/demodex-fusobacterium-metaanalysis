#!/usr/bin/env Rscript
# Classify PRJDB18292 full-length 16S reads to genus (SILVA 138.1) and build a
# genus-by-sample count matrix. Per sample: derepFastq -> uniques+abundances;
# pool global uniques -> assignTaxonomy once -> sum counts by genus per sample.
suppressWarnings(.libPaths(c(path.expand("~/R/library"), .libPaths())))
suppressMessages(library(dada2))

work  <- Sys.getenv("WORK", path.expand("~/db18292"))
trim  <- file.path(work, "trimmed")
out   <- file.path(work, "out"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
ref   <- Sys.getenv("REF_TRAIN", path.expand("~/rosacea_dada2/ref/silva_nr99_v138.1_train_set.fa.gz"))
threads <- TRUE

fqs <- sort(list.files(trim, pattern = "\\.fastq\\.gz$", full.names = TRUE))
samples <- sub("\\.fastq\\.gz$", "", basename(fqs))
cat("samples:", length(samples), "\n")

derep <- list(); allseqs <- character(0)
for (i in seq_along(fqs)) {
  d <- derepFastq(fqs[i])
  ab <- d$uniques            # named vector: seq -> count
  derep[[samples[i]]] <- ab
  allseqs <- union(allseqs, names(ab))
}
cat("global unique full-length seqs:", length(allseqs), "\n")

cat("== assignTaxonomy (SILVA 138.1) ==\n")
tax <- assignTaxonomy(allseqs, ref, multithread = threads, tryRC = TRUE,
                      minBoot = 50)
genus <- tax[, "Genus"]; family <- tax[, "Family"]
names(genus) <- allseqs; names(family) <- allseqs

g_keys <- sort(unique(ifelse(is.na(genus), "Unclassified", genus)))
mat <- matrix(0, nrow = length(samples), ncol = length(g_keys),
              dimnames = list(samples, g_keys))
fam_of <- tapply(family, ifelse(is.na(genus), "Unclassified", genus),
                 function(x) names(sort(table(x), decreasing = TRUE))[1])
for (s in samples) {
  ab <- derep[[s]]
  gk <- ifelse(is.na(genus[names(ab)]), "Unclassified", genus[names(ab)])
  agg <- tapply(ab, gk, sum)
  mat[s, names(agg)] <- agg
}

write.table(data.frame(sample = rownames(mat), mat, check.names = FALSE),
            file.path(out, "genus_by_sample.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(genus = names(fam_of), family = as.character(fam_of)),
            file.path(out, "genus_family.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("CLASSIFY_FL_OK ->", file.path(out, "genus_by_sample.tsv"), "\n")
