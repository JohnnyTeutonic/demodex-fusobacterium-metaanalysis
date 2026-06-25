#!/usr/bin/env Rscript
# Install dada2 (+ BiocManager) into a user library so no sudo is needed.
# System -dev libraries (installed via apt) provide the compile toolchain.
# Bioconductor 3.18 is the release matched to R 4.3.x.

user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
options(Ncpus = max(1L, parallel::detectCores()),
        repos = c(CRAN = "https://cloud.r-project.org"))

cat("R:", R.version.string, "\n")
cat("install target:", .libPaths()[1], "\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = user_lib)
}
BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)

if (!requireNamespace("dada2", quietly = TRUE)) {
  BiocManager::install("dada2", ask = FALSE, update = FALSE, lib = user_lib)
}

suppressMessages(library(dada2))
cat("DADA2 version:", as.character(packageVersion("dada2")), "\n")
cat("INSTALL_OK\n")
