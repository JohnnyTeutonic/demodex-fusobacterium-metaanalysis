#!/usr/bin/env bash
# Recon PRJNA657256 (Frontiers Med 2020 ocular Demodex blepharitis, V3-V4 338F/806R,
# MiSeq 2x250). Pull run table with grouping-relevant fields and read stats, then
# download the head of the smallest run to confirm read length / primer position.
set -eu

PROJECT="PRJNA657256"
WORK="${WORK:-$HOME/rosacea_657256}"
mkdir -p "$WORK"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"

RUNTBL="$WORK/runtable.tsv"
curl -s "$ENA?accession=$PROJECT&result=read_run&fields=run_accession,sample_title,sample_alias,library_layout,read_count,base_count,fastq_ftp&format=tsv" -o "$RUNTBL"
echo "== run table ($PROJECT) =="
cat "$RUNTBL"
echo
echo "== n runs =="
wc -l "$RUNTBL"
echo "RECON_OK"
