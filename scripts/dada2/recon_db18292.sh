#!/usr/bin/env bash
# Recon: confirm PRJDB18292 read structure (length, platform, primer region).
set -eu
WORK="${WORK:-$HOME/db18292}"
RAW="$WORK/raw"
mkdir -p "$RAW"
cd "$RAW"
RUN="DRR574445"   # smallest run
URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/DRR574/$RUN/$RUN.fastq.gz"
[ -s "$RUN.fastq.gz" ] || curl -s -o "$RUN.fastq.gz" "$URL"

echo "=== n reads ==="
expr $(zcat "$RUN.fastq.gz" | wc -l) / 4

echo "=== read-length distribution (100 bp bins) ==="
zcat "$RUN.fastq.gz" | sed -n '2~4p' | awk '{print int(length/100)*100}' | sort -n | uniq -c

echo "=== leading 60 bp of first 5 reads (primer check) ==="
zcat "$RUN.fastq.gz" | sed -n '2~4p' | cut -c1-60 | head -5

echo "=== trailing 30 bp of first 3 reads ==="
zcat "$RUN.fastq.gz" | sed -n '2~4p' | head -3 | rev | cut -c1-30 | rev
echo "RECON_OK"
