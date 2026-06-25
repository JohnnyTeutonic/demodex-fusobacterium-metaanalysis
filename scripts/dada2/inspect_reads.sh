#!/usr/bin/env bash
# Inspect a downloaded run: read-length distribution + primer presence/offset.
# No 'pipefail' so that head closing a zcat pipe (SIGPIPE) is not fatal.
set -eu
RAW="${1:-$HOME/rosacea_dada2/raw}"
R1=$(ls "$RAW"/*_1.fastq.gz | head -1)
R2=$(ls "$RAW"/*_2.fastq.gz | head -1)
F515='GTG[CT]CAGC[AC]GCCGCGGTAA'
R806='GGACTAC[ACGTNVH][ACGTNV]GGGT[AT]TCTAAT'

echo "== R1 length distribution (first 10k reads) =="
zcat "$R1" | head -40000 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c
echo "== R2 length distribution (first 10k reads) =="
zcat "$R2" | head -40000 | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c

echo "== R1 806R primer offset (col1=byte offset, 0-indexed; first 3k reads) =="
zcat "$R1" | head -12000 | awk 'NR%4==2' | grep -obE "$R806" | awk -F: '{print $1}' | sort -n | uniq -c | head
echo "== R2 515F primer offset (first 3k reads) =="
zcat "$R2" | head -12000 | awk 'NR%4==2' | grep -obE "$F515" | awk -F: '{print $1}' | sort -n | uniq -c | head

echo "== primer presence counts (first 2k reads) =="
echo -n "R1 with 806R: "; zcat "$R1" | head -8000 | awk 'NR%4==2' | grep -cE "$R806" || true
echo -n "R2 with 515F: "; zcat "$R2" | head -8000 | awk 'NR%4==2' | grep -cE "$F515" || true
