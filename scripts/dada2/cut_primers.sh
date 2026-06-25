#!/usr/bin/env bash
# Remove 515F/806R V4 primers with cutadapt before DADA2.
# Orientation (verified by inspect_reads.sh):
#   *_2.fastq.gz = forward read, starts with 515F at offset 0
#   *_1.fastq.gz = reverse read, 806R after a ~14 bp heterogeneity spacer
# Non-anchored 5' adapters (-g/-G) tolerate the variable spacer.
set -eu

WORK="${WORK:-$HOME/rosacea_dada2}"
RAW="$WORK/raw"; TRIM="$WORK/trimmed"
mkdir -p "$TRIM"
CUTADAPT="${CUTADAPT:-$HOME/.local/bin/cutadapt}"

FWD_PRIMER="GTGYCAGCMGCCGCGGTAA"      # 515F
REV_PRIMER="GGACTACNVGGGTWTCTAAT"     # 806R

log="$WORK/cutadapt.log"; : > "$log"
runs=$(ls "$RAW"/*_1.fastq.gz | sed -E 's#.*/##; s/_1\.fastq\.gz//')
for r in $runs; do
  fout="$TRIM/${r}_F.fastq.gz"; rout="$TRIM/${r}_R.fastq.gz"
  if [ -s "$fout" ] && [ -s "$rout" ]; then continue; fi
  "$CUTADAPT" -g "$FWD_PRIMER" -G "$REV_PRIMER" \
    --discard-untrimmed -e 0.15 --minimum-length 50 -j 0 \
    -o "$fout" -p "$rout" \
    "$RAW/${r}_2.fastq.gz" "$RAW/${r}_1.fastq.gz" \
    >> "$log" 2>&1
  pass=$(grep -m1 'Pairs written (passing filters):' "$log" | tail -1)
  echo "  $r  $pass"
done
echo "trimmed files: $(ls "$TRIM" | wc -l)"
echo "CUTPRIMERS_OK"
