#!/usr/bin/env bash
# Remove 338F/806R V3-V4 primers with cutadapt before DADA2 (PRJNA657256).
# Frontiers Med 2020: 338F=ACTCCTACGGGAGGCAGCA, 806R=GGACTACHVGGGTWTCTAAT.
# ENA layout for this project: *_1 = forward (R1), *_2 = reverse (R2).
# Non-anchored 5' adapters (-g/-G) tolerate any leading spacer.
set -eu

WORK="${WORK:-$HOME/rosacea_657256}"
RAW="$WORK/raw"; TRIM="$WORK/trimmed"
mkdir -p "$TRIM"
CUTADAPT="${CUTADAPT:-$HOME/.local/bin/cutadapt}"

FWD_PRIMER="ACTCCTACGGGAGGCAGCA"      # 338F
REV_PRIMER="GGACTACHVGGGTWTCTAAT"     # 806R

log="$WORK/cutadapt.log"; : > "$log"
runs=$(ls "$RAW"/*_1.fastq.gz | sed -E 's#.*/##; s/_1\.fastq\.gz//')
for r in $runs; do
  fout="$TRIM/${r}_F.fastq.gz"; rout="$TRIM/${r}_R.fastq.gz"
  if [ -s "$fout" ] && [ -s "$rout" ]; then continue; fi
  "$CUTADAPT" -g "$FWD_PRIMER" -G "$REV_PRIMER" \
    --discard-untrimmed -e 0.15 --minimum-length 50 -j 0 \
    -o "$fout" -p "$rout" \
    "$RAW/${r}_1.fastq.gz" "$RAW/${r}_2.fastq.gz" \
    >> "$log" 2>&1
  pass=$(grep -m1 'Pairs written (passing filters):' "$log" | tail -1)
  echo "  $r  $pass"
done
echo "trimmed files: $(ls "$TRIM" | wc -l)"
echo "CUTPRIMERS_OK"
