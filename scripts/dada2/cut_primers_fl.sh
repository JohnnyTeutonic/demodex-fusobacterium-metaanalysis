#!/usr/bin/env bash
# Strip full-length 16S primers (27F / 1492R) from PRJDB18292 single-end reads.
#   27F      = AGAGTTTGATCMTGGCTCAG          (5' adapter, required)
#   1492R rc = AAGTCGTAACAAGGTARCCGTA        (3' adapter, optional)
# --revcomp lets cutadapt orient any reverse-strand reads to forward.
set -eu
WORK="${WORK:-$HOME/db18292}"
RAW="$WORK/raw"; TRIM="$WORK/trimmed"
mkdir -p "$TRIM"
CUTADAPT="${CUTADAPT:-$HOME/.local/bin/cutadapt}"

FWD="AGAGTTTGATCMTGGCTCAG"
REV_RC="AAGTCGTAACAAGGTARCCGTA"

log="$WORK/cutadapt_fl.log"; : > "$log"
for f in "$RAW"/*.fastq.gz; do
  r=$(basename "$f" .fastq.gz)
  out="$TRIM/${r}.fastq.gz"
  [ -s "$out" ] && continue
  "$CUTADAPT" --revcomp -g "$FWD" -a "$REV_RC" \
    --discard-untrimmed -e 0.2 --minimum-length 1000 --maximum-length 1700 -j 0 \
    -o "$out" "$f" >> "$log" 2>&1
  pass=$(grep -m1 'Reads written (passing filters):' "$log" | tail -1)
  echo "  $r  $pass"
done
echo "trimmed files: $(ls "$TRIM" | wc -l)"
echo "CUTPRIMERS_FL_OK"
