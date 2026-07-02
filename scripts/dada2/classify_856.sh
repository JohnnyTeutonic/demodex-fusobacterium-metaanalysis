#!/usr/bin/env bash
# Classify PRJNA856121 conjunctival shotgun reads with Kraken2 (Standard-8 DB),
# then Bracken genus-level abundance re-estimation (75-mer distrib for 75bp reads).
# --memory-mapping keeps the 7.5GB DB warm in page cache across samples.
set -eu
WORK="${WORK:-$HOME/rosacea_856}"
DB="${DB:-$HOME/k2db}"
K2="$HOME/tools/kraken2/kraken2"
BRACKEN="$HOME/tools/Bracken/bracken"
OUT="$WORK/kraken"; mkdir -p "$OUT/reports" "$OUT/bracken"
THREADS="$(nproc)"

tail -n +2 "$WORK/sample_sheet.tsv" | cut -f1 | while read -r run; do
  fq="$WORK/raw/$run.fastq.gz"
  rep="$OUT/reports/$run.kreport"
  brk="$OUT/bracken/$run.bracken"
  [ -s "$fq" ] || { echo "MISSING $fq"; continue; }
  if [ ! -s "$rep" ]; then
    "$K2" --db "$DB" --memory-mapping --threads "$THREADS" --gzip-compressed \
      --report "$rep" --output /dev/null "$fq" 2> "$OUT/reports/$run.log"
  fi
  if [ ! -s "$brk" ]; then
    "$BRACKEN" -d "$DB" -i "$rep" -o "$brk" -r 75 -l G -t 1 \
      > "$OUT/bracken/$run.blog" 2>&1 || echo "bracken-warn $run"
  fi
  echo "done $run"
done
echo "reports: $(ls "$OUT/reports"/*.kreport 2>/dev/null | wc -l), bracken: $(ls "$OUT/bracken"/*.bracken 2>/dev/null | wc -l)"
echo CLASSIFY_OK
