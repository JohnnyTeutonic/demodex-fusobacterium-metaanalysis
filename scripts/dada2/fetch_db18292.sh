#!/usr/bin/env bash
# Fetch PRJDB18292 (Nakatsuji 2024 JID facial rosacea ivermectin cohort):
#   24 single-end full-length 16S (V1-V9, LoopSeq/synthetic-long-read, ~1.5 kb).
# Builds a 2x2 sample sheet (site x time) from the ENA SAMPLE descriptions that
# were already harvested into data/probe_PRJDB18292_runs.tsv.
#
# Layout under $WORK (default ~/db18292):
#   raw/<run>.fastq.gz
#   sample_sheet.tsv   (run \t site \t time \t group)
set -eu

WORK="${WORK:-$HOME/db18292}"
REPO="${REPO:-/mnt/c/Users/jonat/OneDrive/Documents/research_portfolio_complete/rosacea_demodex}"
PROBE="$REPO/data/probe_PRJDB18292_runs.tsv"
RAW="$WORK/raw"
mkdir -p "$RAW"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"

echo "== building 2x2 sample sheet from $PROBE =="
SS="$WORK/sample_sheet.tsv"
printf 'run\tsite\ttime\tgroup\n' > "$SS"
# probe TSV columns include: run_accession (col1) ... sample_desc_x (the SAMPLE description)
# Find the sample_desc_x column index from the header, then classify.
awk -F'\t' '
  NR==1{ for(i=1;i<=NF;i++){ if($i=="sample_desc_x") dc=i; if($i=="run_accession") rc=i } next }
  {
    run=$rc; d=tolower($(dc));
    site = (d ~ /non-lesional/) ? "nonlesional" : "lesional";
    time = (d ~ /pre-treatment/) ? "pre" : ((d ~ /post-treatment/) ? "post" : "NA");
    print run"\t"site"\t"time"\t"site"_"time
  }' "$PROBE" >> "$SS"

echo "== group counts =="
awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) print "  "k": "c[k]}' "$SS"

echo "== fetching fastq_ftp list for PRJDB18292 =="
RUNTBL="$WORK/runtable.tsv"
curl -s "$ENA?accession=PRJDB18292&result=read_run&fields=run_accession,fastq_ftp&format=tsv" -o "$RUNTBL"

echo "== downloading 24 single-end FASTQs (parallel) =="
awk -F'\t' 'NR>1{print $2}' "$RUNTBL" | tr ';' '\n' | grep -E '\.fastq\.gz$' \
  | sed 's#^#https://#' \
  | xargs -P 6 -I{} bash -c 'f="$1"; b=$(basename "$f"); [ -s "'"$RAW"'/$b" ] || curl -s -o "'"$RAW"'/$b" "$f"' _ {}
echo "  downloaded files: $(ls "$RAW" | wc -l)"
echo "FETCH_OK"
