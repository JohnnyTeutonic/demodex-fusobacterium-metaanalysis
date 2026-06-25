#!/usr/bin/env bash
# Download all PRJNA692647 paired FASTQs + SILVA v138.1 DADA2 refs, and build a
# sample sheet (run_accession -> group) from ENA sample_title.
#
# Layout under $WORK (default ~/rosacea_dada2):
#   raw/<run>_1.fastq.gz (806R/reverse), raw/<run>_2.fastq.gz (515F/forward)
#   ref/silva_nr99_v138.1_train_set.fa.gz, ref/silva_species_assignment_v138.1.fa.gz
#   sample_sheet.tsv  (run_accession \t group \t sample_title)
set -eu

PROJECT="${1:-PRJNA692647}"
WORK="${WORK:-$HOME/rosacea_dada2}"
RAW="$WORK/raw"; REF="$WORK/ref"
mkdir -p "$RAW" "$REF"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"

echo "== fetching run table for $PROJECT =="
RUNTBL="$WORK/runtable.tsv"
curl -s "$ENA?accession=$PROJECT&result=read_run&fields=run_accession,sample_title,fastq_ftp&format=tsv" -o "$RUNTBL"
n=$(($(wc -l < "$RUNTBL") - 1))
echo "  $n runs"

# Build sample sheet from sample_title prefix.
SS="$WORK/sample_sheet.tsv"
echo -e "run_accession\tgroup\tsample_title" > "$SS"
awk -F'\t' 'NR>1{
  t=tolower($2); g="other";
  if (t ~ /demodex/) g="demodex";
  else if (t ~ /control/ || t ~ /healthy/) g="control";
  print $1"\t"g"\t"$2
}' "$RUNTBL" >> "$SS"
echo "== sample sheet group counts =="
awk -F'\t' 'NR>1{c[$2]++} END{for(k in c) print "  "k": "c[k]}' "$SS"

# Download both mates for every run (parallel, resumable).
echo "== downloading FASTQs (parallel) =="
awk -F'\t' 'NR>1{print $3}' "$RUNTBL" | tr ';' '\n' | grep -E '_[12]\.fastq\.gz$' \
  | sed 's#^#https://#' \
  | xargs -P 6 -I{} bash -c 'f="$1"; b=$(basename "$f"); [ -s "'"$RAW"'/$b" ] || curl -s -o "'"$RAW"'/$b" "$f"' _ {}
echo "  downloaded files: $(ls "$RAW" | wc -l)"

# SILVA v138.1 DADA2-formatted references (McLaren 2020, Zenodo 4587955).
echo "== fetching SILVA v138.1 references =="
declare -A URLS=(
  ["silva_nr99_v138.1_train_set.fa.gz"]="https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
  ["silva_species_assignment_v138.1.fa.gz"]="https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz"
)
for f in "${!URLS[@]}"; do
  if [ -s "$REF/$f" ]; then echo "  have $f"; else
    echo "  downloading $f"; curl -L -s -o "$REF/$f" "${URLS[$f]}"
  fi
done
ls -lh "$REF"
echo "FETCH_OK"
