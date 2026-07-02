#!/usr/bin/env bash
# Fetch PRJNA657256 (Frontiers Med 2020 ocular Demodex blepharitis, V3-V4 338F/806R,
# MiSeq 2x250) paired FASTQs and build a sample sheet.
#
# Grouping: sample_alias prefix A_ = Demodex blepharitis (group A), B_ = control (group B).
# Reuses SILVA 138.1 refs from the discovery WORK (~/rosacea_dada2/ref) if present,
# else downloads them.
#
# Layout under $WORK (default ~/rosacea_657256):
#   raw/<run>_1.fastq.gz, raw/<run>_2.fastq.gz
#   ref/silva_nr99_v138.1_train_set.fa.gz, ref/silva_species_assignment_v138.1.fa.gz
#   sample_sheet.tsv  (run_accession \t group \t sample_alias)
set -eu

PROJECT="PRJNA657256"
WORK="${WORK:-$HOME/rosacea_657256}"
RAW="$WORK/raw"; REF="$WORK/ref"
mkdir -p "$RAW" "$REF"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"

echo "== fetching run table for $PROJECT =="
RUNTBL="$WORK/runtable.tsv"
curl -s "$ENA?accession=$PROJECT&result=read_run&fields=run_accession,sample_alias,fastq_ftp&format=tsv" -o "$RUNTBL"
n=$(($(wc -l < "$RUNTBL") - 1))
echo "  $n runs"

# Build sample sheet from sample_alias prefix (A_* demodex, B_* control).
SS="$WORK/sample_sheet.tsv"
echo -e "run_accession\tgroup\tsample_alias" > "$SS"
awk -F'\t' 'NR>1{
  a=$2; g="other";
  if (a ~ /^A_/) g="demodex";
  else if (a ~ /^B_/) g="control";
  print $1"\t"g"\t"a
}' "$RUNTBL" >> "$SS"
echo "== sample sheet group counts =="
awk -F'\t' 'NR>1{c[$2]++} END{for(k in c) print "  "k": "c[k]}' "$SS"

# Download both mates for every run (parallel, resumable).
echo "== downloading FASTQs (parallel) =="
awk -F'\t' 'NR>1{print $3}' "$RUNTBL" | tr ';' '\n' | grep -E '_[12]\.fastq\.gz$' \
  | sed 's#^#https://#' \
  | xargs -P 6 -I{} bash -c 'f="$1"; b=$(basename "$f"); [ -s "'"$RAW"'/$b" ] || curl -s -o "'"$RAW"'/$b" "$f"' _ {}
echo "  downloaded files: $(ls "$RAW" | wc -l)"

# SILVA v138.1 refs: reuse discovery copy if available, else download.
DISC_REF="$HOME/rosacea_dada2/ref"
declare -A URLS=(
  ["silva_nr99_v138.1_train_set.fa.gz"]="https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
  ["silva_species_assignment_v138.1.fa.gz"]="https://zenodo.org/records/4587955/files/silva_species_assignment_v138.1.fa.gz"
)
for f in "${!URLS[@]}"; do
  if [ -s "$REF/$f" ]; then echo "  have $f"; 
  elif [ -s "$DISC_REF/$f" ]; then echo "  reusing $f from discovery"; cp "$DISC_REF/$f" "$REF/$f";
  else echo "  downloading $f"; curl -L -s -o "$REF/$f" "${URLS[$f]}"; fi
done
ls -lh "$REF"
echo "FETCH_OK"
