#!/usr/bin/env bash
# Fetch PRJNA856121 (fcimb 2022 shotgun) CONJUNCTIVAL samples only (25 DB + 11 control),
# to match the conjunctival niche of the other cohorts. 75bp single-end WGS.
# Requires $WORK/sample_sheet_full.tsv from build_856_meta.sh.
#
# Outputs:
#   $WORK/sample_sheet.tsv   (run_accession  group)   conjunctiva only
#   $WORK/raw/<run>.fastq.gz
set -eu

WORK="${WORK:-$HOME/rosacea_856}"
RAW="$WORK/raw"; mkdir -p "$RAW"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
FULL="$WORK/sample_sheet_full.tsv"

# Conjunctival-only sample sheet (run_accession \t group)
SS="$WORK/sample_sheet.tsv"
printf 'run_accession\tgroup\n' > "$SS"
awk -F'\t' 'NR>1 && $3=="conjunctiva" && ($2=="demodex"||$2=="control"){print $1"\t"$2}' "$FULL" >> "$SS"
echo "conjunctival samples: $(($(wc -l < "$SS") - 1))"
awk -F'\t' 'NR>1{c[$2]++} END{for(k in c) print "  "k": "c[k]}' "$SS"

# Map run -> fastq_ftp
FTP="$WORK/fastq_ftp.tsv"
curl -s "$ENA?accession=PRJNA856121&result=read_run&fields=run_accession,fastq_ftp&format=tsv" -o "$FTP"

# Download SE fastq for each conjunctival run (resumable, parallel).
tail -n +2 "$SS" | cut -f1 | while read -r run; do
  url=$(awk -F'\t' -v r="$run" '$1==r{print $2}' "$FTP" | tr ';' '\n' | grep -E "/${run}\.fastq\.gz$|/${run}_1\.fastq\.gz$" | head -1)
  [ -z "$url" ] && url=$(awk -F'\t' -v r="$run" '$1==r{print $2}' "$FTP" | tr ';' '\n' | head -1)
  echo "$url"
done | sed 's#^#https://#' | grep -v '^https://$' \
  | xargs -P 6 -I{} bash -c 'u="$1"; b=$(basename "$u"); [ -s "'"$RAW"'/$b" ] || curl -s -o "'"$RAW"'/$b" "$u"' _ {}

echo "downloaded files: $(ls "$RAW" | wc -l)"
echo "FETCH856_OK"
