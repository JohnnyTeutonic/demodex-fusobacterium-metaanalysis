#!/usr/bin/env bash
# Build the PRJNA856121 (fcimb 2022 shotgun) sample sheet.
# Group label lives in the BioSample custom attribute <TAG>replicate</TAG>,
# e.g. "P_A4" (Patient = Demodex blepharitis) vs a control prefix (N_/C_/H_).
# Site lives in isolation_source (Conjunctival sac vs meibum).
#
# Output: $WORK/sample_sheet_full.tsv  (run  group  site  replicate)
#   group in {demodex, control, other}; site in {conjunctiva, meibum}
set -eu

WORK="${WORK:-$HOME/rosacea_856}"
mkdir -p "$WORK"
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"

RUNS="$WORK/runs.tsv"
curl -s "$ENA?accession=PRJNA856121&result=read_run&fields=run_accession,sample_accession,isolation_source&format=tsv" -o "$RUNS"
echo "runs fetched: $(($(wc -l < "$RUNS") - 1))"

SS="$WORK/sample_sheet_full.tsv"
printf 'run_accession\tgroup\tsite\treplicate\n' > "$SS"

tail -n +2 "$RUNS" | while IFS=$'\t' read -r run samp src; do
  xml=$(curl -s "https://www.ebi.ac.uk/ena/browser/api/xml/$samp")
  rep=$(printf '%s' "$xml" | grep -A1 '<TAG>replicate</TAG>' | grep '<VALUE>' | head -1 | sed -E 's#.*<VALUE>##; s#</VALUE>.*##')
  grp="other"
  case "$rep" in
    P_*|P*) grp="demodex" ;;
    N_*|N*|C_*|C*|H_*|H*) grp="control" ;;
  esac
  site="other"
  case "$src" in
    *onjunctiv*) site="conjunctiva" ;;
    *eibum*) site="meibum" ;;
  esac
  printf '%s\t%s\t%s\t%s\n' "$run" "$grp" "$site" "$rep" >> "$SS"
done

echo "== replicate prefixes seen =="
tail -n +2 "$SS" | cut -f4 | sed -E 's/[0-9].*$//' | sort | uniq -c
echo "== group x site counts =="
tail -n +2 "$SS" | awk -F'\t' '{c[$2"/"$3]++} END{for(k in c) print "  "k": "c[k]}'
echo "SHEET: $SS"
echo "META_OK"
