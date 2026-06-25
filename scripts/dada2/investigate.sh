#!/usr/bin/env bash
# Copy DADA2 outputs into the repo, then diagnose the Snodgrassella question:
#  - is Snodgrassella present in the SILVA reference at all?
#  - how many ASVs are unclassified at genus level?
#  - which genera dominate (sanity)?
set -eu

WORK="${WORK:-$HOME/rosacea_dada2}"
OUT="$WORK/out"
REPO="/mnt/c/Users/jonat/OneDrive/Documents/research_portfolio_complete/rosacea_demodex"
DEST="$REPO/data/dada2_PRJNA692647"

echo "== copy outputs to repo =="
mkdir -p "$DEST"
cp -f "$OUT"/asv_table.tsv "$OUT"/taxonomy.tsv "$OUT"/track.tsv \
      "$OUT"/asv_seqs.fasta "$OUT"/E6_asv_validation.md \
      "$WORK"/sample_sheet.tsv "$HOME"/dada2_run.log "$DEST"/
ls -lh "$DEST"

echo "== Snodgrassella reference headers in SILVA train set =="
zcat "$WORK/ref/silva_nr99_v138.1_train_set.fa.gz" | grep -c 'Snodgrassella' || echo 0
echo "== Neisseriaceae references in SILVA (family of Snodgrassella) =="
zcat "$WORK/ref/silva_nr99_v138.1_train_set.fa.gz" | grep -c 'Neisseriaceae' || echo 0

echo "== ASVs with no genus assignment (col 7 = Genus) =="
awk -F'\t' 'NR>1 && ($7=="NA" || $7=="") {n++} END{print n" / "(NR-1)" ASVs unassigned at genus"}' "$OUT/taxonomy.tsv"

echo "== top 15 genera by ASV count =="
awk -F'\t' 'NR>1{g=$7; if(g==""||g=="NA")g="(none)"; c[g]++} END{for(k in c) print c[k]"\t"k}' "$OUT/taxonomy.tsv" | sort -rn | head -15

echo "== families present that could hide Snodgrassella (Neisseriaceae/Weeksellaceae) =="
awk -F'\t' 'NR>1{print $6}' "$OUT/taxonomy.tsv" | grep -iE 'neisser' | sort | uniq -c || echo 'none'

echo "INVESTIGATE_OK"
