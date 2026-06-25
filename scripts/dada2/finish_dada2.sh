#!/usr/bin/env bash
# Wait for the DADA2 run to finish, then run the ASV analysis and copy all
# outputs into the repo (so results are version-controlled + manuscript-ready).
set -eu

WORK="${WORK:-$HOME/rosacea_dada2}"
OUT="$WORK/out"
REPO="/mnt/c/Users/jonat/OneDrive/Documents/research_portfolio_complete/rosacea_demodex"
DEST="$REPO/data/dada2_PRJNA692647"
ANALYZE="$REPO/scripts/dada2/analyze_asv.py"
TIMEOUT="${TIMEOUT:-5400}"   # 90 min cap

echo "waiting for $OUT/taxonomy.tsv ..."
elapsed=0
while [ ! -s "$OUT/taxonomy.tsv" ]; do
  if grep -q 'DADA2_OK' "$HOME/dada2_run.log" 2>/dev/null && [ -s "$OUT/asv_table.tsv" ]; then break; fi
  if grep -qiE 'Error|Execution halted' "$HOME/dada2_run.log" 2>/dev/null; then
    echo "ERROR detected in dada2_run.log:"; tail -n 15 "$HOME/dada2_run.log"; exit 1
  fi
  sleep 15; elapsed=$((elapsed+15))
  if [ "$elapsed" -ge "$TIMEOUT" ]; then echo "TIMEOUT after ${TIMEOUT}s"; exit 2; fi
done
echo "DADA2 outputs present after ~${elapsed}s."

echo "== running ASV analysis =="
python3 "$ANALYZE" --out-dir "$OUT" --sample-sheet "$WORK/sample_sheet.tsv"

echo "== copying outputs into repo: $DEST =="
mkdir -p "$DEST"
cp -f "$OUT"/asv_table.tsv "$OUT"/taxonomy.tsv "$OUT"/track.tsv \
      "$OUT"/asv_seqs.fasta "$OUT"/E6_asv_validation.md \
      "$WORK"/sample_sheet.tsv "$DEST"/ 2>/dev/null || true
cp -f "$HOME"/dada2_run.log "$DEST"/dada2_run.log 2>/dev/null || true
ls -lh "$DEST"
echo "ALL_DONE"
