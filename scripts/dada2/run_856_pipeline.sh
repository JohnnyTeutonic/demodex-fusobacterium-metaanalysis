#!/usr/bin/env bash
# Master pipeline for the shotgun cohort (PRJNA856121), then 3-cohort meta-analysis.
#   1. gzip integrity check (catches truncated/partial downloads)
#   2. Kraken2 + Bracken genus classification
#   3. aggregate -> asv_table/taxonomy
#   4. cross-cohort meta-analysis (discovery + rep2020 + shotgun)
set -eu
D="$(dirname "$0")"
WORK="$HOME/rosacea_856"

echo "== 1. integrity check =="
bad=0
for f in "$WORK"/raw/*.fastq.gz; do
  if ! gzip -t "$f" 2>/dev/null; then echo "CORRUPT: $f"; rm -f "$f"; bad=$((bad+1)); fi
done
echo "corrupt (removed): $bad"
if [ "$bad" -gt 0 ]; then
  echo "re-fetching missing files ..."
  bash "$D/fetch_856_conj.sh"
  for f in "$WORK"/raw/*.fastq.gz; do gzip -t "$f" 2>/dev/null || { echo "STILL CORRUPT: $f"; exit 1; }; done
fi

echo "== 2. classify (kraken2 + bracken) =="
bash "$D/classify_856.sh"

echo "== 3. aggregate bracken -> genus table =="
python3 "$D/agg_bracken.py" --work "$WORK" --db "$HOME/k2db"

echo "== 4. three-cohort meta-analysis =="
python3 "$D/cross_cohort_concordance.py" \
  --cohort "$HOME/rosacea_dada2/out:$HOME/rosacea_dada2/sample_sheet.tsv:discovery" \
  --cohort "$HOME/rosacea_657256/out:$HOME/rosacea_657256/sample_sheet.tsv:rep2020" \
  --cohort "$WORK/out:$WORK/sample_sheet.tsv:shotgun2022" \
  --min-prev 0.5 --out "$HOME/rosacea_meta3.md"

echo PIPELINE_OK
