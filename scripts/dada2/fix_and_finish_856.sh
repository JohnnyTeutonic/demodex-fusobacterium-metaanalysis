#!/usr/bin/env bash
# Fix Bracken's hard-coded `python` (this system only has python3), then
# (re)run bracken for any reports missing bracken output, aggregate, and
# run the 3-cohort meta-analysis. Skips the slow integrity + kraken2 steps
# (kraken2 reports are already produced).
set -eu
D="$(dirname "$0")"
WORK="$HOME/rosacea_856"
DB="$HOME/k2db"
BRACKEN="$HOME/tools/Bracken/bracken"

# 1. patch bracken wrapper: exec line starts with spaces then 'python '
sed -i 's/^\( *\)python /\1python3 /' "$BRACKEN"
grep -n 'python3 .*est_abundance' "$BRACKEN" || true

# 2. run bracken for every kraken report lacking a bracken output
for rep in "$WORK"/kraken/reports/*.kreport; do
  run="$(basename "$rep" .kreport)"
  brk="$WORK/kraken/bracken/$run.bracken"
  if [ ! -s "$brk" ]; then
    "$BRACKEN" -d "$DB" -i "$rep" -o "$brk" -r 75 -l G -t 1 \
      > "$WORK/kraken/bracken/$run.blog" 2>&1 || echo "bracken-warn $run"
  fi
done
echo "reports: $(ls "$WORK"/kraken/reports/*.kreport 2>/dev/null | wc -l), bracken: $(ls "$WORK"/kraken/bracken/*.bracken 2>/dev/null | wc -l)"

# 3. aggregate + 4. meta-analysis
python3 "$D/agg_bracken.py" --work "$WORK" --db "$DB"
python3 "$D/cross_cohort_concordance.py" \
  --cohort "$HOME/rosacea_dada2/out:$HOME/rosacea_dada2/sample_sheet.tsv:discovery" \
  --cohort "$HOME/rosacea_657256/out:$HOME/rosacea_657256/sample_sheet.tsv:rep2020" \
  --cohort "$WORK/out:$WORK/sample_sheet.tsv:shotgun2022" \
  --min-prev 0.5 --out "$HOME/rosacea_meta3.md"
echo FINISH856_OK
