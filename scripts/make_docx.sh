#!/usr/bin/env bash
# Convert the mSystems manuscript LaTeX to .docx for ASM submission.
# Prefers pandoc (faithful structure); reports tool availability.
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1
SRC="demodex_surrogate_16S.tex"
OUT="demodex_surrogate_16S.docx"

echo "== tool check =="
if command -v pandoc >/dev/null 2>&1; then
  echo "pandoc: $(pandoc --version | head -1)"
  HAVE_PANDOC=1
else
  echo "pandoc: MISSING"
  HAVE_PANDOC=0
fi
python3 -c 'import docx; print("python-docx:", docx.__version__)' 2>/dev/null \
  || echo "python-docx: MISSING"

if [ "$HAVE_PANDOC" = "1" ]; then
  echo "== converting with pandoc =="
  pandoc "$SRC" -o "$OUT" \
    --from=latex --to=docx \
    --resource-path=.:results 2>pandoc.log \
    && echo "WROTE $OUT" \
    || { echo "pandoc FAILED; see pandoc.log"; tail -20 pandoc.log; }
fi
