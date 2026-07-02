#!/usr/bin/env bash
# Build kraken2 + Bracken from source into $HOME/tools (no sudo required).
set -eu
TOOLS="$HOME/tools"; mkdir -p "$TOOLS"

if [ ! -x "$TOOLS/kraken2/kraken2" ]; then
  echo "== building kraken2 =="
  rm -rf "$TOOLS/kraken2-src"
  git clone --depth 1 https://github.com/DerrickWood/kraken2 "$TOOLS/kraken2-src"
  ( cd "$TOOLS/kraken2-src" && ./install_kraken2.sh "$TOOLS/kraken2" )
fi

if [ ! -x "$TOOLS/Bracken/bracken" ]; then
  echo "== building Bracken =="
  rm -rf "$TOOLS/Bracken"
  git clone --depth 1 https://github.com/jenniferlu717/Bracken "$TOOLS/Bracken"
  ( cd "$TOOLS/Bracken" && bash install_bracken.sh )
fi

echo "== versions =="
"$TOOLS/kraken2/kraken2" --version | head -1
ls -l "$TOOLS/Bracken/bracken"
echo BUILD_OK
