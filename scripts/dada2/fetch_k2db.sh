#!/usr/bin/env bash
# Download + extract the prebuilt Kraken2 Standard-8 DB (capped 8GB; loads in RAM).
# Includes Bracken kmer_distrib files for genus-level abundance re-estimation.
set -eu
DB="${DB:-$HOME/k2db}"
mkdir -p "$DB"
URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz"
TARB="$DB/k2_standard_08gb.tar.gz"
if [ ! -s "$DB/hash.k2d" ]; then
  echo "downloading DB ..."
  [ -s "$TARB" ] || curl -L -s -o "$TARB" "$URL"
  echo "download size: $(du -h "$TARB" | cut -f1)"
  echo "extracting ..."
  tar -xzf "$TARB" -C "$DB"
fi
echo "DB contents:"; ls -lh "$DB"
echo "K2DB_OK"
