#!/usr/bin/env bash
# WS2 recon: what tools are available, and what do PRJEB13411 D. folliculorum
# WGA reads look like (length, layout)? Downloads the smallest run only.
set -u
echo "=== available tools ==="
for t in minimap2 bwa bwa-mem2 bowtie2 kraken2 samtools seqkit vsearch blastn makeblastdb Rscript python3 cutadapt; do
  p=$(command -v "$t" 2>/dev/null || true)
  printf '  %-12s %s\n' "$t" "${p:-MISSING}"
done

WORK="${WORK:-$HOME/demodex}"
RAW="$WORK/raw"; mkdir -p "$RAW"; cd "$RAW"
RUN="ERR2338205"   # 249 MB, smallest-ish
URL="https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR233/005/$RUN/$RUN.fastq.gz"
echo "=== downloading $RUN (recon) ==="
[ -s "$RUN.fastq.gz" ] || curl -s -o "$RUN.fastq.gz" "$URL"
echo "  bytes: $(stat -c %s "$RUN.fastq.gz" 2>/dev/null)"
echo "=== n reads ==="
expr $(zcat "$RUN.fastq.gz" | head -400000 | wc -l) / 4
echo "  (from first 100k reads sampled)"
echo "=== read-length distribution (100k sample, 50 bp bins) ==="
zcat "$RUN.fastq.gz" | head -400000 | sed -n '2~4p' | awk '{print int(length/50)*50}' | sort -n | uniq -c
echo "=== first 3 read headers ==="
zcat "$RUN.fastq.gz" | sed -n '1~4p' | head -3
echo "RECON_TOOLS_OK"
