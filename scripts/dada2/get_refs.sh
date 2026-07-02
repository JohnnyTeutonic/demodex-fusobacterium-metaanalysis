#!/usr/bin/env bash
# Fetch authoritative Vancouver-style citations via DOI content negotiation
# (crossref/datacite), so manuscript references are not hallucinated.
set -u
dois=(
  "10.1007/s40123-021-00356-z"     # Liang 2021 discovery cohort
  "10.3389/fmed.2020.592759"       # Yan 2020 rep cohort
  "10.3389/fcimb.2022.922753"      # shotgun 2022 cohort
  "10.1186/s13071-024-06122-x"     # Zhang/Zou 2024 corroboration cohort
  "10.1186/s13059-019-1891-0"      # Kraken2 (Wood 2019)
  "10.7717/peerj-cs.104"           # Bracken (Lu 2017)
  "10.1111/j.1420-9101.2005.00917.x" # Whitlock 2005 weighted Z meta
  "10.1111/j.2517-6161.1995.tb02031.x" # Benjamini-Hochberg 1995 FDR
  "10.1093/nar/gks1219"            # SILVA (Quast 2013)
  "10.1038/nmeth.3869"             # DADA2 (Callahan 2016)
  "10.14806/ej.17.1.200"           # cutadapt (Martin 2011)
  "10.1038/ismej.2017.119"         # Callahan 2017 ASV>OTU
  "10.1016/j.jaad.2017.03.040"     # Chang 2017 Demodex-rosacea meta
)
for d in "${dois[@]}"; do
  echo "### $d"
  curl -sL --max-time 40 -H "Accept: text/x-bibliography; style=vancouver" "https://doi.org/$d"
  echo; echo
done
echo REFS_DONE
