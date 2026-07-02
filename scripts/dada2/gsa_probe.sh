#!/usr/bin/env bash
# Probe GSA/NGDC endpoints for CRA014100 (P&V 2024 eyelash cohort).
# Goal: find (a) a directory listing exposing CRR run accessions + fastq URLs,
# and (b) metadata mapping runs -> group (Demodex vs control).
# Logs reachability + first bytes of each endpoint. Nothing large is downloaded.
set -u
OUT="$HOME/rosacea_gsa"; mkdir -p "$OUT"
CRA="CRA014100"
T=60

probe() {
  local name="$1" url="$2"
  local f="$OUT/$name"
  local code
  code=$(curl -sL --max-time "$T" -o "$f" -w '%{http_code}' "$url" 2>"$OUT/$name.err" || echo "ERR")
  printf '%-28s http=%s size=%s\n' "$name" "$code" "$(wc -c < "$f" 2>/dev/null || echo 0)"
  head -c 200 "$f" 2>/dev/null | tr '\n' ' '; echo; echo "----"
}

echo "== listing endpoints =="
probe "https_gsa"   "https://download.cncb.ac.cn/gsa/$CRA/"
probe "https_gsa2"  "https://download.cncb.ac.cn/gsa2/$CRA/"
probe "https_gsa3"  "https://download.cncb.ac.cn/gsa3/$CRA/"
probe "ftp_big"     "ftp://download.big.ac.cn/gsa/$CRA/"
probe "ftp_big2"    "ftp://download.big.ac.cn/gsa2/$CRA/"

echo "== metadata endpoints =="
probe "browse_html" "https://ngdc.cncb.ac.cn/gsa/browse/$CRA"
probe "meta_json"   "https://ngdc.cncb.ac.cn/gsa/search/data?searchTerm=$CRA"

echo GSA_PROBE_OK
