#!/usr/bin/env bash
# Resumable, parallel per-file download of GEO GSE279452 (Ye 2026 sepsis scRNA-seq + VDJ).
# Files land in <GSM>_<sample>/ subfolders; finished files are checked against the GEO filelist size.
# Usage: bash download_GSE279452.sh [parallel_jobs=8] [time_budget_s=160]
set -u
cd "$(dirname "$0")"
P=${1:-6}; T=${2:-160}
awk -F'\t' 'NR>2 && $1=="File"{print $2"\t"$4}' filelist.txt > .files.tsv
fetch() {
  f="$1"; size="$2"; gsm="${f%%_*}"; rest="${f#*_}"; sample="${rest%%.*}"
  dir="${gsm}_${sample}"; mkdir -p "$dir"; out="$dir/$f"
  if [ -f "$out" ]; then
    cur=$(stat -c %s "$out"); [ "$cur" = "$size" ] && head -c 2 "$out" | od -An -tx1 | grep -q "1f 8b" && return 0
    # discard server error pages (non-gzip start) or oversize files before resuming
    if [ "$cur" -gt "$size" ] || ! head -c 2 "$out" | od -An -tx1 | grep -q "1f 8b"; then : > "$out"; fi
  fi
  url="https://ftp.ncbi.nlm.nih.gov/geo/samples/${gsm:0:7}nnn/${gsm}/suppl/${f}"
  curl -sS --fail -C - --retry 5 --retry-delay 5 -o "$out" "$url" 2>>curl_errors.log || true
}
export -f fetch
timeout "$T" xargs -a .files.tsv -P "$P" -L1 bash -c 'fetch "$0" "$1"' || true
done=0; total=0
while IFS=$'\t' read -r f size; do total=$((total+1)); gsm="${f%%_*}"; rest="${f#*_}"; s="${rest%%.*}"
  [ -f "${gsm}_${s}/$f" ] && [ "$(stat -c %s "${gsm}_${s}/$f")" = "$size" ] && done=$((done+1)); done < .files.tsv
echo "complete files: $done / $total"; du -sh --exclude='*.partial_unused' . | tail -1
