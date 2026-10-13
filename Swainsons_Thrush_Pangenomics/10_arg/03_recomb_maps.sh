#!/bin/bash
# Per-chromosome recombination maps (100 kb bins) from the inland LDhat map of Justen et al. 2025,
# converted with rho = 4 Ne r at a fixed genome-wide Ne

set -euo pipefail

RHOMAP="refs/rhomap_res_win_100kb_chrnames.csv"
NE=200000
BIN=100000
OUT="arg/maps"
mkdir -p "$OUT"

sed 's/"//g' "$RHOMAP" | awk -F',' -v NE=$NE -v BIN=$BIN -v OUT="$OUT" 'NR > 1 {
    mid = $3 + 0; start = mid - BIN / 2; if (start < 0) start = 0
    printf "%d %d %.10g\n", start, mid + BIN / 2, ($2 + 0) / (4 * NE * 1000) > (OUT "/" $5 "_recomb.txt")
}'

for F in "$OUT"/Chromosome_*_recomb.txt; do sort -k1,1n -o "$F" "$F"; done
