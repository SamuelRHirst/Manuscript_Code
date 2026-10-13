#!/bin/bash
# Usage: 02_edta.sh <genome name from genomes.tsv>

set -euo pipefail

NAME=$1
THREADS=16
GENOME=$(readlink -f "$(awk -v n="$NAME" '$1 == n {print $2}' "$(dirname "$0")/genomes.tsv")")
OUT="repeats/edta/$NAME"

# EDTA must start from an empty directory; reused directories silently lose TIR calls
rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"

EDTA.pl \
    --genome "$GENOME" \
    --species others \
    --step all \
    --threads "$THREADS" \
    --anno 1 \
    --force 1
