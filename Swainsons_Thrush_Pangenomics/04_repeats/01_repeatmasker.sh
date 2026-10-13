#!/bin/bash
# RepeatMasker with the curated Catharus library (Termignoni-Garcia et al. 2026)
# Usage: 01_repeatmasker.sh <genome name from genomes.tsv>

set -euo pipefail

NAME=$1
THREADS=32
LIB="refs/Catharus_TE_library.fasta"
GENOME=$(awk -v n="$NAME" '$1 == n {print $2}' "$(dirname "$0")/genomes.tsv")
OUT="repeats/repeatmasker/$NAME"

mkdir -p "$OUT"

RepeatMasker -pa "$THREADS" -lib "$LIB" -gff -xsmall -dir "$OUT" "$GENOME"
