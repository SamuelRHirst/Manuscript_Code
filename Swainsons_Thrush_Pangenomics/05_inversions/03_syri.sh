#!/bin/bash
# Usage: 03_syri.sh <name>   e.g. AH14K04_hap1, coastal_ref

set -euo pipefail

NAME=$1
THREADS=15
FIX="inversions/syri/fixed"
OUT="inversions/syri/results"
mkdir -p "$OUT"

if [ "$NAME" == "coastal_ref" ]; then
    REF="$FIX/inland_ref_for_coastal.fasta"
else
    REF="$FIX/inland_ref_noW.fasta"
fi
QUERY="$FIX/${NAME}.fasta"
BAM="$OUT/${NAME}.bam"

minimap2 -a --eqx -x asm5 --cs -t "$THREADS" "$REF" "$QUERY" | samtools sort -@ 4 -o "$BAM"
samtools index "$BAM"

syri -c "$BAM" -r "$REF" -q "$QUERY" -F B \
    --prefix "${NAME}_" --dir "$OUT" --nc 12 --cigar --no-chrmatch
