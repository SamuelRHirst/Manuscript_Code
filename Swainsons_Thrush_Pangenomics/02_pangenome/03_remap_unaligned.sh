#!/bin/bash
# Remap contigs missed in step 02 with relaxed segment length and identity
# Usage: 03_remap_unaligned.sh <name>

set -euo pipefail

NAME=$1
THREADS=16
QUERY="01fastas_gz/${NAME}.pansn.fa.gz"

cd pangenome
mkdir -p 03remap

comm -23 \
    <(cut -f1 "${QUERY}.fai" | sort) \
    <(cut -f1 "02pafs/${NAME}.vs.inland.paf" | sort -u) \
    > "03remap/${NAME}.unaligned.txt"

if [ ! -s "03remap/${NAME}.unaligned.txt" ]; then
    touch "03remap/${NAME}.unaligned.vs.inland.paf"
    exit 0
fi

xargs samtools faidx "$QUERY" < "03remap/${NAME}.unaligned.txt" > "03remap/${NAME}.unaligned.fa"
samtools faidx "03remap/${NAME}.unaligned.fa"

wfmash 01fastas_gz/inland_ref.pansn.fa.gz "03remap/${NAME}.unaligned.fa" \
    -s 10000 \
    -p 85 \
    -m \
    -t "$THREADS" > "03remap/${NAME}.unaligned.vs.inland.paf"
