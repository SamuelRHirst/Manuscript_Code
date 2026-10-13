#!/bin/bash
# Usage: 07_pggb.sh <chrom>   e.g. 01, Z, W (no chromosome 19 in this reference)

set -euo pipefail

CHR=$1
THREADS=25

cd pangenome
mkdir -p 06pggb

# concurrent runs sharing a temp dir collide on smoothxg mmap files
TMP="temp/Chromosome_${CHR}"
mkdir -p "$TMP"

pggb \
    -i "05parts/Inland_ref#1#Chromosome_${CHR}.pan.fa.gz" \
    -o 06pggb \
    -n 23 \
    -p 90 \
    -s 50000 \
    -k 79 \
    -t "$THREADS" \
    -T "$THREADS" \
    -D "$TMP"

rm -rf "$TMP"
