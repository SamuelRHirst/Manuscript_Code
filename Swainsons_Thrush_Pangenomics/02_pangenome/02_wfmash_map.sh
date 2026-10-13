#!/bin/bash
# Map each assembly to the inland reference, one alignment per contig (-N)
# Usage: 02_wfmash_map.sh <name>   e.g. AH14K04_hap1

set -euo pipefail

NAME=$1
THREADS=24

cd pangenome
mkdir -p 02pafs

wfmash 01fastas_gz/inland_ref.pansn.fa.gz "01fastas_gz/${NAME}.pansn.fa.gz" \
    -s 50000 \
    -p 90 \
    -m \
    -N \
    -t "$THREADS" > "02pafs/${NAME}.vs.inland.paf"
