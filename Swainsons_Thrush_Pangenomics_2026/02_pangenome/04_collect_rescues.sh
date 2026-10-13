#!/bin/bash
# Keep the highest-scoring rescue alignment per contig

set -euo pipefail

cd pangenome
mkdir -p 04partitioning

cat 03remap/*.unaligned.vs.inland.paf \
    | awk -v OFS='\t' '{print $1, $11, $0}' \
    | sort -k1,1 -k2,2nr \
    | awk '$1 != last {print; last = $1}' \
    | cut -f3- > 04partitioning/rescues.paf
