#!/bin/bash
# Usage: 01_deconstruct.sh <chrom>   e.g. 01, Z, W
# vg v1.63.1; later versions do not resolve PanSN reference prefixes with -P

set -euo pipefail

CHR=$1
THREADS=24
NAME="Inland_ref#1#Chromosome_${CHR}"

cd pangenome
mkdir -p 07vcf_raw

GFA=$(ls 06pggb/"${NAME}".*.smooth.final.gfa)

# -P is a prefix, so for Z it picks up all three inland scaffolds
vg deconstruct \
    -P "$NAME" \
    -H '#' \
    -e -a \
    -t "$THREADS" \
    "$GFA" | bgzip -c > "07vcf_raw/${NAME}.vcf.gz"

tabix -p vcf "07vcf_raw/${NAME}.vcf.gz"
