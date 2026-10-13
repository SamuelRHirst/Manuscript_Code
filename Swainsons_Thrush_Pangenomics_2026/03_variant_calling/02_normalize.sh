#!/bin/bash
# Flatten nested bubbles (vcfbub) and realign complex alleles (vcfwave)
# Usage: 02_normalize.sh <chrom>

set -euo pipefail

CHR=$1
THREADS=16
NAME="Inland_ref#1#Chromosome_${CHR}"

cd pangenome
mkdir -p 08vcf_norm

vcfbub --input "07vcf_raw/${NAME}.vcf.gz" -l 0 -a 100000 \
    | vcfwave -I 1000 -t "$THREADS" \
    | bcftools sort -Oz -o "08vcf_norm/${NAME}.norm.vcf.gz"

tabix -p vcf "08vcf_norm/${NAME}.norm.vcf.gz"
