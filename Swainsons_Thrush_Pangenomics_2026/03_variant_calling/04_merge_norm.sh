#!/bin/bash
# Concatenate chromosomes, merge multiallelic records and left-align against the reference

set -euo pipefail

THREADS=8
REF="refs/inland_ref.pansn.fa"

mkdir -p pangenome/10vcf_merged pangenome/tmp_sort

bcftools concat -a --threads "$THREADS" pangenome/09vcf_ploid/*.ploid.vcf.gz \
    | bcftools +fill-tags -- -t AN,AC,AF \
    | bcftools norm -m+both -f "$REF" --threads "$THREADS" \
    | bcftools sort -T pangenome/tmp_sort -m 8G \
    | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING \
    | bgzip -c --threads "$THREADS" > pangenome/10vcf_merged/merged_norm.vcf.gz

tabix -p vcf pangenome/10vcf_merged/merged_norm.vcf.gz
rm -rf pangenome/tmp_sort
