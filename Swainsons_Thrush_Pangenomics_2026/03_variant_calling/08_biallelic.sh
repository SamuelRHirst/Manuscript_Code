#!/bin/bash
# Split multiallelic sites for population-genetic analyses

set -euo pipefail

THREADS=8
REF="refs/inland_ref.pansn.fa"
OUT="pangenome/15vcf_biallelic"
mkdir -p "$OUT/tmp_sort"

bcftools norm -m-any -f "$REF" --threads "$THREADS" pangenome/11vcf_aa/ingroup_aa_gr.vcf.gz \
    | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING \
    | bcftools view -c 1 \
    | bcftools sort -T "$OUT/tmp_sort" -m 8G \
    | bgzip -c --threads "$THREADS" > "$OUT/ingroup_biallelic.vcf.gz"

tabix -p vcf "$OUT/ingroup_biallelic.vcf.gz"
rm -rf "$OUT/tmp_sort"
