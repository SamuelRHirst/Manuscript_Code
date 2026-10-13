#!/bin/bash
# Weir-Cockerham Fst between four inland and four coastal individuals, by variant class.
# Genome-wide weighted Fst is reported in each vcftools log; windows are 500 kb, non-overlapping.

set -euo pipefail

VCF="pangenome/15vcf_biallelic/ingroup_biallelic.vcf.gz"
OUT="fst"
mkdir -p "$OUT"

printf "AH19K17\nAH21K03\nKF14K01\nKF15K02\n" > "$OUT/inland.txt"
printf "AH19K01\nAH19K05\nJF08K03\nJF08K05\n" > "$OUT/coastal.txt"

bcftools view -S <(cat "$OUT/inland.txt" "$OUT/coastal.txt") -t ^Chromosome_W "$VCF" \
    | bcftools +fill-tags -- -t AC,AN,AF \
    | bcftools view --min-ac 2:minor -Oz -o "$OUT/subset_mac2.vcf.gz"
tabix -p vcf "$OUT/subset_mac2.vcf.gz"

bcftools view --type snps,mnps "$OUT/subset_mac2.vcf.gz" -Oz -o "$OUT/snps.vcf.gz"

bcftools view --type indels "$OUT/subset_mac2.vcf.gz" \
    | awk -F'\t' '/^#/ || (length($4) < 50 && length($5) < 50)' | bgzip -c > "$OUT/indels.vcf.gz"

bcftools view --type indels "$OUT/subset_mac2.vcf.gz" \
    | awk -F'\t' '/^#/ || length($4) >= 50 || length($5) >= 50' | bgzip -c > "$OUT/svs.vcf.gz"

for T in snps indels svs; do
    vcftools --gzvcf "$OUT/$T.vcf.gz" \
        --weir-fst-pop "$OUT/inland.txt" \
        --weir-fst-pop "$OUT/coastal.txt" \
        --fst-window-size 500000 \
        --fst-window-step 500000 \
        --out "$OUT/${T}_fst_500kb"
done
