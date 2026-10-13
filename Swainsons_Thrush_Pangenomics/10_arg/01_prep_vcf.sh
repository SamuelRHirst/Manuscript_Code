#!/bin/bash
# Per-autosome SINGER input: the eight pure HiFi individuals, biallelic SNPs, <= 20% missing, PanSN prefix removed

set -euo pipefail

VCF="pangenome/15vcf_biallelic/ingroup_biallelic.vcf.gz"
SAMPLES="AH19K01,AH19K05,AH19K17,AH21K03,JF08K03,JF08K05,KF14K01,KF15K02"
OUT="arg/vcf"
mkdir -p "$OUT"

bcftools view -h "$VCF" | grep -oP '(?<=##contig=<ID=)[^,>]+' \
    | awk -v OFS='\t' '{new = $1; sub(/^Inland_ref#1#/, "", new); print $1, new}' > "$OUT/chr_rename.txt"

awk '$2 ~ /^Chromosome_[0-9]+$/' "$OUT/chr_rename.txt" | while read -r ORIG NEW; do
    bcftools view "$VCF" -s "$SAMPLES" -v snps -m2 -M2 -r "$ORIG" --min-ac 1 -e 'F_MISSING > 0.2' \
        | bcftools annotate --rename-chrs "$OUT/chr_rename.txt" -Oz -o "$OUT/$NEW.vcf.gz"
    tabix -p vcf "$OUT/$NEW.vcf.gz"
done
