#!/bin/bash
# Polarize with the C. minimus outgroup (INFO/AA), then drop the outgroup

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
VCF="pangenome/10vcf_merged/merged_norm.vcf.gz"
OUT="pangenome/11vcf_aa"
mkdir -p "$OUT"

bcftools query -s bCatMin1 -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$VCF" \
    | python3 "$SCRIPTS/ancestral_allele.py" \
    | bgzip -c > "$OUT/aa.tab.gz"
tabix -s1 -b2 -e2 "$OUT/aa.tab.gz"

echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele inferred from bCatMin1">' > "$OUT/aa_hdr.txt"

bcftools annotate -x INFO/CONFLICT,INFO/LV,INFO/PS,INFO/AT,INFO/TYPE,INFO/LEN,INFO/ORIGIN,INFO/NS "$VCF" \
    | bcftools annotate -a "$OUT/aa.tab.gz" -h "$OUT/aa_hdr.txt" -c CHROM,POS,REF,ALT,INFO/AA \
    | bcftools view -a -s ^bCatMin1 \
    | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING \
    | bcftools view -c 1 -Oz -o "$OUT/ingroup_aa.vcf.gz"

tabix -p vcf "$OUT/ingroup_aa.vcf.gz"
