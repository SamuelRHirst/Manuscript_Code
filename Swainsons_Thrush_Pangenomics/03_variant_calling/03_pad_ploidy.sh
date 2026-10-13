#!/bin/bash
# Give every chromosome VCF the same sample columns, then force diploid genotypes
# Usage: 03_pad_ploidy.sh <chrom>

set -euo pipefail

CHR=$1
NAME="Inland_ref#1#Chromosome_${CHR}"
SAMPLES="AH14K04 AH19K01 AH19K05 AH19K17 AH21K03 Coastal_ref JF08K03 JF08K05 KF14K01 KF15K02 MF12K02 bCatMin1"

cd pangenome
mkdir -p 09vcf_ploid

zcat "08vcf_norm/${NAME}.norm.vcf.gz" | awk -v samples="$SAMPLES" '
BEGIN { OFS = "\t"; n = split(samples, s, " ") }
/^##/ { print; next }
/^#CHROM/ {
    for (i = 10; i <= NF; i++) col[$i] = i
    NF = 9
    for (i = 1; i <= n; i++) $(9 + i) = s[i]
    print; next
}
{
    for (i = 1; i <= n; i++) out[i] = (s[i] in col) ? $(col[s[i]]) : "."
    NF = 9
    for (i = 1; i <= n; i++) $(9 + i) = out[i]
    print
}' | bcftools +fixploidy -Oz -o "09vcf_ploid/Chromosome_${CHR}.ploid.vcf.gz" -- -f 2

tabix -p vcf "09vcf_ploid/Chromosome_${CHR}.ploid.vcf.gz"
