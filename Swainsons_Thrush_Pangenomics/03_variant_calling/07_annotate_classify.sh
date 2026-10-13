#!/bin/bash
# Add genomic region (INFO/GR), export to a table and classify variants

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
VCF="pangenome/11vcf_aa/ingroup_aa.vcf.gz"
VCF_GR="pangenome/11vcf_aa/ingroup_aa_gr.vcf.gz"
OUT="pangenome/14vcf_aa_tab"
mkdir -p "$OUT"

SAMPLES=$(cut -f1 "$SCRIPTS/sample_populations.txt" | paste -sd,)

echo '##INFO=<ID=GR,Number=1,Type=String,Description="Genomic region: exon, intron, promoter or intergenic">' > "$OUT/gr_hdr.txt"

bcftools annotate \
    -a pangenome/13genomic_overlaps/genomic_overlaps.bed.gz \
    -h "$OUT/gr_hdr.txt" \
    -c CHROM,FROM,TO,INFO/GR \
    "$VCF" -Oz -o "$VCF_GR"
tabix -p vcf "$VCF_GR"

bcftools query -s "$SAMPLES" \
    -f '%CHROM\t%POS0\t%END\t%TYPE\t%REF\t%ALT\t%INFO/AA\t%INFO/GR[\t%GT]\n' \
    "$VCF_GR" > "$OUT/pggb_variation.tab"

python3 "$SCRIPTS/classify_variants.py" \
    "$OUT/pggb_variation.tab" \
    "$SCRIPTS/sample_populations.txt" \
    "$OUT/variation_final.tab"
