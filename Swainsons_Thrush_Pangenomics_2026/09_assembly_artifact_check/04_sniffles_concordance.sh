#!/bin/bash
# Joint-genotype the ten Sniffles2 samples and compare with the pangenome SV genotypes

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
THREADS=16
OUT="artifact_check/sniffles"

for F in "$OUT"/*.snf; do printf '%s\t%s\n' "$F" "$(basename "$F" .snf)"; done > "$OUT/snf_list.tsv"

sniffles --input "$OUT/snf_list.tsv" --reference refs/inland_ref_NO_UNPLACED.fasta \
    --vcf "$OUT/merged.vcf.gz" --minsvlen 50 --threads "$THREADS"

{
    bcftools query -l "$OUT/merged.vcf.gz" | paste -sd'\t' | sed 's/^/#CHROM\tPOS\tSVTYPE\tSVLEN\tEND\t/'
    bcftools query -i 'FILTER="PASS"' -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END[\t%GT]\n' "$OUT/merged.vcf.gz"
} > "$OUT/merged_calls.tsv"

awk -F'\t' 'NR == 1 || $5 ~ /^SV/' pangenome/14vcf_aa_tab/variation_final.tab > "$OUT/pangenome_svs.tab"

python3 "$SCRIPTS/concordance.py" "$OUT/pangenome_svs.tab" "$OUT/merged_calls.tsv" "$OUT"
