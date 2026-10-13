#!/bin/bash
# Assign each variant one feature, by priority CDS > 5'UTR > 3'UTR > intron > intergenic

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
VARS="pangenome/14vcf_aa_tab/variation_final.tab"
OUT="sfs"
BEDS="$OUT/beds"
mkdir -p "$BEDS"

python3 "$SCRIPTS/parse_gff.py" refs/inland_ref.pansn.gff "$BEDS"

FEATURES="cds five_utr three_utr intron"
for F in $FEATURES; do
    bedtools sort -i "$BEDS/$F.bed" | bedtools merge > "$BEDS/${F}_merged.bed"
done

# variant line numbers, so hits can be matched back to the table
awk -v OFS='\t' 'NR > 1 {print $1, $2, $3, NR}' "$VARS" | sort -k1,1 -k2,2n > "$BEDS/variants.bed"

for F in $FEATURES; do
    bedtools intersect -a "$BEDS/variants.bed" -b "$BEDS/${F}_merged.bed" -u | cut -f4 > "$BEDS/hit_$F.txt"
done

awk -v OFS='\t' -v beds="$BEDS" '
BEGIN {
    n = split("cds five_utr three_utr intron", f, " ")
    split("CDS|5'"'"'UTR|3'"'"'UTR|Intron", label, "|")
    for (i = 1; i <= n; i++)
        while ((getline line < (beds "/hit_" f[i] ".txt")) > 0)
            if (!(line in feat)) feat[line] = label[i]
}
NR == 1 { print $0, "feature"; next }
{ print $0, (NR in feat) ? feat[NR] : "Intergenic" }
' "$VARS" > "$OUT/variants_annotated.tsv"
