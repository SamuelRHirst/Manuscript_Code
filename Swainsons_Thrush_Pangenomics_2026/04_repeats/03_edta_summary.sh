#!/bin/bash
# TE content per genome, recomputed from the EDTA GFFs
#   intact elements: intact.gff3 records without a Parent
#   masked %:        merged TEanno.gff3 bp / genome length
#   intact %:        merged intact.gff3 bp / genome length

set -euo pipefail

OUT="repeats/edta_summary.tsv"

merged_bp() {
    grep -v '^#' "$1" | awk -v OFS='\t' '$5 >= $4 {print $1, $4-1, $5}' \
        | sort -k1,1 -k2,2n | bedtools merge | awk '{s += $3-$2} END {print s+0}'
}

printf "genome\tintact_n\tintact_bp\tteanno_bp\tgenome_bp\tmasked_pct\tintact_pct\n" > "$OUT"

while IFS=$'\t' read -r NAME GENOME; do
    D="repeats/edta/$NAME"
    INTACT=$(ls "$D"/*.mod.EDTA.intact.gff3)
    TEANNO=$(ls "$D"/*.mod.EDTA.TEanno.gff3)

    n=$(grep -v '^#' "$INTACT" | grep -vc 'Parent=')
    ib=$(merged_bp "$INTACT")
    tb=$(merged_bp "$TEANNO")
    gl=$(awk '!/^>/ {s += length($0)} END {print s}' "$GENOME")

    awk -v g="$NAME" -v n="$n" -v ib="$ib" -v tb="$tb" -v gl="$gl" 'BEGIN {
        printf "%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\n", g, n, ib, tb, gl, 100*tb/gl, 100*ib/gl
    }' >> "$OUT"
done < "$(dirname "$0")/genomes.tsv"
