#!/bin/bash
# Reference intervals of TE-annotated SVs from GraffiTE (chrom, 1-based position, REF length)

set -euo pipefail

OUT="methylation/dmr/te_enrichment"
mkdir -p "$OUT"

bcftools query -f '%CHROM\t%POS\t%REF\n' repeats/graffite/graffite_te_svs.vcf.gz \
    | awk -F'\t' -v OFS='\t' '{sub(/^Inland_ref#1#/, "", $1); print $1, $2, length($3)}' > "$OUT/graffite_raw.tsv"
