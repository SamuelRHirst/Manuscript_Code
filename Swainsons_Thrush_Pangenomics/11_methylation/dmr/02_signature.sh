#!/bin/bash
# DMRs between groups with methbat signature (default thresholds), for each comparison in comparisons.tsv

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
THREADS=8
OUT="methylation/dmr/signature"
mkdir -p "$OUT"

tail -n +2 "$SCRIPTS/comparisons.tsv" | while IFS=$'\t' read -r CMP BASELINE COMPARE; do
    {
        printf "identifier\tfilename\tlabels\n"
        awk -v OFS='\t' -v a="$BASELINE" -v b="$COMPARE" \
            '$2 == a || $2 == b {print $1, "methylation/dmr/pileup/" $1 ".5mC.bed.gz", $2}' "$SCRIPTS/samples.tsv"
    } > "$OUT/${CMP}_collection.tsv"

    methbat signature \
        --threads "$THREADS" \
        --baseline-category "$BASELINE" \
        --compare-category "$COMPARE" \
        --input-collection "$OUT/${CMP}_collection.tsv" \
        --output-prefix "$OUT/$CMP"
done
