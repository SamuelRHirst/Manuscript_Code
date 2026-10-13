#!/bin/bash

set -euo pipefail

THREADS=16
OUT="pangenome/06pggb/stats"
SUMMARY="$OUT/graph_stats_summary.tsv"
mkdir -p "$OUT"

printf "chrom\tlength\tnodes\tedges\tpaths\tsteps\tmean_node_bp\tedge_node_ratio\n" > "$SUMMARY"

for CHR in $(printf "%02d\n" {1..18} {20..40}) W Z; do
    NAME="Chromosome_${CHR}"
    GRAPH=$(ls pangenome/06pggb/"Inland_ref#1#${NAME}".*.smooth.final.og)
    odgi stats -i "$GRAPH" -S -t "$THREADS" > "$OUT/${NAME}_stats.txt"
    awk -v c="$NAME" 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=5 {
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%.4f\n", c, $1, $2, $3, $4, $5, $1/$2, $3/$2
    }' "$OUT/${NAME}_stats.txt" >> "$SUMMARY"
done

# W is excluded from genome-wide totals
awk 'BEGIN{FS="\t"} NR>1 && $1 != "Chromosome_W" {len+=$2; nod+=$3; edg+=$4; pth+=$5}
     END{printf "length %.1f Mb\nnodes %d\nedges %d\npaths %d\n", len/1e6, nod, edg, pth}' "$SUMMARY"
