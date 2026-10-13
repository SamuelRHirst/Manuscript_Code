#!/bin/bash
# Lift each pseudo-haplotype's TE copies (with their methylation) to inland reference coordinates,
# through primary minimap2 alignments of at least 10 kb with MAPQ >= 5
# Usage: 07_hap_liftover.sh <sample> <hap1|hap2>

set -euo pipefail

SAMPLE=$1
HAP=$2
H="${SAMPLE}_${HAP}"
THREADS=16
ASM="assembly/$SAMPLE/purge_dups/$HAP/${H}_purgedups.fa"
TE="methylation/te_copies"
OUT="methylation/hap_discordance"
mkdir -p "$OUT/paf" "$OUT/lifted"

minimap2 -cx asm5 --cs -t "$THREADS" refs/inland_ref_NO_UNPLACED.fasta "$ASM" > "$OUT/paf/${H}_to_inland.paf"

cat "$TE/$H.fulllength.tsv" "$TE/$H.fragmented.tsv" \
    | awk -F'\t' -v OFS='\t' '{print $3, $4, $5, $2 "|" $7 "|" $6 "|" $8 "|" $9 "|" $10 "|" $3 ":" $4}' \
    | sort -k1,1 -k2,2n > "$OUT/$H.te.bed"

grep 'tp:A:P' "$OUT/paf/${H}_to_inland.paf" \
    | paftools.js liftover -l 10000 -q 5 - "$OUT/$H.te.bed" \
    | sort -k1,1 -k2,2n > "$OUT/lifted/$H.te_inland.bed"

rm "$OUT/$H.te.bed"
