#!/bin/bash
# Per-copy methylation of RepeatMasker TE copies (non-TE repeat classes removed), split into full-length
# (overlapping an EDTA intact element) and fragmented copies; copies need >= 3 covered CpGs
# Usage: 04_te_copies.sh <sample_hap>

set -euo pipefail

H=$1
RM_OUT="repeats/repeatmasker/$H/${H}_purgedups.fa.out"
EDTA_GFF=$(ls repeats/edta/$H/*.mod.EDTA.intact.gff3)
CPG="methylation/haplotype/$H/$H.combined.bed.gz"
OUT="methylation/te_copies"
TMP="$OUT/tmp_$H"
mkdir -p "$TMP"

tail -n +4 "$RM_OUT" | awk -v OFS='\t' '
    $11 ~ /^(Simple_repeat|Low_complexity|Satellite|rRNA|tRNA|snRNA|scRNA|srpRNA|ARTEFACT)/ { next }
    { print $5, $6 - 1, $7, $2, $11 }' | sort -k1,1 -k2,2n > "$TMP/te_copies.bed"
awk '!/^#/ {print $1 "\t" ($4 - 1) "\t" $5}' "$EDTA_GFF" | sort -k1,1 -k2,2n | bedtools merge > "$TMP/edta.bed"
zcat "$CPG" | awk '!/^#/ {print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7}' | sort -k1,1 -k2,2n > "$TMP/cpg.bed"

bedtools intersect -u -a "$TMP/te_copies.bed" -b "$TMP/edta.bed" > "$TMP/fulllength.bed"
bedtools intersect -v -a "$TMP/te_copies.bed" -b "$TMP/edta.bed" > "$TMP/fragmented.bed"

for STATUS in fulllength fragmented; do
    sort -k1,1 -k2,2n "$TMP/$STATUS.bed" \
        | bedtools map -a - -b "$TMP/cpg.bed" -c 4,5,4 -o sum,sum,count -null 0 \
        | awk -v OFS='\t' -v h="$H" -v s="$STATUS" '$8 >= 3 {print h, s, $1, $2, $3, $4, $5, $6, $7, $8}' \
        > "$OUT/$H.$STATUS.tsv"
done

rm -rf "$TMP"
