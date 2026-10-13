#!/bin/bash
# CpG counts, coverage and modified-read counts per haplotype for:
#   baseline:  all CpGs; inside / outside RepeatMasker TEs; inside / outside EDTA intact elements
#   partition: full-length TE (EDTA intact), fragmented TE (RepeatMasker outside intact), non-TE
# Usage: 02_te_methylation.sh <sample_hap>   e.g. AH14K04_hap1

set -euo pipefail

H=$1
RM_GFF="repeats/repeatmasker/$H/${H}_purgedups.fa.out.gff"
EDTA_GFF=$(ls repeats/edta/$H/*.mod.EDTA.intact.gff3)
CPG="methylation/haplotype/$H/$H.combined.bed.gz"
OUT="methylation/te_classes"
TMP="$OUT/tmp_$H"
mkdir -p "$TMP"

awk '!/^#/ {print $1 "\t" ($4 - 1) "\t" $5}' "$RM_GFF" | sort -k1,1 -k2,2n | bedtools merge > "$TMP/rm.bed"
awk '!/^#/ && NF >= 5 {print $1 "\t" ($4 - 1) "\t" $5}' "$EDTA_GFF" | sort -k1,1 -k2,2n | bedtools merge > "$TMP/edta.bed"
zcat "$CPG" | awk '!/^#/ {print $1 "\t" $2 "\t" $3 "\t" $6 "\t" $7}' | sort -k1,1 -k2,2n > "$TMP/cpg.bed"

sum_cov_mod() { awk '{n++; cov += $4; mod += $5} END {printf "%d\t%d\t%d\n", n+0, cov+0, mod+0}'; }

TOTAL=$(sum_cov_mod < "$TMP/cpg.bed")
RM_TE=$(bedtools intersect -u -a "$TMP/cpg.bed" -b "$TMP/rm.bed" | sum_cov_mod)
RM_NONTE=$(bedtools intersect -v -a "$TMP/cpg.bed" -b "$TMP/rm.bed" | sum_cov_mod)
EDTA_FULL=$(bedtools intersect -u -a "$TMP/cpg.bed" -b "$TMP/edta.bed" | sum_cov_mod)
EDTA_NONFULL=$(bedtools intersect -v -a "$TMP/cpg.bed" -b "$TMP/edta.bed" | sum_cov_mod)
printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$H" "$TOTAL" "$RM_TE" "$RM_NONTE" "$EDTA_FULL" "$EDTA_NONFULL" > "$OUT/$H.baseline.tsv"

bedtools intersect -v -a "$TMP/cpg.bed" -b "$TMP/edta.bed" > "$TMP/not_full.bed"
FRAG=$(bedtools intersect -u -a "$TMP/not_full.bed" -b "$TMP/rm.bed" | sum_cov_mod)
NONTE=$(bedtools intersect -v -a "$TMP/not_full.bed" -b "$TMP/rm.bed" | sum_cov_mod)
printf "%s\t%s\t%s\t%s\t%s\n" "$H" "$TOTAL" "$EDTA_FULL" "$FRAG" "$NONTE" > "$OUT/$H.partition.tsv"

rm -rf "$TMP"
