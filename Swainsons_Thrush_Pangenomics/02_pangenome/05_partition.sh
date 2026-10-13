#!/bin/bash
# One FASTA per reference chromosome: the reference sequence plus every contig assigned to it
# Usage: 05_partition.sh <chrom>   e.g. Inland_ref#1#Chromosome_01 (see chrom_list.txt)

set -euo pipefail

CHR=$1
THREADS=4

cd pangenome
mkdir -p 05parts

LIST="05parts/${CHR}.contigs.txt"
OUT="05parts/${CHR}.pan.fa"

awk -v T="$CHR" '$6 == T {print $1}' 02pafs/*.vs.inland.paf 04partitioning/rescues.paf \
    | grep -vi mitochondria | sort -u > "$LIST"

samtools faidx 01fastas_gz/inland_ref.pansn.fa.gz "$CHR" > "$OUT"

for FASTA in 01fastas_gz/*.pansn.fa.gz; do
    [[ "$FASTA" == *inland_ref* ]] && continue
    comm -12 <(cut -f1 "${FASTA}.fai" | sort) "$LIST" | while read -r CONTIG; do
        samtools faidx "$FASTA" "$CONTIG"
    done >> "$OUT"
done

bgzip -f -@ "$THREADS" "$OUT"
samtools faidx "${OUT}.gz"
rm "$LIST"
