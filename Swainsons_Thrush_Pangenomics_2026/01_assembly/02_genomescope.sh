#!/bin/bash
# k-mer genome size and heterozygosity from raw HiFi reads
# Usage: 02_genomescope.sh <sample>

set -euo pipefail

SAMPLE=$1
THREADS=16
K=21
READS="reads/${SAMPLE}.hifi_reads.fastq.gz"
OUT="assembly/$SAMPLE/genomescope"

mkdir -p "$OUT/kmc_tmp"
echo "$READS" > "$OUT/files.txt"

kmc -k${K} -t"$THREADS" -m100 -ci1 -cs100000 "@$OUT/files.txt" "$OUT/kmc_db" "$OUT/kmc_tmp"
kmc_tools transform "$OUT/kmc_db" histogram "$OUT/${SAMPLE}.histo" -cx100000
genomescope2 -i "$OUT/${SAMPLE}.histo" -o "$OUT" -k $K -p 2 -n "$SAMPLE"

rm -rf "$OUT"/kmc_db.kmc_pre "$OUT"/kmc_db.kmc_suf "$OUT/kmc_tmp"
