#!/bin/bash
# Usage: 01_hifiasm.sh <sample>

set -euo pipefail

SAMPLE=$1
THREADS=50
READS="reads/${SAMPLE}.hifi_reads.fastq.gz"
OUT="assembly/$SAMPLE"

mkdir -p "$OUT"
cd "$OUT"

hifiasm -o "${SAMPLE}.asm" -t "$THREADS" "../../$READS"

for HAP in hap1 hap2; do
    awk '$1=="S"{print ">"$2; print $3}' "${SAMPLE}.asm.bp.${HAP}.p_ctg.gfa" > "${SAMPLE}_${HAP}.fa"
done
