#!/bin/bash
# HiFi reads aligned to the inland reference, then per-individual 5mC pileup with MethBat v1.1.0
# Usage: 01_align_pileup.sh <sample>

set -euo pipefail

SAMPLE=$1
THREADS=16
REF="refs/inland_ref_NO_UNPLACED.fasta"
OUT="methylation/dmr"
mkdir -p "$OUT/align" "$OUT/pileup"

pbmm2 align --preset HIFI --sort -j "$THREADS" --sort-memory 3G \
    "$REF" "reads/${SAMPLE}.hifi_reads.bam" "$OUT/align/${SAMPLE}.inland_ref.bam"
samtools index "$OUT/align/${SAMPLE}.inland_ref.bam"

methbat pileup \
    --threads "$THREADS" \
    --input-bam "$OUT/align/${SAMPLE}.inland_ref.bam" \
    --output-prefix "$OUT/pileup/$SAMPLE" \
    --skip-5hmC \
    --skip-6mA
