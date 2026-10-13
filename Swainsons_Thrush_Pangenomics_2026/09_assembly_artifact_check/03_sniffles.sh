#!/bin/bash
# Read-based SV calls with Sniffles2 from HiFi reads aligned to the inland reference (pbmm2)
# Usage: 03_sniffles.sh <sample>

set -euo pipefail

SAMPLE=$1
THREADS=16
OUT="artifact_check/sniffles"
mkdir -p "$OUT"

sniffles \
    --input "methylation/align/${SAMPLE}.inland_ref.bam" \
    --reference refs/inland_ref_NO_UNPLACED.fasta \
    --sample-id "$SAMPLE" \
    --snf "$OUT/$SAMPLE.snf" \
    --vcf "$OUT/$SAMPLE.vcf.gz" \
    --minsvlen 50 \
    --threads "$THREADS"
