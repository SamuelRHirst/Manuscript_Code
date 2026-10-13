#!/bin/bash
# HiFi reads aligned to the bird's own purged pseudo-haplotype, then CpG methylation with pb-CpG-tools
# Usage: 01_align_cpg.sh <sample> <hap1|hap2>

set -euo pipefail

SAMPLE=$1
HAP=$2
THREADS=16
REF="assembly/$SAMPLE/purge_dups/$HAP/${SAMPLE}_${HAP}_purgedups.fa"
OUT="methylation/haplotype/${SAMPLE}_${HAP}"
mkdir -p "$OUT"

pbmm2 align --preset HIFI --sort --unmapped --sort-memory 3G -j "$THREADS" \
    "$REF" "reads/${SAMPLE}.hifi_reads.bam" "$OUT/${SAMPLE}_${HAP}.bam"

aligned_bam_to_cpg_scores \
    --bam "$OUT/${SAMPLE}_${HAP}.bam" \
    --ref "$REF" \
    --output-prefix "$OUT/${SAMPLE}_${HAP}" \
    --modsites-mode reference \
    --min-coverage 4 \
    --threads "$THREADS"
