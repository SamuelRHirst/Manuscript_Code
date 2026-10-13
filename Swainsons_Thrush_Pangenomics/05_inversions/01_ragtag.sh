#!/bin/bash
# Scaffold an assembly onto the inland reference for SyRI
# Usage: 01_ragtag.sh <name> <assembly.fa>   e.g. AH14K04_hap1, coastal_ref

set -euo pipefail

NAME=$1
QUERY=$2
THREADS=16

ragtag.py scaffold \
    refs/inland_genome_folded_renamed.fasta \
    "$QUERY" \
    -t "$THREADS" \
    -o "inversions/ragtag/$NAME" \
    --mm2-params="-t $THREADS -x asm10"
