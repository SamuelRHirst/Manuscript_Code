#!/bin/bash
# Usage: 05_busco.sh <sample> <hap1|hap2>

set -euo pipefail

SAMPLE=$1
HAP=$2
THREADS=12

busco \
    -i "assembly/$SAMPLE/purge_dups/$HAP/${SAMPLE}_${HAP}_purgedups.fa" \
    -l aves_odb10 \
    -m genome \
    -c "$THREADS" \
    -o "${SAMPLE}_${HAP}.aves_odb10.busco" \
    --out_path busco
