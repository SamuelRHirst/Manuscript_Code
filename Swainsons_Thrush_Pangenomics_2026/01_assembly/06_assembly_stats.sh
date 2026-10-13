#!/bin/bash
# assembly_stats: https://github.com/MikeTrizna/assembly_stats

set -euo pipefail

SAMPLES=(AH14K04 AH19K01 AH19K05 AH19K17 AH21K03 JF08K03 JF08K05 KF14K01 KF15K02 MF12K02)

printf 'sample\thap\traw_bp\tpurged_bp\tremoved_pct\n'
for SAMPLE in "${SAMPLES[@]}"; do
    for HAP in hap1 hap2; do
        awk -v s="$SAMPLE" -v h="$HAP" '
            FNR==NR { if (!/^>/) r += length($0); next }
            !/^>/ { p += length($0) }
            END { printf "%s\t%s\t%d\t%d\t%.2f\n", s, h, r, p, 100*(r-p)/r }
        ' "assembly/$SAMPLE/${SAMPLE}_${HAP}.fa" \
          "assembly/$SAMPLE/purge_dups/$HAP/${SAMPLE}_${HAP}_purgedups.fa"
    done
done

for SAMPLE in "${SAMPLES[@]}"; do
    for HAP in hap1 hap2; do
        echo "$SAMPLE $HAP"
        assembly_stats "assembly/$SAMPLE/purge_dups/$HAP/${SAMPLE}_${HAP}_purgedups.fa"
    done
done
