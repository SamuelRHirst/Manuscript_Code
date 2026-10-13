#!/bin/bash
# Usage: 04_svim_asm.sh <bird>   e.g. AH14K04, coastal_ref (haploid)

set -euo pipefail

BIRD=$1
THREADS=15
REF="pangenome/01fastas_gz/inland_ref.pansn.fa.gz"
OUT="inversions/svim"
mkdir -p "$OUT/bams"

align() {
    minimap2 -a --eqx -x asm5 -t "$THREADS" --cs "$REF" "pangenome/01fastas_gz/$1.pansn.fa.gz" \
        | samtools sort -m4G -@ "$THREADS" -o "$OUT/bams/$1.bam"
    samtools index "$OUT/bams/$1.bam"
}

if [ "$BIRD" == "coastal_ref" ]; then
    align coastal_ref
    svim-asm haploid "$OUT/$BIRD" "$OUT/bams/coastal_ref.bam" "$REF"
else
    align "${BIRD}_hap1"
    align "${BIRD}_hap2"
    svim-asm diploid "$OUT/$BIRD" "$OUT/bams/${BIRD}_hap1.bam" "$OUT/bams/${BIRD}_hap2.bam" "$REF"
fi
