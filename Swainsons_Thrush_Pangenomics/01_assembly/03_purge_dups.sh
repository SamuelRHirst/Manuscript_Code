#!/bin/bash
# Usage: 03_purge_dups.sh <sample> <hap1|hap2>

set -euo pipefail

SAMPLE=$1
HAP=$2
THREADS=24

ASM="$PWD/assembly/$SAMPLE/${SAMPLE}_${HAP}.fa"
READS="$PWD/reads/${SAMPLE}.hifi_reads.fastq.gz"
OUT="assembly/$SAMPLE/purge_dups/$HAP"

mkdir -p "$OUT"
cd "$OUT"

minimap2 -x map-hifi -t "$THREADS" "$ASM" "$READS" | gzip -c > reads.paf.gz
pbcstat reads.paf.gz
calcuts PB.stat > cutoffs

split_fa "$ASM" > asm.split.fa
minimap2 -x asm5 -D -P -t "$THREADS" asm.split.fa asm.split.fa | gzip -c > asm.split.self.paf.gz

purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed
get_seqs -e -p purged dups.bed "$ASM"

mv purged.purged.fa "${SAMPLE}_${HAP}_purgedups.fa"
mv purged.hap.fa "${SAMPLE}_${HAP}_removed.fa"

rm -f asm.split.fa reads.paf.gz
