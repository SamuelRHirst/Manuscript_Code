#!/bin/bash
# calcuts fell into two-peak mode for these three hap1 assemblies.
# Both haplotypes were mapped with the same reads, so hap2 cutoffs are reused.

set -euo pipefail

for SAMPLE in AH19K05 AH21K03 KF14K01; do
    ASM="$PWD/assembly/$SAMPLE/${SAMPLE}_hap1.fa"
    (
        cd "assembly/$SAMPLE/purge_dups/hap1"
        cp cutoffs cutoffs.calcuts
        cp ../hap2/cutoffs cutoffs

        purge_dups -2 -T cutoffs -c PB.base.cov asm.split.self.paf.gz > dups.bed
        get_seqs -e -p purged dups.bed "$ASM"

        mv purged.purged.fa "${SAMPLE}_hap1_purgedups.fa"
        mv purged.hap.fa "${SAMPLE}_hap1_removed.fa"
    )
done
