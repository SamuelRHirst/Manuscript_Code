#!/bin/bash
# PanSN naming (sample#haplotype#contig), bgzip and index for all 23 assemblies

set -euo pipefail

THREADS=8
OUT="pangenome/01fastas_gz"
mkdir -p "$OUT"

pansn() {
    local fasta=$1 sample=$2 hap=$3 dest=$4
    awk -v s="$sample" -v h="$hap" '/^>/{split(substr($0,2),a," "); print ">"s"#"h"#"a[1]; next}{print}' "$fasta" \
        | bgzip -c -@ "$THREADS" > "$dest"
    samtools faidx "$dest"
}

for S in AH14K04 AH19K01 AH19K05 AH19K17 AH21K03 JF08K03 JF08K05 KF14K01 KF15K02 MF12K02; do
    for H in 1 2; do
        pansn "assembly/$S/purge_dups/hap$H/${S}_hap${H}_purgedups.fa" "$S" "$H" "$OUT/${S}_hap${H}.pansn.fa.gz"
    done
done

pansn refs/inland_genome_folded_renamed.fasta Inland_ref 1 "$OUT/inland_ref.pansn.fa.gz"
pansn <(sed 's/_coastal$//' refs/coastal_genome_folded_renamed.fasta) Coastal_ref 1 "$OUT/coastal_ref.pansn.fa.gz"
pansn refs/GCA_052095055.2_bCatMin1_p1.2_genomic.fasta bCatMin1 1 "$OUT/bCatMin1_p.pansn.fa.gz"
