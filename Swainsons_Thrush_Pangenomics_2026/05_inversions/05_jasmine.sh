#!/bin/bash
# Merge per-sample calls across birds with Jasmine, separately for SyRI and SVIM-asm

set -euo pipefail

BIRDS=(AH14K04 AH19K01 AH19K05 AH19K17 AH21K03 coastal_ref JF08K03 JF08K05 KF14K01 KF15K02 MF12K02)
REF="refs/inland_ref_NO_UNPLACED.fasta"
OPTS="--ignore_strand --mutual_distance --max_dist_linear=0.1 --min_dist=50 --use_end --output_genotypes --normalize_type"

# SyRI: inversion records only, one column per haplotype named after its bird
D="inversions/syri/inv"
mkdir -p "$D"
: > "$D/vcf_list.txt"
for B in "${BIRDS[@]}"; do
    for VCF in inversions/syri/results/${B}_*syri.vcf; do
        awk -F'\t' '/^#/ || $5 == "<INV>"' "$VCF" | bcftools reheader -s <(echo "$B") > "$D/$(basename "$VCF")"
        echo "$D/$(basename "$VCF")" >> "$D/vcf_list.txt"
    done
done

jasmine file_list="$D/vcf_list.txt" out_file=inversions/syri/SyRI_combined.vcf \
    genome_file="$REF" $OPTS --allow_intrasample --clique_merging threads=8

# SVIM-asm: one diploid VCF per bird, in BIRDS order
printf "%s\n" "${BIRDS[@]/#/inversions/svim/}" | sed 's|$|/variants.vcf|' > inversions/svim/vcf_list.txt

jasmine file_list=inversions/svim/vcf_list.txt out_file=inversions/svim/SVIM_combined.vcf \
    genome_file=pangenome/01fastas_gz/inland_ref.pansn.fa.gz $OPTS threads=15
