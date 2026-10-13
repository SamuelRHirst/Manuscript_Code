#!/bin/bash
# The inland Z is split across three scaffolds; build it as a single graph

set -euo pipefail

cd pangenome/05parts

cat "Inland_ref#1#Chromosome_Z_1.pan.fa.gz" \
    "Inland_ref#1#Chromosome_Z_2.pan.fa.gz" \
    "Inland_ref#1#Chromosome_Z_3.pan.fa.gz" > "Inland_ref#1#Chromosome_Z.pan.fa.gz"
samtools faidx "Inland_ref#1#Chromosome_Z.pan.fa.gz"

mkdir -p split_parts
mv "Inland_ref#1#Chromosome_Z_"[123].pan.fa.gz* split_parts/
