#!/bin/bash
# SyRI needs matching chromosome names on both sides: keep placed chromosomes and drop the RagTag suffix.
# HiFi birds are male, so W is removed from their reference. The coastal reference lacks several
# microchromosomes, and SyRI aborts if any chromosome has no syntenic region, so those are removed too.

set -euo pipefail

REF="refs/inland_ref_NO_UNPLACED.fasta"
OUT="inversions/syri/fixed"
mkdir -p "$OUT"
samtools faidx "$REF"

for D in inversions/ragtag/*/; do
    NAME=$(basename "$D")
    awk '/^>/ { keep = /^>Chromosome_[A-Za-z0-9_]+_RagTag$/; if (keep) { sub(/_RagTag$/, ""); print } next } keep' \
        "$D/ragtag.scaffold.fasta" > "$OUT/${NAME}.fasta"
done

samtools faidx "$REF" $(cut -f1 "$REF.fai" | grep -v '^Chromosome_W$') > "$OUT/inland_ref_noW.fasta"

DROP='^Chromosome_(W|32|33|34|35|37|38|40)$'
SHARED=$(comm -12 <(cut -f1 "$REF.fai" | sort) <(grep '^>' "$OUT/coastal_ref.fasta" | tr -d '>' | sort) | grep -vE "$DROP")
samtools faidx "$REF" $SHARED > "$OUT/inland_ref_for_coastal.fasta"
