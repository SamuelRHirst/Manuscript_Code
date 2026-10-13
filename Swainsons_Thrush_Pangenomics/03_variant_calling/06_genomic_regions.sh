#!/bin/bash
# Non-overlapping genomic regions from the inland annotation: exon > intron > promoter (2 kb upstream) > intergenic

set -euo pipefail

GFF="refs/inland_ref.pansn.gff"
FAI="refs/inland_ref.pansn.fa.fai"
OUT="pangenome/13genomic_overlaps"
mkdir -p "$OUT"

cut -f1,2 "$FAI" > "$OUT/chrom.sizes"

features() {
    grep -v "^#" "$GFF" | awk -v types="$1" 'BEGIN{n=split(types,t,","); for(i=1;i<=n;i++) keep[t[i]]=1}
        ($3 in keep) {print $1"\t"$4-1"\t"$5"\t.\t.\t"$7}'
}

features exon | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > "$OUT/exon.bed"

features gene,lnc_RNA,pseudogene | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge \
    | bedtools subtract -a - -b "$OUT/exon.bed" | sort -k1,1 -k2,2n | bedtools merge > "$OUT/intron.bed"

features gene | sort -k1,1 -k2,2n \
    | bedtools flank -g "$OUT/chrom.sizes" -l 2000 -r 0 -s | sort -k1,1 -k2,2n \
    | bedtools subtract -a - -b "$OUT/exon.bed" | sort -k1,1 -k2,2n \
    | bedtools subtract -a - -b "$OUT/intron.bed" | sort -k1,1 -k2,2n \
    | bedtools merge > "$OUT/promoter.bed"

cat "$OUT/exon.bed" "$OUT/intron.bed" "$OUT/promoter.bed" | sort -k1,1 -k2,2n | bedtools merge \
    | bedtools complement -i - -g "$OUT/chrom.sizes" > "$OUT/intergenic.bed"

for R in exon intron promoter intergenic; do
    awk -v r="$R" '{print $1"\t"$2"\t"$3"\t"r}' "$OUT/$R.bed"
done | sort -k1,1 -k2,2n | bgzip -c > "$OUT/genomic_overlaps.bed.gz"

tabix -p bed "$OUT/genomic_overlaps.bed.gz"
