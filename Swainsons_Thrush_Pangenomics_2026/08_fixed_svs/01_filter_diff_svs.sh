#!/bin/bash
# Strictly fixed (1 vs 0) and nearly fixed (>= 0.75 vs <= 0.25) SVs between inland and coastal birds.
# Requires AN >= 6 in both panels; undirected complex SVs are excluded.
# Coastal_ref is haploid but encoded as diploid, so its two alleles are removed from the coastal counts.

set -euo pipefail

TAB="pangenome/14vcf_aa_tab/variation_final.tab"
OUT="fixed_svs"
mkdir -p "$OUT"

cut -f1,2 refs/inland_ref_NO_UNPLACED.fasta.fai | sed '1i chr\tlength' > "$OUT/chr_sizes.txt"

awk -F'\t' -v OFS='\t' '
NR == 1 {
    print "chrom", "bedStart", "bedEnd", "type", "subtype", "inland_freq", "coastal_freq",
          "hybrid_freq", "hybrid_AN", "fixed_class", "direction"
    next
}
$5 ~ /^SV/ && $5 != "SV_Complex" {
    if ($16 < 6 || $19 < 8) next
    coastal_AN  = $19 - 2
    coastal_DAC = $18 - 2 * (($8 == $6 && $29 != "0|0") || ($8 == $7 && $29 == "0|0"))

    inland  = $15 / $16
    coastal = coastal_DAC / coastal_AN
    hybrid  = ($22 > 0) ? $21 / $22 : "NA"

    if      (coastal == 1 && inland == 0)        { dir = "coastal_fixed"; cls = "strictly_fixed" }
    else if (inland == 1 && coastal == 0)        { dir = "inland_fixed";  cls = "strictly_fixed" }
    else if (coastal >= 0.75 && inland <= 0.25)  { dir = "coastal_fixed"; cls = "nearly_fixed" }
    else if (inland >= 0.75 && coastal <= 0.25)  { dir = "inland_fixed";  cls = "nearly_fixed" }
    else next

    print $1, $2, $3, $4, $5, inland, coastal, hybrid, $22, cls, dir
}' "$TAB" > "$OUT/diff_svs.tsv"
