#!/bin/bash
# SNP count and nucleotide diversity per 500 kb window, used to record variant density per tile

set -euo pipefail

WIN=500000
OUT="arg/windows"
mkdir -p "$OUT"

printf "chrom\tstart\tend\tn_snp\tpi\n" > "$OUT/win_stats.tsv"
for V in arg/vcf/Chromosome_*.vcf.gz; do
    C=$(basename "$V" .vcf.gz)
    bcftools query -f '%CHROM\t%POS[\t%GT]\n' "$V" | awk -v W=$WIN -v C="$C" -v OFS='\t' '
        {
            an = 0; ac = 0
            for (i = 3; i <= NF; i++) {
                g = $i; gsub(/\|/, "/", g); split(g, a, "/")
                for (k = 1; k <= 2; k++) if (a[k] == "0" || a[k] == "1") { an++; if (a[k] == "1") ac++ }
            }
            if (an < 2) next
            w = int($2 / W) * W
            n[w]++
            if (ac > 0 && ac < an) pi[w] += 2 * ac * (an - ac) / (an * (an - 1))
        }
        END { for (w in n) print C, w, w + W, n[w], pi[w] / W }' >> "$OUT/win_stats.tsv"
done
