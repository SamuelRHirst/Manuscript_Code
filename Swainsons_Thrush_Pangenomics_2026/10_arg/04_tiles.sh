#!/bin/bash
# 1 Mb tiles: contiguous across chr1 25-100 Mb and chr5 21-60 Mb, plus 500 random background tiles
# from the other autosomes. Each tile carries the length-weighted mean recombination rate it covers.

set -euo pipefail

TILE=1000000
NBG=500
SEED=42
MINSNP=1000
OUT="arg/windows/tiles.tsv"

awk -F'\t' -v T=$TILE 'NR > 1 { t = int($2 / T) * T; n[$1 SUBSEP t] += $4; p[$1 SUBSEP t] += $5 * ($3 - $2) }
    END { for (k in n) { split(k, a, SUBSEP); printf "%s\t%d\t%d\t%d\t%.6g\n", a[1], a[2], a[2] + T, n[k], p[k] / T } }' \
    arg/windows/win_stats.tsv | sort -k1,1 -k2,2n > arg/windows/tile_counts.tsv

printf "label\tclass\tchrom\tstart\tend\tn_snp\tpi\tmean_r\n" > "$OUT"

add_tile() {
    local c=$1 s=$2 e=$3 cls=$4
    local map="arg/maps/${c}_recomb.txt"
    [ -s "$map" ] || return 0
    local r row
    r=$(awk -v s="$s" -v e="$e" '$2 > s && $1 < e { a = ($1 > s ? $1 : s); b = ($2 < e ? $2 : e); w += $3 * (b - a); t += b - a }
        END { if (t > 0) printf "%.6g", w / t }' "$map")
    row=$(awk -F'\t' -v c="$c" -v s="$s" '$1 == c && $2 == s { print $4 "\t" $5; exit }' arg/windows/tile_counts.tsv)
    [ -n "$r" ] && [ -n "$row" ] && [ "$(cut -f1 <<< "$row")" -ge $MINSNP ] || return 0
    printf "%s_%d-%dMb\t%s\t%s\t%d\t%d\t%s\t%s\n" "$c" $((s / 1000000)) $((e / 1000000)) "$cls" "$c" "$s" "$e" "$row" "$r" >> "$OUT"
}

for ((s = 25000000; s < 100000000; s += TILE)); do add_tile Chromosome_01 $s $((s + TILE)) focal_chr1; done
for ((s = 21000000; s < 60000000; s += TILE)); do add_tile Chromosome_05 $s $((s + TILE)) focal_chr5; done

awk -F'\t' -v MIN=$MINSNP "BEGIN { srand($SEED) } \$1 != \"Chromosome_01\" && \$1 != \"Chromosome_05\" && \$4 >= MIN { print rand() \"\t\" \$0 }" \
    arg/windows/tile_counts.tsv | sort -k1,1g | head -$NBG | cut -f2- | sort -k1,1 -k2,2n \
    | while IFS=$'\t' read -r c s e n p; do add_tile "$c" "$s" "$e" background; done
