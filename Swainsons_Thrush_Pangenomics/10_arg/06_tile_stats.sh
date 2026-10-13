#!/bin/bash
# Per-tile statistics in 100 kb windows. Tree-sequence coordinates are first placed on the
# chromosome by matching their site positions to the input SNPs.

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
OUT="arg/post"
mkdir -p "$OUT/per_tile"

printf "tile\toffset\tmatch_frac\tn_sites\n" > "$OUT/tile_offsets.tsv"
tail -n +2 arg/windows/tiles.tsv | while IFS=$'\t' read -r L C CHR ST EN rest; do
    LO=$(( ST > 1000000 ? ST - 1000000 : 1 ))
    bcftools query -r "$CHR:$LO-$((EN + 1000000))" -f '%POS\n' "arg/vcf/$CHR.vcf.gz" > "$OUT/pos.txt"
    python3 "$SCRIPTS/tile_offset.py" "arg/tiles/$L/ts_80.trees" "$OUT/pos.txt" "$L" >> "$OUT/tile_offsets.tsv"
done
rm "$OUT/pos.txt"

tail -n +2 arg/windows/tiles.tsv | while IFS=$'\t' read -r L C CHR rest; do
    O=$(awk -F'\t' -v l="$L" '$1 == l {print $2}' "$OUT/tile_offsets.tsv")
    python3 "$SCRIPTS/tile_stats.py" "$L" "$C" "$CHR" "$O" "arg/tiles/$L" "$OUT/per_tile/$L.tsv"
done

python3 "$SCRIPTS/convergence.py" arg "$OUT"
python3 "$SCRIPTS/window_summary.py" arg
