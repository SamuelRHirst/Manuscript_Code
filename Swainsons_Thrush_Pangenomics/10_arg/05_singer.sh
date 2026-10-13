#!/bin/bash
# SINGER v0.1.9 on one tile, keeping posterior samples 80-99 as tree sequences
# Usage: 05_singer.sh <tile number (1-based row of tiles.tsv)> [seed]
# A tile that hits a SINGER numerical assertion fails identically on rerun; a different seed avoids it.

set -euo pipefail

N=$1
SEED=${2:-}
SINGER_DIR="singer-0.1.9-beta-linux-x86_64"
NE=200000
MU=3.44e-9

read -r LABEL CLASS CHR START END NSNP PI R < <(awk -F'\t' -v n=$((N + 1)) 'NR == n' arg/windows/tiles.tsv)
D="arg/tiles/$LABEL"
rm -rf "$D"
mkdir -p "$D"

bcftools view "arg/vcf/$CHR.vcf.gz" -r "$CHR:$START-$END" -Ov -o "$D/$LABEL.vcf"
python "$SINGER_DIR/index_vcf.py" "$D/$LABEL" 100000
SSTART=$(head -1 "$D/$LABEL.index" | cut -f1)
SEND=$(tail -1 "$D/$LABEL.index" | cut -f1)

"$SINGER_DIR/singer" -Ne $NE -m $MU -r "$R" -input "$D/$LABEL" -output "$D/out" \
    -start "$SSTART" -end "$SEND" -n 100 -thin 20 -ploidy 2 ${SEED:+-seed $SEED} > "$D/singer.log" 2>&1

python "$SINGER_DIR/convert_to_tskit" -input "$D/out" -output "$D/ts" -start 80 -end 100

rm -f "$D"/out_*.txt "$D/$LABEL.vcf" "$D/$LABEL.index"
