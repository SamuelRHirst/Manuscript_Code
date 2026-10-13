#!/bin/bash
# Genotype tables per caller, within-caller merge, then cross-caller loci

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
OUT="inversions/tables"
mkdir -p "$OUT"

python3 "$SCRIPTS/inversion_table.py" syri "$OUT/syri.tsv" inversions/syri/SyRI_combined.vcf
python3 "$SCRIPTS/inversion_table.py" svim "$OUT/svim.tsv" inversions/svim/SVIM_combined.vcf
python3 "$SCRIPTS/inversion_table.py" pggb "$OUT/pggb.tsv" pangenome/08vcf_norm/*.norm.vcf.gz

for C in syri svim pggb; do
    python3 "$SCRIPTS/merge_within_caller.py" "$OUT/$C.tsv" "$OUT/${C}_merged.tsv"
done

python3 "$SCRIPTS/combine_inversions.py" \
    "$OUT/syri_merged.tsv" "$OUT/svim_merged.tsv" "$OUT/pggb_merged.tsv" \
    inversions/combined_inversions.tsv
