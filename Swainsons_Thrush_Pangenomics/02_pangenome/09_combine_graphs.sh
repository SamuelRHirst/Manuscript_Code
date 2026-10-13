#!/bin/bash
# Whole-genome GFA from the 41 chromosome graphs (vg v1.63.1)

set -euo pipefail

cd pangenome
mkdir -p 06b_combined_pggb

vg combine 06pggb/*.smooth.final.gfa > 06b_combined_pggb/combined.vg
vg convert -f 06b_combined_pggb/combined.vg > 06b_combined_pggb/swth_pangenome.gfa
vg stats -lz 06b_combined_pggb/combined.vg

rm 06b_combined_pggb/combined.vg
