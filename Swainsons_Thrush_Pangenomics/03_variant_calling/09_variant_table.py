#!/usr/bin/env python3
# Variant counts by allele count, genomic region and subtype, overall and per population
# (a variant counts toward a population if it carries the derived allele there). W excluded.
# Usage: 09_variant_table.py <variation_final.tab> <out.csv>

import csv
import sys
from collections import defaultdict

REGIONS = ["Genome_wide", "exon", "intron", "promoter", "intergenic"]
SUBTYPES = ["SNP", "MNP", "DEL", "INS", "SVDEL", "SVINS", "SV_Complex",
            "DEL_Complex", "INS_Complex", "INDEL_Complex", "SVDEL_Complex", "SVINS_Complex"]

tab, out = sys.argv[1:3]
counts = defaultdict(lambda: [0, 0, 0, 0])

with open(tab) as fh:
    next(fh)
    for line in fh:
        e = line.split("\t")
        if e[0].endswith("Chromosome_W"):
            continue
        group = "Bi-allelic" if int(e[13]) == 2 else "Multi-allelic"
        dac = [int(e[14]), int(e[17]), int(e[20])]
        for region in (e[8], "Genome_wide"):
            c = counts[(group, region, e[4])]
            c[0] += 1
            for i, d in enumerate(dac, 1):
                if d > 0:
                    c[i] += 1

with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["allele_count", "genomic_region", "subtype", "All", "inland", "coastal", "hybrid"])
    for group in ("Bi-allelic", "Multi-allelic"):
        for region in REGIONS:
            for st in SUBTYPES:
                if (group, region, st) in counts:
                    w.writerow([group, region, st] + counts[(group, region, st)])
