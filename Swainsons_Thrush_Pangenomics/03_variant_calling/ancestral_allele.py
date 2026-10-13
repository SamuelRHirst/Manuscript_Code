#!/usr/bin/env python3
# Ancestral allele from a single outgroup genotype.
# Homozygous or half-missing calls give the allele; heterozygous or fully missing give '.'
# stdin/stdout: CHROM POS REF ALT GT  ->  CHROM POS REF ALT AA

import re
import sys


def ancestral(ref, alt, gt):
    a, b = re.split(r"[/|]", gt)
    if a == "." and b == ".":
        return "."
    if a != "." and b != "." and a != b:
        return "."
    idx = int(a if a != "." else b)
    return ([ref] + alt.split(","))[idx]


for line in sys.stdin:
    chrom, pos, ref, alt, gt = line.rstrip("\n").split("\t")
    print(chrom, pos, ref, alt, ancestral(ref, alt, gt), sep="\t")
