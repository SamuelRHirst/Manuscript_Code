#!/usr/bin/env python3
# CDS, 5' UTR, 3' UTR and intron BEDs from a GFF3, per mRNA
# Usage: parse_gff.py <in.gff> <outdir>

import sys, os
from collections import defaultdict

def merge(ivs):
    if not ivs: return []
    s = sorted(ivs)
    m = [list(s[0])]
    for a, b in s[1:]:
        if a <= m[-1][1]: m[-1][1] = max(m[-1][1], b)
        else: m.append([a, b])
    return m

def subtract(a_list, b_list):
    res = []
    for a0, a1 in a_list:
        cur = a0
        j = 0
        while j < len(b_list) and b_list[j][1] < cur: j += 1
        while j < len(b_list) and b_list[j][0] <= a1:
            if b_list[j][0] > cur: res.append((cur, b_list[j][0] - 1))
            cur = max(cur, b_list[j][1] + 1)
            j += 1
        if cur <= a1: res.append((cur, a1))
    return res

def parse_attrs(s):
    d = {}
    for f in s.split(";"):
        f = f.strip()
        if "=" in f:
            k, v = f.split("=", 1)
            d[k] = v
    return d

gff_file, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)

mrna_info      = {}
exons_by_mrna  = defaultdict(list)
cds_by_mrna    = defaultdict(list)
mrna_order     = []

with open(gff_file) as fh:
    for line in fh:
        if line.startswith("#"): continue
        # maxsplit=8 so embedded tabs in attributes do not break unpacking
        p = line.rstrip("\n").split("\t", 8)
        if len(p) < 9: continue
        chrom, _, feat, s, e, _, strand, _, attr_str = p
        s, e = int(s), int(e)
        attrs = parse_attrs(attr_str)
        if feat == "mRNA":
            mid = attrs.get("ID", "")
            mrna_info[mid] = (chrom, strand)
            mrna_order.append(mid)
        elif feat == "exon":
            for p2 in attrs.get("Parent", "").split(","):
                exons_by_mrna[p2].append((s, e))
        elif feat == "CDS":
            for p2 in attrs.get("Parent", "").split(","):
                cds_by_mrna[p2].append((s, e))

cds_fh  = open(os.path.join(outdir, "cds.bed"),       "w")
utr5_fh = open(os.path.join(outdir, "five_utr.bed"),  "w")
utr3_fh = open(os.path.join(outdir, "three_utr.bed"), "w")
intr_fh = open(os.path.join(outdir, "intron.bed"),    "w")

for mid in mrna_order:
    if mid not in mrna_info: continue
    chrom, strand = mrna_info[mid]
    exons = merge(exons_by_mrna.get(mid, []))
    cdss  = merge(cds_by_mrna.get(mid, []))
    if not exons: continue

    for s, e in cdss:
        cds_fh.write(f"{chrom}\t{s-1}\t{e}\n")

    if cdss:
        utrs    = subtract(exons, cdss)
        cds_min = cdss[0][0]
        cds_max = cdss[-1][1]
        for s, e in utrs:
            if strand == "+":
                fh2 = utr5_fh if e < cds_min else utr3_fh
            else:
                fh2 = utr5_fh if s > cds_max else utr3_fh
            fh2.write(f"{chrom}\t{s-1}\t{e}\n")

    if len(exons) > 1:
        body = [(exons[0][0], exons[-1][1])]
        for s, e in subtract(body, exons):
            intr_fh.write(f"{chrom}\t{s-1}\t{e}\n")

for fh2 in [cds_fh, utr5_fh, utr3_fh, intr_fh]:
    fh2.close()
print("Done")
