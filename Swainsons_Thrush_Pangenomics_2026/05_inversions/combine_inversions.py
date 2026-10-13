#!/usr/bin/env python3
# Combine SyRI, SVIM-asm and PGGB inversions into loci by complete-linkage clustering:
# a call joins a cluster only if it reciprocally overlaps every member by >= 50%.
# Calls with no thrush carrier (e.g. outgroup-only PGGB alleles) are dropped.
# Usage: combine_inversions.py <syri.tsv> <svim.tsv> <pggb.tsv> <out.tsv>

import sys
from collections import defaultdict

RECIP = 0.5
BIRDS = ["AH19K17", "AH21K03", "KF14K01", "KF15K02",
         "coastal_ref", "AH19K01", "AH19K05", "JF08K03", "JF08K05",
         "AH14K04", "MF12K02"]
GROUP = {b: g for g, bs in [("inland", BIRDS[:4]), ("coastal", BIRDS[4:9]), ("hybrid", BIRDS[9:])] for b in bs}
CALLERS = ["syri", "svim", "pggb"]


def bird_from_header(col):
    col = col.split("_", 1)[1]
    return "coastal_ref" if col.lower() == "coastal_ref" else col


def load(path, caller):
    recs = []
    with open(path) as fh:
        birds = [bird_from_header(c) for c in fh.readline().rstrip("\n").split("\t")[5:]]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gts = {b: None if v == "NA" else int(v) for b, v in zip(birds, f[5:])}
            if not any(v for v in gts.values()):
                continue
            recs.append(dict(chrom=f[0].split("#")[-1], pos=int(f[1]), end=int(f[2]), gts=gts, src=caller))
    return recs


def compatible(a, b):
    ov = min(a["end"], b["end"]) - max(a["pos"], b["pos"])
    la, lb = a["end"] - a["pos"], b["end"] - b["pos"]
    return ov > 0 and la > 0 and lb > 0 and ov / la >= RECIP and ov / lb >= RECIP


def main(syri, svim, pggb, out):
    by_chrom = defaultdict(list)
    for path, caller in zip((syri, svim, pggb), CALLERS):
        for r in load(path, caller):
            by_chrom[r["chrom"]].append(r)

    loci = []
    for chrom in sorted(by_chrom):
        clusters = []
        for r in sorted(by_chrom[chrom], key=lambda r: (r["pos"], r["end"])):
            for c in clusters:
                if all(compatible(r, m) for m in c):
                    c.append(r)
                    break
            else:
                clusters.append([r])
        loci.extend(clusters)

    with open(out, "w") as fh:
        fh.write("\t".join(["chrom", "pos", "end", "svlen", "total_copies"]
                           + [f"{GROUP[b]}_{b}" for b in BIRDS] + CALLERS) + "\n")
        for c in loci:
            pos, end = min(r["pos"] for r in c), max(r["end"] for r in c)
            gts = []
            for b in BIRDS:
                vals = [r["gts"][b] for r in c if r["gts"].get(b) is not None]
                gts.append(max(vals) if vals else None)
            total = sum(v for v in gts if v is not None)
            row = [c[0]["chrom"], pos, end, end - pos, total]
            row += ["NA" if v is None else v for v in gts]
            row += [int(any(r["src"] == s for r in c)) for s in CALLERS]
            fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:5])
