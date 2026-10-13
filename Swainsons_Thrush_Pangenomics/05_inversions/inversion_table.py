#!/usr/bin/env python3
# Inversion genotype table (copies per bird: 0, 1, 2 or NA) from one caller's VCF(s).
#   svim: Jasmine-merged SVIM-asm VCF, one diploid column per bird in 05_jasmine.sh order
#   syri: Jasmine-merged SyRI VCF, one haploid column per haplotype, named after the bird
#   pggb: vg deconstruct VCFs (normalized), inversions flagged INV=YES
# Usage: inversion_table.py <svim|syri|pggb> <out.tsv> <vcf> [vcf ...]

import gzip
import re
import sys

GROUP = {
    "AH19K17": "inland", "AH21K03": "inland", "KF14K01": "inland", "KF15K02": "inland",
    "coastal_ref": "coastal", "AH19K01": "coastal", "AH19K05": "coastal",
    "JF08K03": "coastal", "JF08K05": "coastal",
    "AH14K04": "hybrid", "MF12K02": "hybrid",
}
OUT_ORDER = list(GROUP)
SVIM_ORDER = ["AH14K04", "AH19K01", "AH19K05", "AH19K17", "AH21K03", "coastal_ref",
              "JF08K03", "JF08K05", "KF14K01", "KF15K02", "MF12K02"]


def info_dict(field):
    return dict(x.split("=", 1) if "=" in x else (x, True) for x in field.split(";"))


def alt_count(gt):
    gt = gt.split(":")[0]
    alleles = [a for a in re.split(r"[/|]", gt) if a != "."]
    if not alleles:
        return None
    return sum(a != "0" for a in alleles)


def bird_name(sample):
    name = re.sub(r"^\d+_", "", sample)
    name = re.sub(r"_hap[12]$", "", name)
    return "coastal_ref" if name.lower() == "coastal_ref" else name


def read_svim(path):
    for f in open(path):
        if f.startswith("#"):
            continue
        f = f.rstrip("\n").split("\t")
        info = info_dict(f[7])
        if info.get("SVTYPE") != "INV":
            continue
        start = int(float(info.get("AVG_START", f[1])))
        end = int(float(info.get("AVG_END", f[1])))
        if end <= start:
            end = int(info.get("END", f[1]))
        gts = {b: alt_count(g) for b, g in zip(SVIM_ORDER, f[9:])}
        yield f[0], start, end, gts


def read_syri(path):
    birds = []
    for f in open(path):
        if f.startswith("#CHROM"):
            birds = [bird_name(s) for s in f.rstrip("\n").split("\t")[9:]]
            continue
        if f.startswith("#"):
            continue
        f = f.rstrip("\n").split("\t")
        info = info_dict(f[7])
        if info.get("SVTYPE") != "INV":
            continue
        pos = int(f[1])
        end = int(info.get("END", pos))
        haps = {}
        for b, g in zip(birds, f[9:]):
            haps.setdefault(b, []).append(alt_count(g))
        gts = {b: None if all(h is None for h in v) else sum(h or 0 for h in v) for b, v in haps.items()}
        yield f[0], pos, end, gts


def read_pggb(path):
    samples = []
    for f in gzip.open(path, "rt"):
        if f.startswith("#CHROM"):
            samples = f.rstrip("\n").split("\t")[9:]
            continue
        if f.startswith("#"):
            continue
        f = f.rstrip("\n").split("\t")
        info = info_dict(f[7])
        inv = info.get("INV")
        if inv is None or (isinstance(inv, str) and inv.upper() not in ("YES", "TRUE", "1")):
            continue
        pos = int(f[1])
        if "SVLEN" in info:
            svlen = abs(int(info["SVLEN"]))
        elif "END" in info:
            svlen = abs(int(info["END"]) - pos)
        else:
            svlen = max(len(f[3]), len(f[4]))
        if svlen < 50:
            continue
        gts = {bird_name(s): alt_count(g) for s, g in zip(samples, f[9:])
               if s not in ("bCatMin1", "Inland_ref")}
        yield f[0].split("#")[-1], pos, pos + svlen, gts


def main(caller, out, vcfs):
    reader = {"svim": read_svim, "syri": read_syri, "pggb": read_pggb}[caller]
    rows = []
    for vcf in vcfs:
        for chrom, pos, end, gts in reader(vcf):
            total = sum(v for v in gts.values() if v is not None)
            rows.append([chrom, pos, end, abs(end - pos), total]
                        + ["NA" if gts.get(b) is None else gts[b] for b in OUT_ORDER])
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out, "w") as fh:
        fh.write("\t".join(["chrom", "pos", "end", "svlen", "total_copies"]
                           + [f"{GROUP[b]}_{b}" for b in OUT_ORDER]) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
