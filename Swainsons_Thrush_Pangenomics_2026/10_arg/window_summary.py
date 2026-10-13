#!/usr/bin/env python3
# Per 100 kb window, mean over the 20 posterior samples of reciprocal monophyly, between-subspecies
# coalescence, within-subspecies coalescence (coastal, inland and their mean) and site density,
# with each window assigned to a core block, its flank, or the background.
# Usage: window_summary.py <arg_dir>

import glob
import os
import sys

A = sys.argv[1]
P = os.path.join(A, "post")
CORE = {"Chromosome_01": (46.0e6, 59.5e6), "Chromosome_05": (21.0e6, 39.0e6)}

tiles = {}
for line in open(os.path.join(A, "windows", "tiles.tsv")).readlines()[1:]:
    f = line.rstrip("\n").split("\t")
    tiles[f[0]] = dict(cls=f[1], chrom=f[2], r=float(f[7]))

agg = {}
for fn in glob.glob(os.path.join(P, "per_tile", "*.tsv")):
    lines = open(fn).read().splitlines()
    ix = {c: i for i, c in enumerate(lines[0].split("\t"))}
    for line in lines[1:]:
        w = line.split("\t")
        a = agg.setdefault((w[ix["tile"]], int(w[ix["win_start"]])),
                           dict(n=0, mono=0, tb=0, tw=0, co=0, inl=0, sites=0, end=int(w[ix["win_end"]])))
        co, inl = float(w[ix["T_within_co"]]), float(w[ix["T_within_in"]])
        a["n"] += 1
        a["mono"] += float(w[ix["mono_recip"]])
        a["tb"] += float(w[ix["T_between"]])
        a["tw"] += (co + inl) / 2
        a["co"] += co
        a["inl"] += inl
        a["sites"] += int(w[ix["n_sites"]])

with open(os.path.join(P, "window_summary.tsv"), "w") as out:
    out.write("tile\tclass\tregion\tchrom\twin_start\twin_end\tsites_per_kb\tr\tmono_recip\t"
              "T_between\tT_within\tratio\tT_within_co\tT_within_in\n")
    for (tile, s), a in sorted(agg.items()):
        t, n, mid = tiles[tile], a["n"], (s + a["end"]) / 2
        if t["cls"] == "background":
            region = "background"
        elif t["chrom"] in CORE and CORE[t["chrom"]][0] <= mid < CORE[t["chrom"]][1]:
            region = "core_" + t["chrom"][-2:]
        else:
            region = "flank_" + t["chrom"][-2:]
        out.write("%s\t%s\t%s\t%s\t%d\t%d\t%.3f\t%.3e\t%.4f\t%.1f\t%.1f\t%.3f\t%.1f\t%.1f\n" % (
            tile, t["cls"], region, t["chrom"], s, a["end"],
            a["sites"] / n / ((a["end"] - s) / 1000), t["r"], a["mono"] / n,
            a["tb"] / n, a["tw"] / n, a["tb"] / a["tw"], a["co"] / n, a["inl"] / n))
