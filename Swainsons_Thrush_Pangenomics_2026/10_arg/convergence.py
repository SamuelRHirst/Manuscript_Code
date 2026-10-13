#!/usr/bin/env python3
# SINGER convergence per tile: number of trees at each MCMC iteration, change between iterations 70-79
# and 90-99, and rank correlation with iteration across the kept samples (80-99)
# Usage: convergence.py <arg_dir> <out_dir>

import sys, os, re
import numpy as np
A = sys.argv[1]; outdir = sys.argv[2]
def ranks(x):
    o = np.argsort(x, kind="mergesort"); r = np.empty(len(x)); r[o] = np.arange(len(x))
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1: r[m] = r[m].mean()
    return r
def spearman(x, y):
    rx, ry = ranks(np.asarray(x, float)), ranks(np.asarray(y, float))
    return float(np.corrcoef(rx, ry)[0, 1])
trace = open(os.path.join(outdir, "convergence_trace.tsv"), "w")
trace.write("tile\tclass\titeration\tn_trees\n")
summ = open(os.path.join(outdir, "convergence_summary.tsv"), "w")
summ.write("tile\tclass\tn_iter\tmean_70_79\tmean_90_99\trel_change_pct\trho_80_99\tcv_80_99_pct\n")
with open(os.path.join(A, "windows", "tiles.tsv")) as f:
    f.readline()
    for ln in f:
        L, C = ln.split("\t")[:2]
        vals = [int(m.group(1)) for m in re.finditer(r"Number of trees: (\d+)", open(os.path.join(A, "tiles", L, "singer.log")).read())]
        for i, v in enumerate(vals):
            trace.write("%s\t%s\t%d\t%d\n" % (L, C, i, v))
        if len(vals) != 100:
            summ.write("%s\t%s\t%d\tNA\tNA\tNA\tNA\tNA\n" % (L, C, len(vals))); continue
        v = np.array(vals, float)
        a, b, k = v[70:80].mean(), v[90:100].mean(), v[80:100]
        summ.write("%s\t%s\t%d\t%.1f\t%.1f\t%.2f\t%.3f\t%.2f\n" % (L, C, len(vals), a, b, 100 * (b - a) / a,
                   spearman(np.arange(20), k), 100 * k.std() / k.mean()))
trace.close(); summ.close()
