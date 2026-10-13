#!/usr/bin/env python3
# Per-tile statistics from the 20 SINGER posterior samples, per (sample, 100 kb window):
#   mean pairwise coalescence time within coastal, within inland and between subspecies (tskit branch mode),
#   most recent cross-subspecies coalescence, and the fraction of sequence where each subspecies, and both,
#   are monophyletic.
# With -ploidy 2, VCF sample i is on nodes 2i and 2i+1, so the VCF order
# (AH19K01 AH19K05 AH19K17 AH21K03 JF08K03 JF08K05 KF14K01 KF15K02) puts coastal on 0-3, 8-11 and inland on 4-7, 12-15.
# Usage: tile_stats.py <tile> <class> <chrom> <offset> <tile_dir> <out.tsv>

import sys
import os
import numpy as np
import tskit

COASTAL = [0, 1, 2, 3, 8, 9, 10, 11]
INLAND = [4, 5, 6, 7, 12, 13, 14, 15]
WIN = 100_000
SAMPLES = range(80, 100)


def window_grid(L):
    edges = list(np.arange(0, L, WIN))
    edges.append(L)
    return np.array(edges, dtype=float)


def tree_topology_stats(ts, edges, coastal_set):
    nw = len(edges) - 1
    w_tbmin = np.zeros(nw)
    w_monoco = np.zeros(nw)
    w_monoin = np.zeros(nw)
    w_recip = np.zeros(nw)
    w_span = np.zeros(nw)
    n_co = len(COASTAL)
    n_in = len(INLAND)
    for tree in ts.trees(tracked_samples=COASTAL):
        left, right = tree.interval
        tbmin = np.inf
        mono_co = False
        mono_in = False
        for u in tree.nodes():
            if tree.is_sample(u):
                continue
            c = tree.num_tracked_samples(u)
            n = tree.num_samples(u)
            if 0 < c < n:
                t = tree.time(u)
                if t < tbmin:
                    tbmin = t
            if c == n_co and n == n_co:
                mono_co = True
            if c == 0 and n == n_in:
                mono_in = True
        if not np.isfinite(tbmin):
            raise RuntimeError("tree at %.0f has no mixed node; is it multi-rooted?" % left)
        i = np.searchsorted(edges, left, side="right") - 1
        while i < nw and edges[i] < right:
            span = min(right, edges[i + 1]) - max(left, edges[i])
            if span > 0:
                w_tbmin[i] += tbmin * span
                w_monoco[i] += mono_co * span
                w_monoin[i] += mono_in * span
                w_recip[i] += (mono_co and mono_in) * span
                w_span[i] += span
            i += 1
    return (w_tbmin / w_span, w_monoco / w_span, w_monoin / w_span, w_recip / w_span, w_span)


def main():
    label, klass, chrom, start, tiledir, outfile = sys.argv[1:7]
    start = int(start)
    rows = []
    for k in SAMPLES:
        f = os.path.join(tiledir, "ts_%d.trees" % k)
        ts = tskit.load(f)
        if ts.num_samples != 16:
            raise RuntimeError("%s: %d samples, expected 16" % (f, ts.num_samples))
        if sorted(ts.samples()) != list(range(16)):
            raise RuntimeError("%s: sample node ids are not 0..15" % f)
        L = ts.sequence_length
        edges = window_grid(L)
        tw_co = np.ravel(ts.diversity(sample_sets=[COASTAL], windows=edges, mode="branch")) / 2
        tw_in = np.ravel(ts.diversity(sample_sets=[INLAND], windows=edges, mode="branch")) / 2
        tb = np.ravel(ts.divergence(sample_sets=[COASTAL, INLAND], windows=edges, mode="branch")) / 2
        tbmin, mco, min_, rec, span = tree_topology_stats(ts, edges, set(COASTAL))
        pos = ts.tables.sites.position
        nsites, _ = np.histogram(pos, bins=edges)
        for w in range(len(edges) - 1):
            rows.append((label, klass, chrom, start + int(edges[w]), start + int(edges[w + 1]), k,
                         int(nsites[w]), tw_co[w], tw_in[w], tb[w], tbmin[w],
                         mco[w], min_[w], rec[w]))
    bad = [r for r in rows if r[9] + 1e-6 < r[10]]
    if bad:
        raise RuntimeError("%d windows with T_between < T_between_min, which is impossible" % len(bad))
    tmp = outfile + ".tmp"
    with open(tmp, "w") as out:
        out.write("tile\tclass\tchrom\twin_start\twin_end\tsample\tn_sites\t"
                  "T_within_co\tT_within_in\tT_between\tT_between_min\t"
                  "mono_co\tmono_in\tmono_recip\n")
        for r in rows:
            out.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.5f\t%.5f\t%.5f\n" % r)
    os.replace(tmp, outfile)


if __name__ == "__main__":
    main()
