#!/usr/bin/env python3
# Collapse one caller's overlapping calls: walking along each chromosome, a call joins the current
# cluster if it reciprocally overlaps the cluster span by >= 50%. Each bird takes its maximum copy number.
# Usage: merge_within_caller.py <in.tsv> <out.tsv>

import sys


def recip_overlap(a0, a1, b0, b1):
    ov = max(0, min(a1, b1) - max(a0, b0))
    return min(ov / (a1 - a0), ov / (b1 - b0)) if ov else 0.0


def collapse(cluster, n_birds):
    pos = min(int(r[1]) for r in cluster)
    end = max(int(r[2]) for r in cluster)
    gts = []
    for i in range(5, 5 + n_birds):
        vals = [int(r[i]) for r in cluster if r[i] != "NA"]
        gts.append(max(vals) if vals else "NA")
    total = sum(v for v in gts if v != "NA")
    return [cluster[0][0], pos, end, end - pos, total] + gts


def main(inpath, outpath, min_recip=0.5):
    with open(inpath) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    rows.sort(key=lambda r: (r[0], int(r[1])))
    n_birds = len(header) - 5

    clusters, cur = [], []
    for r in rows:
        p, e = int(r[1]), int(r[2])
        if cur and r[0] == cur[0][0] and recip_overlap(cs, ce, p, e) >= min_recip:
            cur.append(r)
            cs, ce = min(cs, p), max(ce, e)
        else:
            if cur:
                clusters.append(cur)
            cur, cs, ce = [r], p, e
    if cur:
        clusters.append(cur)

    with open(outpath, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for c in clusters:
            fh.write("\t".join(map(str, collapse(c, n_birds))) + "\n")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
