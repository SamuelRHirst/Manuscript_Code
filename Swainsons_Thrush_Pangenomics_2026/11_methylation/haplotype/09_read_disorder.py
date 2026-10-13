#!/usr/bin/env python3
# Read-level methylation disorder at TE copies for one pseudo-haplotype. Each primary read (MAPQ >= 20) in
# the haplotype's own alignment contributes its confident CpG calls inside a copy (methylated if p >= 0.8,
# unmethylated if p <= 0.2). Per copy, over reads with >= 4 calls in it: reads, discordant reads (both
# states present), intermediate reads (25-75% methylated) and the summed methylated fraction.
# Usage: 09_read_disorder.py <sample_hap> <out.tsv>

import sys, bisect, collections
import pysam

TE = "methylation/te_copies"
ALN = "methylation/haplotype"
MIN_CALLS, MIN_MAPQ = 4, 20
METH_HI, METH_LO = 204, 51

sh, out = sys.argv[1], sys.argv[2]

copies = collections.defaultdict(list)
meta = []
for status in ("fulllength", "fragmented"):
    with open(f"{TE}/{sh}.{status}.tsv") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            chrom, s, e = f[2], int(f[3]), int(f[4])
            copies[chrom].append((s, e, len(meta)))
            meta.append((chrom, s, e, f[1], f[5], f[6]))
for c in copies:
    copies[c].sort()
starts = {c: [x[0] for x in v] for c, v in copies.items()}
n_reads = [0] * len(meta); disc = [0] * len(meta); inter = [0] * len(meta); sfrac = [0.0] * len(meta)

bam = pysam.AlignmentFile(f"{ALN}/{sh}/{sh}.bam")
for chrom in copies:
    cs, cv = starts[chrom], copies[chrom]
    for r in bam.fetch(chrom):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < MIN_MAPQ:
            continue
        mb = r.modified_bases
        if not mb:
            continue
        refpos = r.get_reference_positions(full_length=True)
        per = collections.defaultdict(lambda: [0, 0])
        for key, calls in mb.items():
            if key[2] != "m":
                continue
            for qpos, q in calls:
                if q < 0 or METH_LO < q < METH_HI:
                    continue
                p = refpos[qpos]
                if p is None:
                    continue
                k = bisect.bisect_right(cs, p) - 1
                if k < 0 or not (cv[k][0] <= p + 1 <= cv[k][1]):
                    continue
                a = per[cv[k][2]]; a[0] += 1; a[1] += q >= METH_HI
        for i, (n, m) in per.items():
            if n < MIN_CALLS:
                continue
            fr = m / n
            n_reads[i] += 1; sfrac[i] += fr
            if 0 < m < n:
                disc[i] += 1
            if 0.25 <= fr <= 0.75:
                inter[i] += 1

with open(out, "w") as o:
    o.write("sample_hap\tchrom\tstart\tend\tstatus\tdivergence\tclass_family\tn_reads\tdiscordant\tintermediate\tsum_frac\n")
    for i, (chrom, s, e, st, dv, cf) in enumerate(meta):
        if n_reads[i]:
            o.write(f"{sh}\t{chrom}\t{s}\t{e}\t{st}\t{dv}\t{cf}\t{n_reads[i]}\t{disc[i]}\t{inter[i]}\t{sfrac[i]:.3f}\n")
