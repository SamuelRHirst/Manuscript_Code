#!/usr/bin/env python3
# Chromosome offset of a tile's tree sequence, found by matching its site positions to the input SNPs
# Usage: tile_offset.py <ts.trees> <snp_positions.txt> <tile>

import sys, numpy as np, tskit
ts = tskit.load(sys.argv[1]); label = sys.argv[3]
V = np.unique(np.loadtxt(sys.argv[2], dtype=np.int64, ndmin=1))
P = np.array([int(round(s.position)) for s in ts.sites()], dtype=np.int64)
if len(P) == 0 or len(V) == 0:
    print("%s\tNA\t0\t%d" % (label, len(P))); sys.exit()
probe = P[:20]
cands = V - probe[0]
M = probe[None, :] + cands[:, None]
idx = np.clip(np.searchsorted(V, M), 0, len(V) - 1)
score = (V[idx] == M).sum(axis=1)
best = int(cands[int(np.argmax(score))])
idx = np.clip(np.searchsorted(V, P + best), 0, len(V) - 1)
frac = float(np.mean(V[idx] == P + best))
print("%s\t%d\t%.5f\t%d" % (label, best, frac, len(P)))
