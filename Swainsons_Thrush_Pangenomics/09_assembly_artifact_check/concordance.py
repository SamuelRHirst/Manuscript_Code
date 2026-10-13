#!/usr/bin/env python3
# Do read-based Sniffles2 calls confirm the pangenome SVs? A pangenome SV is confirmed by a call of
# compatible type (DEL-DEL, INS-INS or INS-DUP) within 1 kb with size ratio >= 0.5. Each bird is then
# compared for presence, and differentiated SVs are checked for the same subspecies split.
# Test sets: differentiated SVs, and other eligible SVs (all in the blocks, 20,000 sampled in the background).
# Usage: concordance.py <pangenome_svs.tab> <merged_calls.tsv> <outdir>

import sys, bisect, random, csv
from collections import defaultdict

pg_path, calls_path, outdir = sys.argv[1:4]

INLAND = ["AH19K17", "AH21K03", "KF14K01", "KF15K02"]
COASTAL = ["AH19K01", "AH19K05", "JF08K03", "JF08K05"]
HYBRID = ["AH14K04", "MF12K02"]
BIRDS = INLAND + COASTAL + HYBRID
POS_TOL, SIZE_RATIO, N_BG = 1000, 0.5, 20000


def region_of(chrom, pos):
    mb = pos / 1e6
    if chrom == "Chromosome_01" and 46.0 <= mb < 59.5: return "chr1_core"
    if chrom == "Chromosome_01" and 25.0 <= mb < 100:  return "chr1_flank"
    if chrom == "Chromosome_05" and 21 <= mb < 39:     return "chr5_core"
    if chrom == "Chromosome_05" and 21 <= mb < 60:     return "chr5_flank"
    if chrom in ("Chromosome_W",) or "_Z_" in chrom:   return None
    return "background"


def alleles(gt):
    a = [x for x in gt.replace("/", "|").split("|") if x not in (".", "")]
    return a or None


pg = []
with open(pg_path) as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ix = {h: i for i, h in enumerate(hdr)}
    for b in BIRDS + ["Coastal_ref"]:
        assert b in ix, f"missing genotype column {b}"
    for line in fh:
        f = line.rstrip("\n").split("\t")
        sub = f[ix["subtype"]]
        if not sub.startswith("SV"): continue
        chrom = f[0].replace("Inland_ref#1#", "")
        pos = int(f[1])
        reg = region_of(chrom, pos)
        if reg is None: continue
        iDAC, iAN = int(f[ix["inland_DAC"]]), int(f[ix["inland_AN"]])
        cDAC, cAN = int(f[ix["coastal_DAC"]]), int(f[ix["coastal_AN"]])
        elig = sub != "SV_Complex" and iAN >= 6 and cAN >= 8
        if not elig: continue
        ref, alt, aa, cref = f[ix["ref"]], f[ix["alt"]], f[ix["aa"]], f[ix["Coastal_ref"]]
        cAN2 = cAN - 2
        cDAC2 = cDAC - 2 * ((aa == ref and cref != "0|0") or (aa == alt and cref == "0|0"))
        fi, fc = iDAC / iAN, cDAC2 / cAN2
        cls = "none"
        if (fc == 1 and fi == 0) or (fi == 1 and fc == 0): cls = "strictly_fixed"
        elif (fc >= 0.75 and fi <= 0.25) or (fi >= 0.75 and fc <= 0.25): cls = "nearly_fixed"
        rl, amax, amin = int(f[ix["base_allele_len"]]), int(f[ix["alt_len_max"]]), int(f[ix["alt_len_min"]])
        ins, dele = amax - rl, rl - amin
        if max(ins, dele) < 50: svtype, size = "OTHER", max(ins, dele)
        elif ins >= dele: svtype, size = "INS", ins
        else: svtype, size = "DEL", dele
        gts = {b: alleles(f[ix[b]]) for b in BIRDS}
        pg.append(dict(chrom=chrom, pos=pos, region=reg, subtype=sub, fixed_class=cls,
                       svtype=svtype, size=size, gts=gts))

random.seed(1)
bg_nondiff = [r for r in pg if r["fixed_class"] == "none" and r["region"] == "background"]
bg_sample = set(map(id, random.sample(bg_nondiff, min(N_BG, len(bg_nondiff)))))
tests = []
for r in pg:
    if r["fixed_class"] != "none":
        r["set"] = "differentiated"; tests.append(r)
    elif r["region"] != "background" or id(r) in bg_sample:
        r["set"] = "eligible_nondiff"; tests.append(r)

calls = defaultdict(list)
with open(calls_path) as fh:
    chdr = fh.readline().rstrip("\n").split("\t")
    samples = chdr[5:]
    missing = [b for b in BIRDS if b not in samples]
    if missing:
        sys.exit(f"Sniffles VCF lacks samples {missing}")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        t = f[2]
        if t not in ("INS", "DEL", "DUP"): continue
        pos = int(f[1])
        try: svlen = abs(int(f[3]))
        except ValueError: svlen = 0
        if t == "DUP" or svlen == 0:
            try: svlen = int(f[4]) - pos
            except ValueError: pass
        gts = {s: alleles(g) for s, g in zip(samples, f[5:])}
        calls[f[0]].append((pos, t, svlen, gts))
for c in calls: calls[c].sort(key=lambda x: x[0])
cpos = {c: [x[0] for x in v] for c, v in calls.items()}


def nonref(a):
    return None if a is None else any(x != "0" for x in a)


def dosage(a):
    return None if a is None else sum(x != "0" for x in a), (0 if a is None else len(a))


def freq(gts, birds):
    nr = tot = 0
    for b in birds:
        a = gts.get(b)
        if a is None: continue
        nr += sum(x != "0" for x in a); tot += len(a)
    return (nr / tot if tot else None), tot


def best_match(r):
    if r["svtype"] == "OTHER" or r["chrom"] not in calls: return None
    ok_types = {"DEL"} if r["svtype"] == "DEL" else {"INS", "DUP"}
    lo = bisect.bisect_left(cpos[r["chrom"]], r["pos"] - POS_TOL)
    hi = bisect.bisect_right(cpos[r["chrom"]], r["pos"] + POS_TOL)
    best, key = None, None
    for c in calls[r["chrom"]][lo:hi]:
        if c[1] not in ok_types or c[2] <= 0: continue
        ratio = min(c[2], r["size"]) / max(c[2], r["size"])
        if ratio < SIZE_RATIO: continue
        k = (abs(c[0] - r["pos"]), -ratio)
        if key is None or k < key: best, key = c, k
    return best


rows = []
for r in tests:
    carriers = [b for b in BIRDS if b in samples and nonref(r["gts"][b])]
    m = best_match(r)
    agree = n_cmp = 0
    if m:
        for b in BIRDS:
            if b not in samples: continue
            p, s = nonref(r["gts"][b]), nonref(m[3].get(b))
            if p is None or s is None: continue
            n_cmp += 1; agree += (p == s)
    rep = ""
    if r["set"] == "differentiated":
        pi, _ = freq(r["gts"], INLAND); pc, _ = freq(r["gts"], COASTAL)
        rep = "no_call"
        if m:
            si, ti = freq(m[3], INLAND); sc, tc = freq(m[3], COASTAL)
            if si is None or sc is None or ti < 6 or tc < 6: rep = "insufficient_genotypes"
            else:
                hi_i = pi is not None and pc is not None and pi > pc
                (sh, sl) = (si, sc) if hi_i else (sc, si)
                rep = "replicated" if (sh >= 0.75 and sl <= 0.25) else \
                      ("same_direction" if sh - sl >= 0.5 else "not_replicated")
    rows.append(dict(set=r["set"], region=r["region"], chrom=r["chrom"], pos=r["pos"], subtype=r["subtype"],
                     fixed_class=r["fixed_class"], svtype=r["svtype"], size=r["size"], n_carriers=len(carriers),
                     matched=int(m is not None), match_pos=m[0] if m else "", match_type=m[1] if m else "",
                     match_len=m[2] if m else "", birds_compared=n_cmp, birds_agree=agree, replication=rep))

with open(f"{outdir}/concordance_per_sv.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)

summ = defaultdict(lambda: defaultdict(int))
for x in rows:
    if x["svtype"] == "OTHER" or x["n_carriers"] == 0: key_skip = True
    else: key_skip = False
    k = (x["set"], x["region"])
    s = summ[k]
    s["n"] += 1
    if key_skip: s["not_testable"] += 1; continue
    s["testable"] += 1; s["matched"] += x["matched"]
    s["birds_compared"] += x["birds_compared"]; s["birds_agree"] += x["birds_agree"]
    if x["replication"]: s[x["replication"]] += 1

order = ["chr1_core", "chr1_flank", "chr5_core", "chr5_flank", "background"]
with open(f"{outdir}/concordance_summary.tsv", "w") as fh:
    cols = ["set", "region", "n", "testable", "pct_matched", "pct_bird_agreement",
            "pct_replicated", "pct_same_direction", "pct_not_replicated", "pct_insufficient", "pct_no_call"]
    fh.write("\t".join(cols) + "\n")
    for k in sorted(summ, key=lambda k: (k[0], order.index(k[1]))):
        s = summ[k]; t = s["testable"] or 1
        pct = lambda v: f"{100 * v / t:.1f}"
        out = [k[0], k[1], s["n"], s["testable"], pct(s["matched"]),
               f"{100 * s['birds_agree'] / s['birds_compared']:.1f}" if s["birds_compared"] else "NA"]
        out += [pct(s[c]) if k[0] == "differentiated" else "" for c in
                ("replicated", "same_direction", "not_replicated", "insufficient_genotypes", "no_call")]
        fh.write("\t".join(map(str, out)) + "\n")
