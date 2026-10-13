#!/usr/bin/env python3
# Variant classification and per-population derived allele counts.
# Adapted from compute_var_table.py in harvardinformatics/scrub-jay-genomics.
# Usage: classify_variants.py <pggb_variation.tab> <sample_populations.txt> <out.tab>

import re
import sys

POPS = ["inland", "coastal", "hybrid"]


def read_samples(path):
    samples, pop = [], {}
    with open(path) as fh:
        for line in fh:
            s, p = line.split()
            samples.append(s)
            pop[s] = p
    return samples, pop


def alt_length(ref, alts, aa, longest):
    base = aa if aa != "." else ref
    var = [x for x in [ref] + alts if x != base] if aa != "." else alts
    if not var:
        var = alts
    lengths = [len(x) for x in var]
    return max(lengths) if longest else min(lengths)


def pop_counts(genotypes, samples, alleles, aa, pop):
    counts = {p: [0, 0, 0] for p in POPS}  # derived, called, missing
    for sample, gt in zip(samples, genotypes):
        c = counts[pop[sample]]
        for hap in re.split(r"[/|]", gt):
            if hap == ".":
                c[2] += 1
                continue
            c[1] += 1
            if (aa != "." and alleles[int(hap)] != aa) or (aa == "." and int(hap) != 0):
                c[0] += 1
    return counts


def subtype(vtype, polarized, base_len, alt_min, alt_max, n_alleles):
    if vtype in ("SNP", "MNP"):
        return vtype
    if not polarized:
        return "SV_Complex" if max(alt_min, alt_max, base_len) >= 50 else "INDEL_Complex"

    alt = alt_max if n_alleles == 2 else alt_min
    sv = max(alt_min, alt_max, base_len) >= 50
    prefix = "SV" if sv else ""
    threshold = 50 if sv else 2

    if base_len == 1 and alt >= threshold:
        return prefix + "INS"
    if base_len >= threshold and alt == 1:
        return prefix + "DEL"
    if base_len > alt:
        return prefix + "DEL_Complex"
    if base_len < alt_max:
        return prefix + "INS_Complex"
    return "SV_Complex" if sv else "INDEL_Complex"


def main(tab, sample_file, out):
    samples, pop = read_samples(sample_file)
    header = ["chrom", "bedStart", "bedEnd", "type", "subtype", "ref", "alt", "aa", "gr",
              "polarized", "base_allele_len", "alt_len_max", "alt_len_min", "allele_count"]
    for p in POPS:
        header += [f"{p}_DAC", f"{p}_AN", f"{p}_MISS"]
    header += samples

    with open(tab) as fin, open(out, "w") as fout:
        fout.write("\t".join(header) + "\n")
        for line in fin:
            e = line.rstrip("\n").split("\t")
            chrom, start, end, vtype, ref, alt, aa, gr = e[:8]
            vtype = vtype.upper()
            genotypes = e[8:]

            alts = alt.split(",")
            alleles = [ref] + alts
            polarized = aa != "."
            alt_max = alt_length(ref, alts, aa, True)
            alt_min = alt_length(ref, alts, aa, False)
            base_len = len(aa) if polarized else len(ref)

            row = [chrom, start, end, vtype,
                   subtype(vtype, polarized, base_len, alt_min, alt_max, len(alleles)),
                   ref, alt, aa, gr, str(polarized), str(base_len),
                   str(alt_max), str(alt_min), str(len(alleles))]
            counts = pop_counts(genotypes, samples, alleles, aa, pop)
            for p in POPS:
                row += [str(x) for x in counts[p]]
            fout.write("\t".join(row + genotypes) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
