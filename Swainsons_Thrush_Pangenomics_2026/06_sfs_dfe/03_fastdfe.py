#!/usr/bin/env python3
# DFE per feature and variant class with fastDFE, from the pure-subspecies SFS (n = 16).
# Neutral reference: intergenic SNPs for SNPs, intergenic INDELs for INDELs and SVs.
# Ancestral misidentification (eps) is estimated; 100 bootstrap replicates.

from pathlib import Path

import numpy as np
import pandas as pd
from fastdfe import BaseInference, Spectrum

SFS = Path("sfs")
N = 16
GENOME_BP = 1_130_000_000
INTERVALS = np.array([-np.inf, -100, -10, -1, 0])
BINS = ["strong", "moderate", "weak", "neutral"]


def bed_bp(path):
    with open(path) as fh:
        return sum(int(l.split("\t")[2]) - int(l.split("\t")[1]) for l in fh if l.strip())


feature_bp = {
    "CDS": bed_bp(SFS / "beds/cds_merged.bed"),
    "5UTR": bed_bp(SFS / "beds/five_utr_merged.bed"),
    "3UTR": bed_bp(SFS / "beds/three_utr_merged.bed"),
    "Intron": bed_bp(SFS / "beds/intron_merged.bed"),
}
feature_bp["Intergenic"] = GENOME_BP - sum(feature_bp.values())

counts = pd.read_csv(SFS / "sfs_counts.csv")
counts = counts[counts["group"] == "Pure"]


def spectrum(feature, vtype):
    sub = counts[(counts["feature"] == feature) & (counts["variant_type"] == vtype)]
    vec = np.zeros(N + 1)
    for dac, n in zip(sub["dac"], sub["count"]):
        if 1 <= dac < N:
            vec[dac] = n
    vec[0] = max(0, feature_bp[feature] - vec[1:].sum())
    return vec


neutral = {
    "SNP": Spectrum(spectrum("Intergenic", "SNP")),
    "INDEL": Spectrum(spectrum("Intergenic", "INDEL")),
    "SV": Spectrum(spectrum("Intergenic", "INDEL")),
}

results, eps = [], []
for feature in ["CDS", "Intron", "5UTR", "3UTR"]:
    for vtype in ["SNP", "INDEL", "SV"]:
        sel = spectrum(feature, vtype)
        if sel[1:].sum() < 50:
            continue
        inf = BaseInference(
            sfs_neut=neutral[vtype],
            sfs_sel=Spectrum(sel),
            do_bootstrap=True,
            n_bootstraps=100,
            fixed_params={"all": {"h": 0.5, "p_b": 0, "S_b": 1}},
        )
        inf.run()
        props, ci = inf.get_discretized(intervals=INTERVALS, confidence_intervals=True)
        mle = dict(inf.params_mle)
        eps.append({"feature": feature, "variant_type": vtype,
                    "eps": mle.get("eps"), "S_d": mle.get("S_d"), "b": mle.get("b")})
        for i, b in enumerate(BINS):
            results.append({"feature": feature, "variant_type": vtype, "bin": b,
                            "proportion": props[i],
                            "ci_halfwidth_lower": ci[0][i], "ci_halfwidth_upper": ci[1][i]})

pd.DataFrame(results).to_csv(SFS / "dfe_results.csv", index=False)
pd.DataFrame(eps).to_csv(SFS / "dfe_eps_estimates.csv", index=False)
