#!/bin/bash
# Is core-block enrichment of differentiated SVs explained by SV density or SNP differentiation?
#   A. differentiated SVs as a share of all SVs, and of SVs eligible for the filter in 01_filter_diff_svs.sh
#   B. differentiated SVs per 500 kb window in the cores vs background windows with SNP Fst > 0.8

set -euo pipefail

TAB="pangenome/14vcf_aa_tab/variation_final.tab"
FST="fst/snps_fst_500kb.windowed.weir.fst"
D="fixed_svs"

REGION='function region(c, p) {
    if (c == "Chromosome_01" && p >= 46e6 && p < 59.5e6) return "core_chr1"
    if (c == "Chromosome_05" && p >= 21e6 && p < 39e6)   return "core_chr5"
    if (c == "Chromosome_01" && p >= 25e6 && p < 100e6)  return "zone"
    if (c == "Chromosome_05" && p >= 21e6 && p < 60e6)   return "zone"
    return "background"
}'

{
echo "== per-SV normalization =="
awk -F'\t' -v OFS='\t' "$REGION"'
    FNR == NR { if (NR > 1 && $1 !~ /Chromosome_(W|Z)/) auto[$1] = $2; next }
    FNR == 1 { next }
    { c = $1; sub(/^.*#/, "", c); if (!(c in auto)) next; p = $2 }
    FILENAME ~ /variation_final/ {
        if ($5 !~ /^SV/) next
        r = region(c, p); all[r]++
        if ($5 != "SV_Complex" && $16 >= 6 && $19 >= 8) elig[r]++
        next
    }
    { diff[region(c, p)]++ }
    END {
        len["core_chr1"] = 13.5; len["core_chr5"] = 18
        for (k in auto) bg += auto[k] / 1e6
        len["background"] = bg - 75 - 39
        print "region", "Mb", "all_SV", "all_SV_per_Mb", "eligible_SV", "diff_SV", "diff_per_Mb", "diff_per_1000_all", "diff_per_1000_elig"
        for (r in len)
            printf "%s\t%.1f\t%d\t%.0f\t%d\t%d\t%.2f\t%.3f\t%.3f\n", r, len[r], all[r], all[r]/len[r], elig[r], diff[r],
                diff[r]/len[r], 1000*diff[r]/all[r], 1000*diff[r]/elig[r]
        for (r in len) if (r != "background")
            printf "fold per SV %s: all=%.1f eligible=%.1f (per Mb fold %.2f)\n", r,
                (diff[r]/all[r]) / (diff["background"]/all["background"]),
                (diff[r]/elig[r]) / (diff["background"]/elig["background"]),
                (diff[r]/len[r]) / (diff["background"]/len["background"])
    }' "$D/chr_sizes.txt" "$TAB" "$D/diff_svs.tsv"

echo
echo "== differentiated SVs per 500 kb window =="
awk -F'\t' "$REGION"'
    FILENAME ~ /diff_svs/ { if (FNR == 1) next; c = $1; sub(/^.*#/, "", c); n[c ":" int($2 / 500000)]++; next }
    FNR == 1 { next }
    {
        c = $1; sub(/^.*#/, "", c)
        if (c !~ /^Chromosome_/ || c ~ /Chromosome_(W|Z)/) next
        s = $2 - 1; k = n[c ":" int(s / 500000)] + 0; r = region(c, s + 250000)
        if (r ~ /^core/) { cs += k; cn++ }
        else if (r == "background") { if ($5 > 0.8) { hs += k; hn++ }; bs += k; bn++ }
    }
    END {
        printf "core windows: n=%d mean=%.2f\n", cn, cs / cn
        printf "background windows with SNP Fst > 0.8: n=%d mean=%.2f\n", hn, hs / hn
        printf "all background windows: n=%d mean=%.2f\n", bn, bs / bn
    }' "$D/diff_svs.tsv" "$FST"
} > "$D/sv_normalization.txt"
