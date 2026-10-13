#!/bin/bash
# Reference-mapping check: per individual, mean pileup coverage in each DMR set, in hybrid-hypomethylated
# DMRs overlapping TEs, and across RepeatMasker TE sequence (every 50th interval), relative to that
# individual's genome-wide mean coverage

set -euo pipefail

SIGDIR="methylation/dmr/signature"
PILEUP="methylation/dmr/pileup"
OUT="methylation/dmr/coverage_bias"
mkdir -p "$OUT/regions"

COASTAL="AH19K01 AH19K05 JF08K03 JF08K05"
INLAND="AH19K17 AH21K03 KF14K01 KF15K02"
HYBRID="AH14K04 MF12K02"
taxon_of() {
    case " $COASTAL " in *" $1 "*) echo coastal; return;; esac
    case " $INLAND " in *" $1 "*) echo inland; return;; esac
    echo hybrid
}

grep -v '^#' repeats/repeatmasker/inland_ref/inland_genome_folded_renamed.fasta.out.gff \
    | awk -F'\t' -v OFS='\t' 'NF >= 5 {print $1, $4 - 1, $5}' | sort -k1,1 -k2,2n | bedtools merge > "$OUT/te_RepeatMasker.bed"
awk -F'\t' -v OFS='\t' '{print $1, $2 - 1, $2 - 1 + $3}' methylation/dmr/te_enrichment/graffite_raw.tsv \
    | sort -k1,1 -k2,2n | bedtools merge > "$OUT/te_GraffiTE.bed"

for CMP in coastal_vs_hybrid inland_vs_hybrid coastal_vs_inland; do
    for LAB in HypoMethylated HyperMethylated; do
        grep -v '^##' "$SIGDIR/$CMP.signature_stats.tsv" | tail -n +2 \
            | awk -F'\t' -v L="$LAB" -v OFS='\t' '$4 == L {print $1, $2, $3}' \
            | sort -k1,1 -k2,2n > "$OUT/regions/$CMP.$LAB.bed"
    done
done
for CMP in coastal_vs_hybrid inland_vs_hybrid; do
    for SRC in RepeatMasker GraffiTE; do
        bedtools intersect -u -a "$OUT/regions/$CMP.HypoMethylated.bed" -b "$OUT/te_$SRC.bed" \
            > "$OUT/regions/$CMP.HypoMethylated.TE_$SRC.bed"
    done
done
awk 'NR % 50 == 0' "$OUT/te_RepeatMasker.bed" > "$OUT/regions/baseline_TE_RepeatMasker.bed"

RES="$OUT/coverage_by_region.tsv"
printf 'individual\ttaxon\tregion_set\tn_sites\tmean_cov\n' > "$RES"
for IND in $COASTAL $INLAND $HYBRID; do
    TAX=$(taxon_of "$IND")
    BG="$PILEUP/$IND.5mC.bed.gz"
    read -r NS MC < <(zcat "$BG" | awk -F'\t' '!/^#/ {n++; s += $9} END {printf "%d %.4f\n", n, s / n}')
    printf '%s\t%s\tGENOME_WIDE\t%s\t%s\n' "$IND" "$TAX" "$NS" "$MC" >> "$RES"
    for F in "$OUT"/regions/*.bed; do
        RS=$(basename "$F" .bed)
        if [ ! -s "$F" ]; then printf '%s\t%s\t%s\t0\tNA\n' "$IND" "$TAX" "$RS" >> "$RES"; continue; fi
        read -r N M < <(tabix -R "$F" "$BG" \
            | awk -F'\t' '!/^#/ {n++; s += $9} END {if (n > 0) printf "%d %.4f\n", n, s / n; else print "0 NA"}')
        printf '%s\t%s\t%s\t%s\t%s\n' "$IND" "$TAX" "$RS" "$N" "$M" >> "$RES"
    done
done

# taxon means of coverage relative to each individual's genome-wide mean
awk -F'\t' -v OFS='\t' 'NR > 1 {
    if ($3 == "GENOME_WIDE") { gw[$1] = $5; next }
    cov[$3, $1] = $5; tax[$1] = $2; sets[$3]; inds[$1]
} END {
    print "region_set", "coastal", "inland", "hybrid"
    for (s in sets) {
        delete sum; delete cnt
        for (i in inds) if (cov[s, i] != "NA" && cov[s, i] != "") { sum[tax[i]] += cov[s, i] / gw[i]; cnt[tax[i]]++ }
        printf "%s", s
        n = split("coastal inland hybrid", T, " ")
        for (q = 1; q <= n; q++) printf "\t%s", (cnt[T[q]] > 0) ? sprintf("%.3f", sum[T[q]] / cnt[T[q]]) : "NA"
        printf "\n"
    }
}' "$RES" > "$OUT/relative_coverage.tsv"
