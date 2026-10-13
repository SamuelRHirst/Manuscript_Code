#!/bin/bash
# Assembly quality per 100 kb window and per SV site on the autosomes:
#   windows.tsv   reference gaps, satellite, RepeatMasker and EDTA bp, haplotypes covering the window,
#                 contig ends, graph depth
#   sv_sites.tsv  each SV site with its differentiation class, missingness, contig ends within 10 and
#                 50 kb, and overlap with satellite or gaps (+-1 kb)
# Inputs: reference N-gaps and satellite-array BEDs for the inland reference.

set -euo pipefail

THREADS=8
FAI="refs/inland_ref_NO_UNPLACED.fasta.fai"
NGAPS="refs/inland_ngaps.bed"
SAT="refs/inland_satellite.bed"
EDTA="repeats/edta/inland_ref/inland_genome_folded_renamed.fasta.mod.EDTA.TEanno.gff3"
RMGFF="repeats/repeatmasker/inland_ref/inland_genome_folded_renamed.fasta.out.gff"
VARTAB="pangenome/14vcf_aa_tab/variation_final.tab"
OUT="artifact_check"
TMP="$OUT/tmp"
mkdir -p "$TMP"

HAPS=$(for s in AH14K04 AH19K01 AH19K05 AH19K17 AH21K03 JF08K03 JF08K05 KF14K01 KF15K02 MF12K02; do echo ${s}_hap1 ${s}_hap2; done)

cut -f1 "$FAI" | grep -vE '_Z_|_W$' > "$TMP/chroms.txt"
grep -wFf "$TMP/chroms.txt" "$FAI" | cut -f1,2 > "$TMP/genome.txt"
keep() { awk 'NR == FNR {k[$1]; next} $1 in k' "$TMP/chroms.txt" -; }

bedtools makewindows -g "$TMP/genome.txt" -w 100000 | sort -k1,1 -k2,2n > "$TMP/win.bed"

tail -n +2 "$NGAPS" | cut -f1-3 | keep | sort -k1,1 -k2,2n | bedtools merge > "$TMP/gaps.bed"
cut -f1-3 "$SAT" | keep | sort -k1,1 -k2,2n | bedtools merge > "$TMP/sat.bed"
grep -v '^#' "$EDTA"  | awk -v OFS='\t' '{print $1, $4-1, $5}' | keep | sort -k1,1 -k2,2n | bedtools merge > "$TMP/edta.bed"
grep -v '^#' "$RMGFF" | awk -v OFS='\t' '{print $1, $4-1, $5}' | keep | sort -k1,1 -k2,2n | bedtools merge > "$TMP/rm.bed"
for F in gaps sat edta rm; do
    bedtools coverage -a "$TMP/win.bed" -b "$TMP/$F.bed" -sorted -g "$TMP/genome.txt" | cut -f1-3,5 > "$TMP/win_$F.tsv"
done

# haplotype coverage and contig ends, from the wfmash mappings used to partition the pangenome
: > "$TMP/contig_ends.bed"
: > "$OUT/hapcov.tsv"
for H in $HAPS; do
    cat pangenome/02pafs/${H}.vs.inland.paf pangenome/03remap/${H}.unaligned.vs.inland.paf 2>/dev/null \
        | awk -v OFS='\t' '{sub(/^Inland_ref#1#/, "", $6); print}' > "$TMP/$H.paf"
    awk -v OFS='\t' '{print $6, $8, $9}' "$TMP/$H.paf" | keep | sort -k1,1 -k2,2n | bedtools merge > "$TMP/$H.cov.bed"
    bedtools coverage -a "$TMP/win.bed" -b "$TMP/$H.cov.bed" -sorted -g "$TMP/genome.txt" \
        | awk -v OFS='\t' -v h="$H" '{print $1, $2, $3, h, $7}' >> "$OUT/hapcov.tsv"
    # a mapping starting or ending within 1 kb of a contig end marks where that contig stops
    awk -v OFS='\t' -v h="$H" '{
        qlen = $2; qs = $3; qe = $4; st = $5; ts = $8; te = $9
        if (qs < 1000)        { p = (st == "+") ? ts : te; print $6, (p > 0 ? p-1 : 0), (p > 0 ? p : 1), h }
        if (qe > qlen - 1000) { p = (st == "+") ? te : ts; print $6, (p > 0 ? p-1 : 0), (p > 0 ? p : 1), h }
    }' "$TMP/$H.paf" | keep >> "$TMP/contig_ends.bed"
done
sort -k1,1 -k2,2n "$TMP/contig_ends.bed" > "$TMP/contig_ends.sorted.bed"

awk -v OFS='\t' '{k = $1 "\t" $2 "\t" $3; n[k]++; s[k] += $5; if ($5 >= 0.5) c[k]++}
    END {for (k in n) print k, c[k]+0, s[k]/n[k]}' "$OUT/hapcov.tsv" | sort -k1,1 -k2,2n > "$TMP/win_haps.tsv"
bedtools intersect -a "$TMP/win.bed" -b "$TMP/contig_ends.sorted.bed" -c -sorted -g "$TMP/genome.txt" > "$TMP/win_ends.tsv"

: > "$TMP/win_depth.tsv"
while read -r C; do
    OG=$(ls pangenome/06pggb/*"#${C}."*smooth.final.og)
    awk -v OFS='\t' -v c="$C" '$1 == c {print "Inland_ref#1#" $1, $2, $3}' "$TMP/win.bed" > "$TMP/win_$C.bed"
    odgi depth -i "$OG" -b "$TMP/win_$C.bed" -t "$THREADS" \
        | grep -v '^#' | awk -v OFS='\t' '{sub(/^Inland_ref#1#/, "", $1); print $1, $2, $3, $4}' >> "$TMP/win_depth.tsv"
done < "$TMP/chroms.txt"
sort -k1,1 -k2,2n "$TMP/win_depth.tsv" > "$TMP/win_depth.sorted.tsv"

paste "$TMP/win_gaps.tsv" <(cut -f4 "$TMP/win_sat.tsv") <(cut -f4 "$TMP/win_edta.tsv") \
      <(cut -f4 "$TMP/win_rm.tsv") <(cut -f4 "$TMP/win_ends.tsv") > "$TMP/win_part.tsv"
{
    printf "chrom\tstart\tend\tgap_bp\tsat_bp\tedta_bp\trm_bp\thaps_covered50\tmean_hap_cov\tcontig_ends\todgi_depth\n"
    awk -v OFS='\t' '
        FILENAME == ARGV[1] {h[$1 "\t" $2 "\t" $3] = $4 "\t" $5; next}
        FILENAME == ARGV[2] {d[$1 "\t" $2 "\t" $3] = $4; next}
        { k = $1 "\t" $2 "\t" $3
          print $1, $2, $3, $4, $5, $6, $7, (k in h) ? h[k] : "0\t0", $8, (k in d) ? d[k] : "NA" }
    ' "$TMP/win_haps.tsv" "$TMP/win_depth.sorted.tsv" "$TMP/win_part.tsv"
} > "$OUT/windows.tsv"

# SV sites, with the same eligibility and differentiation rules as 08_fixed_svs/01_filter_diff_svs.sh
awk -v OFS='\t' '
NR == 1 {next}
$5 ~ /^SV/ {
    c = $1; sub(/^Inland_ref#1#/, "", c)
    span = ($11 > 1) ? $11 : 1
    elig = ($5 != "SV_Complex" && $16 >= 6 && $19 >= 8)
    cls = "none"; dir = "none"
    if (elig) {
        cAN  = $19 - 2
        cDAC = $18 - 2 * (($8 == $6 && $29 != "0|0") || ($8 == $7 && $29 == "0|0"))
        fi = $15 / $16; fc = cDAC / cAN
        if      (fc == 1 && fi == 0)       { cls = "strictly_fixed"; dir = "coastal_fixed" }
        else if (fi == 1 && fc == 0)       { cls = "strictly_fixed"; dir = "inland_fixed" }
        else if (fc >= 0.75 && fi <= 0.25) { cls = "nearly_fixed";   dir = "coastal_fixed" }
        else if (fi >= 0.75 && fc <= 0.25) { cls = "nearly_fixed";   dir = "inland_fixed" }
    }
    print c, $2, $2 + span, $5, $11, $12, $13, elig, cls, dir, $17, $20, $23
}' "$VARTAB" | keep | sort -k1,1 -k2,2n > "$TMP/sv.bed"

bedtools window -a "$TMP/sv.bed" -b "$TMP/contig_ends.sorted.bed" -w 10000 -c \
    | bedtools window -a - -b "$TMP/contig_ends.sorted.bed" -w 50000 -c > "$TMP/sv_ends.bed"
bedtools slop -i "$TMP/sat.bed"  -g "$TMP/genome.txt" -b 1000 | bedtools merge > "$TMP/sat1k.bed"
bedtools slop -i "$TMP/gaps.bed" -g "$TMP/genome.txt" -b 1000 | bedtools merge > "$TMP/gaps1k.bed"
{
    printf "chrom\tstart\tend\tsubtype\tref_len\talt_len_max\talt_len_min\teligible\tfixed_class\tdirection\t"
    printf "inland_miss\tcoastal_miss\thybrid_miss\tcontig_ends_10kb\tcontig_ends_50kb\tsat_1kb\tgap_1kb\n"
    bedtools intersect -a "$TMP/sv_ends.bed" -b "$TMP/sat1k.bed" -c -sorted -g "$TMP/genome.txt" \
        | bedtools intersect -a - -b "$TMP/gaps1k.bed" -c -sorted -g "$TMP/genome.txt"
} > "$OUT/sv_sites.tsv"

rm -rf "$TMP"
