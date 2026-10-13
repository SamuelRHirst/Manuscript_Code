#!/bin/bash
# Repeat content by class for each genome (RepeatMasker) and for TE-derived SVs (GraffiTE), W excluded.
# RepeatMasker hits are swept in position order so each masked base is counted once.

set -euo pipefail

SCRIPTS=$(dirname "$(readlink -f "$0")")
OUT="repeats/class_tables"
mkdir -p "$OUT"

printf "genome\tclass\tbp\n" > "$OUT/repeatmasker_class_bp.tsv"

while IFS=$'\t' read -r NAME GENOME; do
    awk 'NR > 3 && NF >= 15 {print $5 "\t" $6 "\t" $7 "\t" $11}' "repeats/repeatmasker/$NAME/$(basename "$GENOME").out" \
        | sort -t $'\t' -k1,1 -k2,2n \
        | awk -F'\t' -v g="$NAME" -f "$SCRIPTS/class_map.awk" -e '
            $1 != chr { chr = $1; covered = 0 }
            $3 > covered {
                start = ($2 > covered + 1) ? $2 : covered + 1
                bp[cls($4)] += $3 - start + 1
                covered = $3
            }
            END { for (c in bp) printf "%s\t%s\t%d\n", g, c, bp[c] }' >> "$OUT/repeatmasker_class_bp.tsv"
done < "$SCRIPTS/genomes.tsv"

# SVs matching more than one class are counted as "Multiple classes"
zcat repeats/graffite/graffite_te_svs.vcf.gz | awk -F'\t' -f "$SCRIPTS/class_map.awk" -e '
    /^#/ || $1 ~ /Chromosome_W$/ { next }
    {
        match($8, /matching_classes=[^;]+/)
        n = split(substr($8, RSTART + 17, RLENGTH - 17), labels, ",")
        delete seen; k = 0
        for (i = 1; i <= n; i++) { c = cls(labels[i]); if (!(c in seen)) { seen[c] = 1; k++; only = c } }
        c = (k == 1) ? only : "Multiple classes"
        len = (length($4) > length($5)) ? length($4) : length($5)
        nsv[c]++; bp[c] += len
    }
    END {
        print "class\tn_sv\tsv_bp"
        for (c in nsv) printf "%s\t%d\t%d\n", c, nsv[c], bp[c]
    }' > "$OUT/graffite_class_counts.tsv"
