#!/bin/bash
# TE annotation of SVs with GraffiTE (https://github.com/cgroza/GraffiTE), then transferred to the biallelic VCF

set -euo pipefail

VCF="pangenome/15vcf_biallelic/ingroup_biallelic.vcf.gz"
REF="refs/inland_ref.pansn.fa"
LIB="refs/Catharus_TE_library.fasta"
OUT="repeats/graffite"
mkdir -p "$OUT"

# SVs only; GraffiTE also rejects variant IDs over 50 characters
bcftools view -i 'STRLEN(ALT)>=50 || STRLEN(REF)>=50' "$VCF" \
    | bcftools annotate --set-id '%CHROM\_%POS' -Oz -o "$OUT/graffite_input.vcf.gz"
tabix -p vcf "$OUT/graffite_input.vcf.gz"

nextflow run cgroza/GraffiTE \
    --vcf "$OUT/graffite_input.vcf.gz" \
    --reference "$REF" \
    --TE_library "$LIB" \
    --genotype false \
    --tsd_batch_size 5000 \
    --repeatmasker_threads 4 \
    --repeatmasker_memory 10G \
    --repeatmasker_time 24h \
    --out "$OUT/graffite_out"

bgzip -c "$OUT/graffite_out/3_TSD_search/pangenome.vcf" > "$OUT/pangenome.vcf.gz"
tabix -p vcf "$OUT/pangenome.vcf.gz"

bcftools annotate \
    -a "$OUT/pangenome.vcf.gz" \
    -c INFO/GR,INFO/n_hits,INFO/repeat_ids,INFO/matching_classes,INFO/fragmts,INFO/total_match_span \
    "$VCF" \
    | bcftools view -i 'INFO/matching_classes!="."' -Oz -o "$OUT/graffite_te_svs.vcf.gz"
tabix -p vcf "$OUT/graffite_te_svs.vcf.gz"
