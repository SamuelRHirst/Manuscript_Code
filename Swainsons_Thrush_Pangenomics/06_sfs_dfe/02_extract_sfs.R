# Unfolded SFS per feature and variant class from biallelic, polarized sites with no missing calls
#   All:     inland + coastal + hybrid haplotypes (n = 20)
#   Pure:    inland + coastal (n = 16), input to the DFE
#   Inland, Coastal: n = 8 each

library(data.table)

N_INLAND  <- 8L
N_COASTAL <- 8L
N_HYBRID  <- 4L

vars <- fread("sfs/variants_annotated.tsv", sep = "\t", header = TRUE, fill = TRUE)
vars <- vars[allele_count == 2 & polarized == TRUE]

num_cols <- c("inland_DAC", "coastal_DAC", "hybrid_DAC", "inland_MISS", "coastal_MISS", "hybrid_MISS")
vars[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]

# Coastal_ref is haploid but encoded as diploid; remove its derived alleles from the coastal count
vars[, coastal_DAC := coastal_DAC - 2L * as.integer(
    (aa == ref & Coastal_ref == "1|1") | (aa == alt & Coastal_ref == "0|0"))]

vars[, DAC_pure := inland_DAC + coastal_DAC]
vars[, DAC_all  := inland_DAC + coastal_DAC + hybrid_DAC]
vars[, type_class := fcase(
    subtype == "SNP", "SNP",
    subtype %in% c("DEL", "INS", "DEL_Complex", "INS_Complex", "INDEL_Complex"), "INDEL",
    subtype %in% c("SVDEL", "SVINS", "SV_Complex", "SVDEL_Complex", "SVINS_Complex"), "SV",
    default = NA_character_)]
vars <- vars[!is.na(type_class)]
vars[, feature := gsub("'", "", feature)]

features <- c("Intergenic", "CDS", "5UTR", "3UTR", "Intron")
types    <- c("SNP", "INDEL", "SV")

make_sfs <- function(data, dac_col, n_hap, group_name) {
    sub <- data[get(dac_col) > 0 & get(dac_col) < n_hap]
    rows <- list()
    for (feat in features) {
        for (typ in types) {
            d <- sub[feature == feat & type_class == typ]
            if (nrow(d) < 50) next
            counts <- d[, .N, by = .(dac = get(dac_col))][order(dac)]
            rows[[length(rows) + 1]] <- data.table(group = group_name, feature = feat, variant_type = typ,
                                                   dac = counts$dac, count = counts$N)
        }
    }
    rbindlist(rows)
}

sfs <- rbindlist(list(
    make_sfs(vars[inland_MISS == 0 & coastal_MISS == 0 & hybrid_MISS == 0 & !is.na(DAC_all)],
             "DAC_all", N_INLAND + N_COASTAL + N_HYBRID, "All"),
    make_sfs(vars[inland_MISS == 0 & coastal_MISS == 0 & !is.na(DAC_pure)],
             "DAC_pure", N_INLAND + N_COASTAL, "Pure"),
    make_sfs(vars[inland_MISS == 0 & !is.na(inland_DAC)], "inland_DAC", N_INLAND, "Inland"),
    make_sfs(vars[coastal_MISS == 0 & !is.na(coastal_DAC)], "coastal_DAC", N_COASTAL, "Coastal")
))

fwrite(sfs, "sfs/sfs_counts.csv")
