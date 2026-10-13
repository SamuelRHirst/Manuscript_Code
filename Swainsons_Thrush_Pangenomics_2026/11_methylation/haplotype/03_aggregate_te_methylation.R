# Per-individual CpG methylation (hap1 + hap2 pooled by summing coverage and modified reads) and tests:
#   Kruskal-Wallis across taxon for CpG count and each sequence class
#   exact paired Wilcoxon across individuals for TE vs non-TE and fragmented vs full-length TE

library(data.table)

DIR <- "methylation/te_classes"
TAXON <- c(AH14K04 = "hybrid", MF12K02 = "hybrid",
           AH19K01 = "coastal", AH19K05 = "coastal", JF08K03 = "coastal", JF08K05 = "coastal",
           AH19K17 = "inland", AH21K03 = "inland", KF14K01 = "inland", KF15K02 = "inland")

pool <- function(pattern, cols) {
    dt <- rbindlist(lapply(list.files(DIR, pattern, full.names = TRUE), fread, header = FALSE, col.names = cols))
    dt[, individual := sub("_hap[12]$", "", sample_hap)]
    dt[, taxon := TAXON[individual]]
    dt[, lapply(.SD, sum), by = .(individual, taxon),
       .SDcols = setdiff(names(dt), c("sample_hap", "individual", "taxon"))]
}
pct <- function(mod, cov) 100 * mod / cov
kw <- function(d, cols) rbindlist(lapply(cols, function(m)
    data.table(metric = m, kw_p = kruskal.test(d[[m]] ~ factor(d$taxon))$p.value)))

base <- pool("\\.baseline\\.tsv$", c("sample_hap", "total_n", "total_cov", "total_mod",
                                     "rmTE_n", "rmTE_cov", "rmTE_mod", "rmNonTE_n", "rmNonTE_cov", "rmNonTE_mod",
                                     "edtaFull_n", "edtaFull_cov", "edtaFull_mod",
                                     "edtaNonFull_n", "edtaNonFull_cov", "edtaNonFull_mod"))
base[, `:=`(total_pct = pct(total_mod, total_cov), rmTE_pct = pct(rmTE_mod, rmTE_cov),
            rmNonTE_pct = pct(rmNonTE_mod, rmNonTE_cov), edtaFull_pct = pct(edtaFull_mod, edtaFull_cov),
            edtaNonFull_pct = pct(edtaNonFull_mod, edtaNonFull_cov))]
fwrite(base, file.path(DIR, "baseline_per_individual.tsv"), sep = "\t")
fwrite(kw(base, c("total_n", "total_pct", "rmTE_pct", "rmNonTE_pct", "edtaFull_pct", "edtaNonFull_pct")),
       file.path(DIR, "baseline_stats.tsv"), sep = "\t")

part <- pool("\\.partition\\.tsv$", c("sample_hap", "total_n", "total_cov", "total_mod", "full_n", "full_cov", "full_mod",
                                      "frag_n", "frag_cov", "frag_mod", "nonte_n", "nonte_cov", "nonte_mod"))
stopifnot(all(part$total_n == part$full_n + part$frag_n + part$nonte_n))
part[, `:=`(genome_pct = pct(total_mod, total_cov), nonte_pct = pct(nonte_mod, nonte_cov),
            te_pct = pct(frag_mod + full_mod, frag_cov + full_cov),
            frag_pct = pct(frag_mod, frag_cov), full_pct = pct(full_mod, full_cov))]
fwrite(part, file.path(DIR, "partition_per_individual.tsv"), sep = "\t")

paired <- function(x, y) wilcox.test(x, y, paired = TRUE, exact = TRUE)$p.value
stats <- rbind(kw(part, c("genome_pct", "nonte_pct", "te_pct", "frag_pct", "full_pct")),
               data.table(metric = "te_vs_nonte_paired_wilcoxon", kw_p = paired(part$te_pct, part$nonte_pct)),
               data.table(metric = "frag_vs_full_paired_wilcoxon", kw_p = paired(part$frag_pct, part$full_pct)))
setnames(stats, "kw_p", "p")
fwrite(stats, file.path(DIR, "partition_stats.tsv"), sep = "\t")
