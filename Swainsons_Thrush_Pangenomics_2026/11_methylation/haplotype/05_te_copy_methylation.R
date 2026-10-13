# Mean per-copy methylation of full-length and fragmented TE copies per individual, compared across taxon

library(data.table)

DIR <- "methylation/te_copies"
COLS <- c("sample_hap", "status", "chrom", "start", "end", "divergence", "class_family", "sum_cov", "sum_mod", "n_cpg")
TAXON <- c(AH14K04 = "hybrid", MF12K02 = "hybrid",
           AH19K01 = "coastal", AH19K05 = "coastal", JF08K03 = "coastal", JF08K05 = "coastal",
           AH19K17 = "inland", AH21K03 = "inland", KF14K01 = "inland", KF15K02 = "inland")

dt <- rbindlist(lapply(list.files(DIR, "\\.(fulllength|fragmented)\\.tsv$", full.names = TRUE),
                       fread, header = FALSE, col.names = COLS))
dt[, individual := sub("_hap[12]$", "", sample_hap)]
dt[, taxon := TAXON[individual]]
dt[, pct_meth := 100 * sum_mod / sum_cov]

ind <- dt[, .(mean_pct = mean(pct_meth), n_copies = .N), by = .(individual, taxon, status)]
fwrite(dcast(ind, individual + taxon ~ status, value.var = c("mean_pct", "n_copies")),
       file.path(DIR, "te_copy_meth_per_individual.tsv"), sep = "\t")

stats <- rbindlist(lapply(c("fulllength", "fragmented"), function(st) {
    sub <- ind[status == st]
    pw <- pairwise.wilcox.test(sub$mean_pct, sub$taxon, p.adjust.method = "BH", exact = FALSE)
    pw <- setnames(as.data.table(as.table(pw$p.value)), c("group1", "group2", "p_adj"))[!is.na(p_adj)]
    pw[, `:=`(status = st, kw_p = kruskal.test(sub$mean_pct ~ sub$taxon)$p.value)]
}))
fwrite(stats, file.path(DIR, "te_copy_meth_stats.tsv"), sep = "\t")
