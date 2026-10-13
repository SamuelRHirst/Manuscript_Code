# Within each comparison, do hypomethylated DMRs overlap TE sequence more often than hypermethylated DMRs?
# One-sided Fisher exact test on direction x TE overlap, BH-corrected across comparisons within each annotation.
# Both groups come from the same caller, thresholds and callable sequence, so no genomic null is needed.
library(data.table)

DIR <- "methylation/dmr/te_enrichment"

rows <- list()

# whole-TE sources
a <- fread(file.path(DIR, "dmr_te_enrichment_callable.tsv"))
a <- a[null_model == "callable" & subset != "All" & te_source != "Union"]
for (cmp in unique(a$comparison)) for (src in unique(a$te_source)) {
  hypo  <- a[comparison == cmp & te_source == src & subset == "HypoMethylated"]
  hyper <- a[comparison == cmp & te_source == src & subset == "HyperMethylated"]
  if (!nrow(hypo) || !nrow(hyper)) next
  rows[[length(rows) + 1]] <- data.table(comparison = cmp, annotation = src,
    hypo_hit = hypo$n_overlap, hypo_n = hypo$n_regions,
    hyper_hit = hyper$n_overlap, hyper_n = hyper$n_regions)
}

# TE classes and families
b <- fread(file.path(DIR, "dmr_te_class_enrichment.tsv"))
b <- b[null_model == "callable"]
for (cmp in unique(b$comparison)) for (ann in unique(b$annotation)) {
  hypo  <- b[comparison == cmp & annotation == ann & subset == "HypoMethylated"]
  hyper <- b[comparison == cmp & annotation == ann & subset == "HyperMethylated"]
  if (!nrow(hypo) || !nrow(hyper)) next
  rows[[length(rows) + 1]] <- data.table(comparison = cmp, annotation = ann,
    hypo_hit = hypo$n_overlap, hypo_n = hypo$n_regions,
    hyper_hit = hyper$n_overlap, hyper_n = hyper$n_regions)
}

dt <- rbindlist(rows)
dt[, `:=`(
  hypo_pct  = 100 * hypo_hit / hypo_n,
  hyper_pct = 100 * hyper_hit / hyper_n)]

ft <- function(hh, hn, yh, yn) {
  m <- matrix(c(hh, hn - hh, yh, yn - yh), nrow = 2)
  f <- fisher.test(m, alternative = "greater")
  list(or = unname(f$estimate), p = f$p.value,
       lo = f$conf.int[1], hi = f$conf.int[2])
}
res <- dt[, {
  r <- ft(hypo_hit, hypo_n, hyper_hit, hyper_n)
  .(odds_ratio = r$or, ci_lo = r$lo, p_value = r$p)
}, by = .(comparison, annotation, hypo_hit, hypo_n, hyper_hit, hyper_n, hypo_pct, hyper_pct)]
res[, fdr := p.adjust(p_value, method = "BH"), by = annotation]
setorder(res, comparison, p_value)
fwrite(res, file.path(DIR, "dmr_direction_fisher.tsv"), sep = "\t")

