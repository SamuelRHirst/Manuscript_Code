# Age of the TE copies under DMRs. For each DMR, the youngest overlapping RepeatMasker TE copy gives its divergence
# from consensus and age class (young < 5%, intermediate 5-15%, old > 15%), and overlap with EDTA intact elements is
# recorded. Hypo- and hypermethylated DMRs are compared within each comparison (Wilcoxon on youngest divergence,
# Fisher on carrying a young copy), along with the share of TE-overlapping DMRs that are hypomethylated.
suppressPackageStartupMessages(library(data.table))
SIGDIR <- "methylation/dmr/signature"
REPEATS <- "repeats"
OUT  <- "methylation/dmr/te_age"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

rm_dt <- fread(cmd = sprintf("awk 'NR>3 && NF>=11 {print $5\"\\t\"$6-1\"\\t\"$7\"\\t\"$2\"\\t\"$11}' %s",
  file.path(REPEATS, "repeatmasker/inland_ref/inland_genome_folded_renamed.fasta.out")),
  header = FALSE, col.names = c("chrom", "start", "end", "div", "classfam"), fill = TRUE)
rm_dt <- rm_dt[!is.na(start) & !is.na(end) & end > start & !is.na(div)]
rm_te <- rm_dt[!grepl("^(Simple_repeat|Low_complexity|Satellite|rRNA|tRNA|snRNA|scRNA|srpRNA|ARTEFACT|Unknown)", classfam)]
rm_te[, class := sub("/.*$", "", classfam)]
rm_te[, age := cut(div, c(-Inf, 5, 15, Inf), labels = c("young", "mid", "old"), right = FALSE)]
stopifnot(nrow(rm_te) > 100000)

ed <- fread(cmd = sprintf("awk -F'\\t' 'NF>=9 && $1 !~ /^#/ && $9 !~ /Parent=/ {print $1\"\\t\"$4-1\"\\t\"$5}' %s",
  file.path(REPEATS, "edta/inland_ref/inland_genome_folded_renamed.fasta.mod.EDTA.intact.gff3")),
  header = FALSE, col.names = c("chrom", "start", "end"))
edta_te <- fread(cmd = sprintf("awk -F'\\t' 'NF>=9 && $1 !~ /^#/ {print $1\"\\t\"$4-1\"\\t\"$5}' %s",
  file.path(REPEATS, "edta/inland_ref/inland_genome_folded_renamed.fasta.mod.EDTA.TEanno.gff3")),
  header = FALSE, col.names = c("chrom", "start", "end"))
stopifnot(nrow(ed) > 1000, nrow(edta_te) > 100000)

ovl <- function(a, b) {
  a <- copy(a)[, `:=`(s = start + 1L, e = end)]; b <- copy(b)[, `:=`(s = start + 1L, e = end)]
  setkey(b, chrom, s, e)
  foverlaps(a, b, by.x = c("chrom", "s", "e"), by.y = c("chrom", "s", "e"), nomatch = 0L)
}

res <- list(); per <- list(); frame <- list()
for (cmp in c("coastal_vs_inland", "coastal_vs_hybrid", "inland_vs_hybrid")) {
  st <- fread(cmd = paste("grep -v '^##'", file.path(SIGDIR, paste0(cmp, ".signature_stats.tsv"))), sep = "\t", header = TRUE)
  d <- st[summary_label %in% c("HypoMethylated", "HyperMethylated"), .(chrom, start, end, label = summary_label)]
  d[, id := .I]
  o <- ovl(d, rm_te[, .(chrom, start, end, div, classfam, class, age)])
  y <- o[, .SD[which.min(div)], by = id][, .(id, youngest_div = div, youngest_class = classfam, youngest_age = age)]
  d <- merge(d, y, by = "id", all.x = TRUE)
  d[, rm_te := !is.na(youngest_div)]
  d[, edta_te := id %in% ovl(d, edta_te)$id]
  d[, edta_intact := id %in% ovl(d, ed)$id]
  d[, comparison := cmp]
  per[[cmp]] <- d
  for (lab in c("HypoMethylated", "HyperMethylated")) {
    x <- d[label == lab]
    res[[length(res) + 1]] <- data.table(comparison = cmp, direction = lab, dmrs = nrow(x), with_rm_te = sum(x$rm_te),
      median_youngest_div = median(x$youngest_div, na.rm = TRUE),
      young_n = sum(x$youngest_age == "young", na.rm = TRUE), mid_n = sum(x$youngest_age == "mid", na.rm = TRUE),
      old_n = sum(x$youngest_age == "old", na.rm = TRUE), edta_intact = sum(x$edta_intact))
  }
  a <- d[label == "HypoMethylated" & rm_te]; b <- d[label == "HyperMethylated" & rm_te]
  w <- if (nrow(a) > 1 && nrow(b) > 1) wilcox.test(a$youngest_div, b$youngest_div, exact = FALSE)$p.value else NA
  f <- fisher.test(matrix(c(sum(a$youngest_age == "young"), sum(a$youngest_age != "young"),
                            sum(b$youngest_age == "young"), sum(b$youngest_age != "young")), 2))
  res[[length(res)]][, `:=`(p_wilcox_div_hypo_vs_hyper = w, p_fisher_young_hypo_vs_hyper = f$p.value)]
  for (src in c("rm_te", "edta_te")) {
    te <- d[get(src) == TRUE]; nte <- d[get(src) == FALSE]
    frame[[length(frame) + 1]] <- data.table(comparison = cmp, annotation = fifelse(src == "rm_te", "RepeatMasker", "EDTA"),
      te_dmrs = nrow(te), te_hypo = sum(te$label == "HypoMethylated"), te_pct_hypo = 100 * mean(te$label == "HypoMethylated"),
      other_dmrs = nrow(nte), other_hypo = sum(nte$label == "HypoMethylated"), other_pct_hypo = 100 * mean(nte$label == "HypoMethylated"))
  }
}
res <- rbindlist(res, fill = TRUE); per <- rbindlist(per); frame <- rbindlist(frame)
fwrite(res, file.path(OUT, "dmr_te_age_summary.tsv"), sep = "\t")
fwrite(per, file.path(OUT, "dmr_te_age_per_dmr.tsv"), sep = "\t")
fwrite(frame, file.path(OUT, "te_dmr_hypo_share.tsv"), sep = "\t")

