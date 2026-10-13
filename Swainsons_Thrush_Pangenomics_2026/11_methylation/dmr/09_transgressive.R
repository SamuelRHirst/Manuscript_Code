# Hybrid DMRs classified against both parents. DMRs from the two hybrid comparisons are merged into
# non-overlapping intervals, and coverage-weighted methylation is computed per individual over each interval
# from the 5mC pileups and averaged per group. Intervals where hybrids fall below both parental means are
# transgressively demethylated; their TE overlap is compared with all other hybrid DMRs (one-sided Fisher).
suppressMessages(library(data.table))

SIGDIR <- "methylation/dmr/signature"
PILEUP <- "methylation/dmr/pileup"
OUTDIR <- "methylation/dmr/transgressive"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
REPEATS <- "repeats"
GROUP <- c(AH19K01 = "coastal", AH19K05 = "coastal", JF08K03 = "coastal", JF08K05 = "coastal",
           AH19K17 = "inland",  AH21K03 = "inland",  KF14K01 = "inland",  KF15K02 = "inland",
           AH14K04 = "hybrid",  MF12K02 = "hybrid")

read_dmr <- function(cmp) {
  d <- fread(cmd = paste("grep -v '^##'", shQuote(file.path(SIGDIR, paste0(cmp, ".signature_regions.bed")))))
  setnames(d, 1:4, c("chrom", "start", "end", "label"))
  d[label %in% c("HypoMethylated", "HyperMethylated")][, cmp := cmp][]
}
dmr <- rbindlist(list(read_dmr("coastal_vs_hybrid"), read_dmr("inland_vs_hybrid")))
setorder(dmr, chrom, start, end)
# merge overlapping intervals across the two comparisons
dmr[, grp := cumsum(start > shift(cummax(end), fill = -1L)), by = chrom]
reg <- dmr[, .(start = min(start), end = max(end), labels = paste(sort(unique(label)), collapse = "/"),
               cmps = paste(sort(unique(cmp)), collapse = "+")), by = .(chrom, grp)]

# coverage-weighted mean methylation per individual over each interval
setkey(reg, chrom, start, end)
per_ind <- rbindlist(lapply(names(GROUP), function(s) {
  f <- file.path(PILEUP, paste0(s, ".5mC.bed.gz"))
  p <- fread(cmd = paste("zcat", shQuote(f), "| grep -v '^#'"), select = c(1, 2, 3, 8, 9, 10),
             col.names = c("chrom", "start", "end", "type", "cov", "mod"))
  p <- p[type == "Total"]
  ov <- foverlaps(p, reg, by.x = c("chrom", "start", "end"), nomatch = NULL)
  ov[, .(pct = 100 * sum(mod) / sum(cov), sites = .N), by = .(chrom, grp)][, `:=`(sample = s, group = GROUP[s])]
}))
wide <- dcast(per_ind[, .(pct = mean(pct)), by = .(chrom, grp, group)], chrom + grp ~ group, value.var = "pct")
m <- merge(reg, wide, by = c("chrom", "grp"))
m <- m[!is.na(hybrid) & !is.na(coastal) & !is.na(inland)]

m[, class := fifelse(hybrid < pmin(coastal, inland), "transgressive_low",
             fifelse(hybrid > pmax(coastal, inland), "transgressive_high", "intermediate"))]

# TE overlap: transgressive-low against all other hybrid DMRs
load_gff <- function(f) {
  g <- fread(cmd = paste("awk '!/^#/ {print $1\"\\t\"($4-1)\"\\t\"$5}'", shQuote(f)),
             col.names = c("chrom", "start", "end"))
  setorder(g, chrom, start, end)
  g[, grp2 := cumsum(start > shift(cummax(end), fill = -1L)), by = chrom]
  g[, .(start = min(start), end = max(end)), by = .(chrom, grp2)][, grp2 := NULL][]
}
te <- list(EDTA         = load_gff(file.path(REPEATS, "edta/inland_ref/inland_genome_folded_renamed.fasta.mod.EDTA.TEanno.gff3")),
           RepeatMasker = load_gff(file.path(REPEATS, "repeatmasker/inland_ref/inland_genome_folded_renamed.fasta.out.gff")))
res <- rbindlist(lapply(names(te), function(nm) {
  a <- copy(te[[nm]]); setkey(a, chrom, start, end)
  h <- foverlaps(m[, .(chrom, start, end, class)], a, by.x = c("chrom", "start", "end"), nomatch = NULL, which = TRUE)
  m[, te := .I %in% unique(h$xid)]
  tab <- table(factor(m$class == "transgressive_low", c(TRUE, FALSE)), factor(m$te, c(TRUE, FALSE)))
  ft <- fisher.test(tab, alternative = "greater")
  data.table(annotation = nm, trans_lo_with_te = tab[1, 1], trans_lo_n = sum(tab[1, ]),
             other_with_te = tab[2, 1], other_n = sum(tab[2, ]),
             odds = round(unname(ft$estimate), 2), p = signif(ft$p.value, 3))
}))
fwrite(m, file.path(OUTDIR, "hybrid_dmr_classification.tsv"), sep = "\t")
fwrite(res, file.path(OUTDIR, "transgressive_te_test.tsv"), sep = "\t")
