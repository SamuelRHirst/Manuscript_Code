# DMR overlap with RepeatMasker TE classes and the 25 most abundant families (>= 2 Mb), against two nulls that
# keep DMR widths and place DMRs on the callable footprint of that comparison: position only, and additionally
# requiring >= 10 reference CpGs per placed interval (methbat's own minimum).
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })

SIGDIR <- "methylation/dmr/signature"
OUTDIR <- "methylation/dmr/te_enrichment"
FASTA  <- "refs/inland_ref_NO_UNPLACED.fasta"
FAI    <- "refs/inland_ref_NO_UNPLACED.fasta.fai"
RM_OUT <- "repeats/repeatmasker/inland_ref/inland_genome_folded_renamed.fasta.out"

N_PERM <- 1000; MIN_CPG <- 10L; MIN_FAMILY_BP <- 2e6; OVERSAMPLE <- 120L
set.seed(42)
FALLBACK_COUNT <- 0L

genome <- fread(FAI, header = FALSE, select = 1:2, col.names = c("chrom", "length"))
chrom_names <- genome$chrom

merge_intervals <- function(g) {
  if (!nrow(g)) return(data.table(chrom = character(), start = integer(), end = integer()))
  setorder(g, chrom, start, end)
  g[, grp := cumsum(start > data.table::shift(cummax(end), fill = -1L)), by = chrom]
  m <- g[, .(start = min(start), end = max(end)), by = .(chrom, grp)]
  m[, grp := NULL]; setorder(m, chrom, start); m[]
}

dna <- readDNAStringSet(FASTA); names(dna) <- sub("\\s.*$", "", names(dna))
cpg_pos <- list()
for (ch in chrom_names) if (ch %in% names(dna)) cpg_pos[[ch]] <- start(matchPattern("CG", dna[[ch]]))
rm(dna); invisible(gc())

rm_dt <- fread(cmd = sprintf("awk 'NR>3 && NF>=11 {print $5\"\\t\"$6-1\"\\t\"$7\"\\t\"$11}' %s", RM_OUT),
               header = FALSE, col.names = c("chrom", "start", "end", "classfam"), fill = TRUE)
rm_dt <- rm_dt[chrom %in% chrom_names & !is.na(start) & !is.na(end) & end > start]
rm_te <- rm_dt[!grepl("^(Simple_repeat|Low_complexity|Satellite|rRNA|tRNA|snRNA|scRNA|srpRNA|ARTEFACT|Unknown)", classfam)]
rm_te[, class := sub("/.*$", "", classfam)]

te_sets <- list()
for (cl in unique(rm_te$class)) {
  b <- merge_intervals(rm_te[class == cl, .(chrom, start, end)])
  if (b[, sum(as.numeric(end - start))] >= MIN_FAMILY_BP) te_sets[[paste0("class:", cl)]] <- b
}
for (fm in rm_te[, .N, by = classfam][order(-N)][seq_len(min(25, .N)), classfam]) {
  b <- merge_intervals(rm_te[classfam == fm, .(chrom, start, end)])
  if (b[, sum(as.numeric(end - start))] >= MIN_FAMILY_BP) te_sets[[paste0("family:", fm)]] <- b
}

index_set <- function(bed) lapply(split(bed, bed$chrom), function(x) list(s = x$start, e = x$end))
set_idx <- lapply(te_sets, index_set)

overlap_flag <- function(regions, tidx) {
  hits <- logical(nrow(regions)); rc <- regions$chrom; rs <- regions$start; re <- regions$end
  for (ch in unique(rc)) {
    k <- which(rc == ch); t <- tidx[[ch]]
    if (is.null(t)) next
    j <- findInterval(re[k] - 1L, t$s); ok <- j >= 1L
    hits[k[ok]] <- t$e[j[ok]] > rs[k[ok]]
  }
  hits
}
count_cpg <- function(regions) {
  out <- integer(nrow(regions)); rc <- regions$chrom; rs <- regions$start; re <- regions$end
  for (ch in unique(rc)) {
    k <- which(rc == ch); p <- cpg_pos[[ch]]
    if (is.null(p)) next
    out[k] <- findInterval(re[k], p) - findInterval(rs[k], p)
  }
  out
}

build_ci <- function(callable) {
  w <- as.numeric(callable$end - callable$start)
  list(chrom = callable$chrom, start = callable$start, cum = cumsum(w), total = sum(w))
}
draw_raw <- function(widths, ci) {
  u <- runif(length(widths)) * ci$total
  idx <- pmin(findInterval(u, ci$cum) + 1L, length(ci$start))
  st <- ci$start[idx] + as.integer(u - c(0, ci$cum)[idx])
  data.table(chrom = ci$chrom[idx], start = st, end = st + widths)
}
draw_cpg <- function(widths, ci) {
  n <- length(widths)
  i_all <- rep(seq_len(n), each = OVERSAMPLE)
  w_all <- rep(widths,     each = OVERSAMPLE)
  cand  <- draw_raw(w_all, ci)
  cand[, i := i_all]
  cand <- cand[count_cpg(cand) >= MIN_CPG]

  out <- data.table(chrom = rep(NA_character_, n), start = rep(NA_integer_, n), end = rep(NA_integer_, n))
  if (nrow(cand)) {
    sel <- cand[, .SD[1L], by = i]
    out$chrom[sel$i] <- sel$chrom; out$start[sel$i] <- sel$start; out$end[sel$i] <- sel$end
  }
  miss <- which(is.na(out$chrom))
  if (length(miss)) {
    fb <- draw_raw(widths[miss], ci)
    out$chrom[miss] <- fb$chrom; out$start[miss] <- fb$start; out$end[miss] <- fb$end
    FALLBACK_COUNT <<- FALLBACK_COUNT + length(miss)
  }
  out[]
}

comparisons <- c("coastal_vs_inland", "coastal_vs_hybrid", "inland_vs_hybrid")
res <- list()
for (cmp in comparisons) {
  stats <- fread(cmd = paste("grep -v '^##'", file.path(SIGDIR, paste0(cmp, ".signature_stats.tsv"))),
                 sep = "\t", header = TRUE)
  setnames(stats, "summary_label", "label"); stats <- stats[chrom %in% chrom_names]
  ci  <- build_ci(merge_intervals(stats[, .(chrom, start, end)]))
  dmr <- stats[label %in% c("HyperMethylated", "HypoMethylated"), .(chrom, start, end, label)]

  for (subname in c("HyperMethylated", "HypoMethylated")) {
    reg <- dmr[label == subname]; if (!nrow(reg)) next
    widths <- reg$end - reg$start

    perms <- list(
      callable     = lapply(seq_len(N_PERM), function(i) draw_raw(widths, ci)),
      callable_cpg = lapply(seq_len(N_PERM), function(i) draw_cpg(widths, ci)))

    for (nm in names(te_sets)) {
      obs <- overlap_flag(reg, set_idx[[nm]]); of <- mean(obs)
      for (nullname in names(perms)) {
        nf <- vapply(perms[[nullname]], function(p) mean(overlap_flag(p, set_idx[[nm]])), numeric(1))
        res[[length(res) + 1]] <- data.table(comparison = cmp, subset = subname, annotation = nm,
          null_model = nullname, n_regions = nrow(reg), n_overlap = sum(obs),
          pct_overlap = 100 * of, pct_null_mean = 100 * mean(nf),
          fold_enrichment = of / mean(nf),
          perm_pvalue = (sum(nf >= of) + 1) / (N_PERM + 1))
      }
    }
  }
}

final <- rbindlist(res)
final[, fdr := p.adjust(perm_pvalue, method = "BH"), by = .(null_model, annotation, subset)]
fwrite(final, file.path(OUTDIR, "dmr_te_class_enrichment.tsv"), sep = "\t")

