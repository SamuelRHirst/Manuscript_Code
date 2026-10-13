# Permutation enrichment of DMR overlap with TE sequence (EDTA, EDTA intact, RepeatMasker, their union, GraffiTE).
# DMR widths are kept and each DMR is placed at random 1,000 times, either uniformly along the chromosomes or
# uniformly along the sequence methbat signature actually tested in that comparison (the callable footprint).
library(data.table)

SIGDIR <- "methylation/dmr/signature"
OUTDIR <- "methylation/dmr/te_enrichment"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

EDTA_GFF <- "repeats/edta/inland_ref/inland_genome_folded_renamed.fasta.mod.EDTA.TEanno.gff3"
EDTA_INTACT_GFF <- "repeats/edta/inland_ref/inland_genome_folded_renamed.fasta.mod.EDTA.intact.gff3"
RM_GFF <- "repeats/repeatmasker/inland_ref/inland_genome_folded_renamed.fasta.out.gff"
GRAFFITE_RAW <- file.path(OUTDIR, "graffite_raw.tsv")
FAI <- "refs/inland_ref_NO_UNPLACED.fasta.fai"
N_PERM <- 1000
set.seed(42)

genome <- fread(FAI, header = FALSE, select = 1:2, col.names = c("chrom", "length"))
chrom_names <- genome$chrom

merge_intervals <- function(g) {
  setorder(g, chrom, start, end)
  g[, grp := cumsum(start > shift(cummax(end), fill = -1L)), by = chrom]
  merged <- g[, .(start = min(start), end = max(end)), by = .(chrom, grp)]
  merged[, grp := NULL]
  setorder(merged, chrom, start)
  merged[]
}

load_merged_gff <- function(gff_file) {
  g <- fread(cmd = paste("grep -v '^#'", gff_file), header = FALSE, sep = "\t", fill = TRUE)
  g <- g[V1 %in% chrom_names, .(chrom = V1, start = as.integer(V4) - 1L, end = as.integer(V5))]
  merge_intervals(g)
}

load_merged_graffite <- function(raw_tsv) {
  g <- fread(raw_tsv, header = FALSE, col.names = c("chrom", "pos", "reflen"))
  g <- g[chrom %in% chrom_names, .(chrom, start = pos - 1L, end = pos - 1L + reflen)]
  merge_intervals(g)
}

te_sources <- list(
  EDTA          = load_merged_gff(EDTA_GFF),
  EDTA_intact   = load_merged_gff(EDTA_INTACT_GFF),
  RepeatModeler = load_merged_gff(RM_GFF),
  GraffiTE      = load_merged_graffite(GRAFFITE_RAW)
)
te_sources$Union <- merge_intervals(rbindlist(list(copy(te_sources$EDTA),
                                                   copy(te_sources$RepeatModeler))))

index_te <- function(bed) {
  split_bed <- split(bed, bed$chrom)
  lapply(split_bed, function(x) list(s = x$start, e = x$end))
}
te_idx <- lapply(te_sources, index_te)

overlap_flag <- function(regions, tidx) {
  hits <- logical(nrow(regions))
  rc <- regions$chrom; rs <- regions$start; re <- regions$end
  for (ch in unique(rc)) {
    k <- which(rc == ch)
    t <- tidx[[ch]]
    if (is.null(t)) next
    j <- findInterval(re[k] - 1L, t$s)
    ok <- j >= 1L
    hits[k[ok]] <- t$e[j[ok]] > rs[k[ok]]
  }
  hits
}

# null 1: uniform along the chromosomes
chrom_len <- setNames(as.numeric(genome$length), genome$chrom)
draw_uniform <- function(widths) {
  n_r <- length(widths)
  pc <- sample(genome$chrom, n_r, replace = TRUE, prob = genome$length)
  ms <- chrom_len[pc] - widths
  ms[ms < 0] <- 0
  st <- as.integer(runif(n_r) * ms)
  data.table(chrom = pc, start = st, end = st + widths)
}

# null 2: uniform along the callable footprint
build_callable_index <- function(callable) {
  w <- as.numeric(callable$end - callable$start)
  list(chrom = callable$chrom, start = callable$start,
       cum = cumsum(w), total = sum(w))
}
draw_callable <- function(widths, ci) {
  u <- runif(length(widths)) * ci$total
  idx <- findInterval(u, ci$cum) + 1L
  idx <- pmin(idx, length(ci$start))
  offset <- as.integer(u - c(0, ci$cum)[idx])
  st <- ci$start[idx] + offset
  data.table(chrom = ci$chrom[idx], start = st, end = st + widths)
}

run_test <- function(regions, tidx, ci, label, cmp, te_name) {
  if (nrow(regions) == 0) return(NULL)
  widths <- regions$end - regions$start
  obs_hits <- overlap_flag(regions, tidx)
  obs_frac <- mean(obs_hits)
  nulls <- list(
    uniform  = replicate(N_PERM, mean(overlap_flag(draw_uniform(widths), tidx))),
    callable = replicate(N_PERM, mean(overlap_flag(draw_callable(widths, ci), tidx)))
  )
  rbindlist(lapply(names(nulls), function(nn) {
    nf <- nulls[[nn]]
    data.table(comparison = cmp, te_source = te_name, subset = label, null_model = nn,
               n_regions = nrow(regions), n_overlap = sum(obs_hits),
               mean_dmr_width = round(mean(widths)),
               pct_overlap = 100 * obs_frac, pct_null_mean = 100 * mean(nf),
               fold_enrichment = obs_frac / mean(nf),
               perm_pvalue = (sum(nf >= obs_frac) + 1) / (N_PERM + 1))
  }))
}

comparisons <- c("coastal_vs_inland", "coastal_vs_hybrid", "inland_vs_hybrid")
all_results <- list()

for (cmp in comparisons) {
  stats <- fread(cmd = paste("grep -v '^##'",
                             file.path(SIGDIR, paste0(cmp, ".signature_stats.tsv"))),
                 sep = "\t", header = TRUE)
  setnames(stats, "summary_label", "label")
  stats <- stats[chrom %in% chrom_names]

  # callable footprint = every segment methbat tested in this comparison
  callable <- merge_intervals(stats[, .(chrom, start, end)])
  ci <- build_callable_index(callable)

  dmr <- stats[label %in% c("HyperMethylated", "HypoMethylated"), .(chrom, start, end, label)]

  subsets <- list(All = dmr,
                  HyperMethylated = dmr[label == "HyperMethylated"],
                  HypoMethylated  = dmr[label == "HypoMethylated"])

  for (subname in names(subsets)) {
    reg <- subsets[[subname]]
    if (nrow(reg) == 0) next
    for (te_name in names(te_sources)) {
      res <- run_test(reg, te_idx[[te_name]], ci, subname, cmp, te_name)
      if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
    }
  }
}

final <- rbindlist(all_results)
final[, fdr := p.adjust(perm_pvalue, method = "BH"), by = .(null_model, te_source, subset)]
fwrite(final, file.path(OUTDIR, "dmr_te_enrichment_callable.tsv"), sep = "\t")
