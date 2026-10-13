# DMRs by gene feature: each DMR is assigned one category by priority promoter (TSS -2000/+50 bp) > 5'UTR > CDS >
# 3'UTR > other exon > intron > intergenic, and CpG-island overlap is recorded. UTRs are exon sequence outside the
# coding span of each transcript. Overlap is compared with 1,000 random placements of the same DMRs within the
# callable sequence, requiring >= 10 reference CpGs per placement, and hypo- vs hypermethylated DMRs are contrasted
# with Fisher exact tests. MODE=run processes the comparisons in COMPS; MODE=combine pools them.
suppressPackageStartupMessages({ library(data.table); library(Biostrings) })

SIGDIR <- "methylation/dmr/signature"
OUTDIR <- "methylation/dmr/gene_features"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
REFDIR <- "refs"
FASTA  <- file.path(REFDIR, "inland_ref_NO_UNPLACED.fasta")
FAI    <- file.path(REFDIR, "inland_ref_NO_UNPLACED.fasta.fai")
GFF    <- file.path(REFDIR, "inland_genome_folded_renamed.gff")

N_PERM <- as.integer(Sys.getenv("N_PERM", "1000"))
COMPS  <- strsplit(Sys.getenv("COMPS", "coastal_vs_inland,coastal_vs_hybrid,inland_vs_hybrid"), ",")[[1]]
MODE   <- Sys.getenv("MODE", "run")
PERDIR <- file.path(OUTDIR, "per_comparison"); dir.create(PERDIR, showWarnings = FALSE)
PROM_UP <- 2000L; PROM_DOWN <- 50L; MIN_CPG <- 10L
CGI_WIN <- 200L; CGI_GC <- 0.50; CGI_OE <- 0.60; CGI_MERGEGAP <- 100L
PRIORITY <- c("promoter", "UTR5", "CDS", "UTR3", "other_exon", "intron", "intergenic")
set.seed(42)
if (MODE == "run") {

genome <- fread(FAI, header = FALSE, select = 1:2, col.names = c("chrom", "length"))
chrom_names <- genome$chrom
GENOME_BP <- sum(as.numeric(genome$length))

merge_intervals <- function(g, gap = 0L) {
  if (nrow(g) == 0) return(data.table(chrom = character(), start = integer(), end = integer()))
  g <- copy(g)[, .(chrom, start, end)]
  setorder(g, chrom, start, end)
  g[, grp := cumsum(start > data.table::shift(cummax(end), fill = -1L) + gap), by = chrom]
  m <- g[, .(start = min(start), end = max(end)), by = .(chrom, grp)]
  m[, grp := NULL]; setorder(m, chrom, start); m[]
}
# a minus b, both merged; sweep per chromosome
subtract_intervals <- function(a, b) {
  a <- merge_intervals(a); b <- merge_intervals(b); out <- list()
  for (ch in unique(a$chrom)) {
    ai <- a[chrom == ch]; bi <- b[chrom == ch]
    if (!nrow(bi)) { out[[length(out) + 1]] <- ai; next }
    bs <- bi$start; be <- bi$end; j0 <- 1L; res <- vector("list", nrow(ai))
    for (k in seq_len(nrow(ai))) {
      cs <- ai$start[k]; ce <- ai$end[k]
      while (j0 <= length(be) && be[j0] <= cs) j0 <- j0 + 1L
      cur <- cs; pieces <- list(); j <- j0
      while (j <= length(bs) && bs[j] < ce) {
        if (bs[j] > cur) pieces[[length(pieces) + 1]] <- c(cur, bs[j])
        cur <- max(cur, be[j]); j <- j + 1L
      }
      if (cur < ce) pieces[[length(pieces) + 1]] <- c(cur, ce)
      if (length(pieces)) { p <- do.call(rbind, pieces); res[[k]] <- data.table(chrom = ch, start = p[, 1], end = p[, 2]) }
    }
    out[[length(out) + 1]] <- rbindlist(res)
  }
  merge_intervals(rbindlist(out))
}
bp <- function(x) x[, sum(as.numeric(end - start))]

# CpG positions + CpG islands
dna <- readDNAStringSet(FASTA); names(dna) <- sub("\\s.*$", "", names(dna))
cpg_pos <- list(); cgi_list <- list()
for (ch in chrom_names) {
  if (!ch %in% names(dna)) next
  s <- dna[[ch]]; L <- length(s)
  cpg_pos[[ch]] <- start(matchPattern("CG", s))
  n <- L %/% CGI_WIN; if (n < 1) next
  v <- Views(s, start = seq(1L, by = CGI_WIN, length.out = n), width = CGI_WIN)
  bf <- letterFrequency(v, c("C", "G")); df <- dinucleotideFrequency(v)[, "CG"]
  nC <- bf[, "C"]; nG <- bf[, "G"]
  oe <- ifelse(nC > 0 & nG > 0, df * CGI_WIN / (nC * nG), 0)
  keep <- which((nC + nG) / CGI_WIN >= CGI_GC & oe >= CGI_OE)
  if (length(keep)) cgi_list[[length(cgi_list) + 1]] <- data.table(chrom = ch,
      start = as.integer((keep - 1L) * CGI_WIN), end = as.integer(keep * CGI_WIN))
}
rm(dna); invisible(gc())
cgi <- merge_intervals(rbindlist(cgi_list), gap = CGI_MERGEGAP)[end - start >= 200L]

# gene features
awk_cmd <- sprintf(paste0("awk -F'\\t' 'NF>=9 && $1 !~ /^[!#]/ {a=$9; for(i=10;i<=NF;i++) a=a\" \"$i; ",
  "id=\"\"; par=\"\"; if (match(a,/ID=[^;]+/)) id=substr(a,RSTART+3,RLENGTH-3); ",
  "if (match(a,/Parent=[^;]+/)) par=substr(a,RSTART+7,RLENGTH-7); ",
  "print $1\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$7\"\\t\"id\"\\t\"par}' %s"), GFF)
g <- fread(cmd = awk_cmd, sep = "\t", header = FALSE, quote = "", fill = TRUE,
           col.names = c("chrom", "feat", "start", "end", "strand", "id", "parent"))
g <- g[chrom %in% chrom_names]
g[, `:=`(start = as.integer(start) - 1L, end = as.integer(end))]
g <- g[!is.na(start) & !is.na(end) & end > start]
stopifnot(g[feat == "gene", .N] > 10000, g[feat == "exon", .N] > 100000, g[feat == "CDS", .N] > 100000)

genes <- merge_intervals(g[feat %in% c("gene", "pseudogene")])
exons <- merge_intervals(g[feat == "exon"])
cds   <- merge_intervals(g[feat == "CDS"])
tx    <- g[feat %in% c("mRNA", "transcript", "lnc_RNA")]
tss   <- tx[, .(chrom, pos = fifelse(strand == "-", end, start), strand)]
promoters <- merge_intervals(tss[, .(chrom,
  start = pmax(0L, fifelse(strand == "-", pos - PROM_DOWN, pos - PROM_UP)),
  end   = fifelse(strand == "-", pos + PROM_UP, pos + PROM_DOWN))])

# UTRs: per coding transcript, exon sequence left/right of the CDS span; orientation by strand
cds_span <- g[feat == "CDS" & parent != "", .(cds_lo = min(start), cds_hi = max(end)), by = parent]
ex <- merge(g[feat == "exon" & parent != "", .(chrom, start, end, strand, parent)], cds_span, by = "parent")
stopifnot(nrow(ex) > 100000)
left  <- ex[start < cds_lo, .(chrom, start, end = pmin(end, cds_lo), strand)]
right <- ex[end > cds_hi,   .(chrom, start = pmax(start, cds_hi), end, strand)]
utr5 <- merge_intervals(rbind(left[strand == "+"], right[strand == "-"]))
utr3 <- merge_intervals(rbind(right[strand == "+"], left[strand == "-"]))
other_exon <- subtract_intervals(exons, rbind(cds, utr5, utr3))
introns    <- subtract_intervals(genes, exons)
whole      <- data.table(chrom = genome$chrom, start = 0L, end = as.integer(genome$length))
intergenic <- subtract_intervals(whole, rbind(genes, promoters))

features <- list(promoter = promoters, UTR5 = utr5, CDS = cds, UTR3 = utr3, other_exon = other_exon,
                 intron = introns, intergenic = intergenic, CpG_island = cgi)
sizes <- rbindlist(lapply(names(features), function(nm) data.table(feature = nm,
  n_intervals = nrow(features[[nm]]), Mb = bp(features[[nm]]) / 1e6, pct_genome = 100 * bp(features[[nm]]) / GENOME_BP)))
stopifnot(all(sizes$Mb > 0.5))
fwrite(sizes, file.path(OUTDIR, "feature_sizes.tsv"), sep = "\t")

# overlap machinery
index_set <- function(bed) lapply(split(bed, bed$chrom), function(x) list(s = x$start, e = x$end))
set_idx <- lapply(features, index_set)
overlap_flag <- function(regions, tidx) {
  hits <- logical(nrow(regions)); rc <- regions$chrom; rs <- regions$start; re <- regions$end
  for (ch in unique(rc)) {
    k <- which(rc == ch); t <- tidx[[ch]]; if (is.null(t)) next
    j <- findInterval(re[k] - 1L, t$s); ok <- j >= 1L
    hits[k[ok]] <- t$e[j[ok]] > rs[k[ok]]
  }
  hits
}
flag_matrix <- function(regions) sapply(names(features), function(nm) overlap_flag(regions, set_idx[[nm]]))
primary_of <- function(fm) {
  p <- rep(NA_character_, nrow(fm))
  for (cat_ in rev(PRIORITY)) p[fm[, cat_]] <- cat_
  p[is.na(p)] <- "intergenic"
  p
}
count_cpg <- function(regions) {
  out <- integer(nrow(regions)); rc <- regions$chrom; rs <- regions$start; re <- regions$end
  for (ch in unique(rc)) {
    k <- which(rc == ch); p <- cpg_pos[[ch]]; if (is.null(p)) next
    out[k] <- findInterval(re[k], p) - findInterval(rs[k], p)
  }
  out
}
build_ci <- function(callable) {
  w <- as.numeric(callable$end - callable$start)
  list(chrom = callable$chrom, start = callable$start, cum = cumsum(w), total = sum(w))
}
draw_callable <- function(widths, ci) {
  u <- runif(length(widths)) * ci$total
  idx <- pmin(findInterval(u, ci$cum) + 1L, length(ci$start))
  st <- ci$start[idx] + as.integer(u - c(0, ci$cum)[idx])
  data.table(chrom = ci$chrom[idx], start = st, end = st + widths)
}
draw_callable_cpg <- function(widths, ci) {
  out <- draw_callable(widths, ci); bad <- which(count_cpg(out) < MIN_CPG); tries <- 0L
  while (length(bad) && tries < 500L) {
    d <- draw_callable(widths[bad], ci)
    out[bad, `:=`(chrom = d$chrom, start = d$start, end = d$end)]
    bad <- bad[count_cpg(out[bad]) < MIN_CPG]; tries <- tries + 1L
  }
  if (length(bad)) warning(sprintf("%d intervals could not meet MIN_CPG", length(bad)))
  out
}

# per comparison
assign_out <- list(); enr_out <- list(); dir_out <- list()
for (cmp in COMPS) {
  stats <- fread(cmd = paste("grep -v '^##'", file.path(SIGDIR, paste0(cmp, ".signature_stats.tsv"))), sep = "\t", header = TRUE)
  setnames(stats, "summary_label", "label")
  stats <- stats[chrom %in% chrom_names]
  ci  <- build_ci(merge_intervals(stats))
  dmr <- stats[label %in% c("HyperMethylated", "HypoMethylated"), .(chrom, start, end, label)]
  stopifnot(nrow(dmr) > 20)
  fm <- flag_matrix(dmr)
  dmr[, primary := primary_of(fm)]
  assign_out[[cmp]] <- cbind(data.table(comparison = cmp), dmr, as.data.table(fm))

  # direction contrast: hypo vs hyper, any-overlap per feature and primary category
  for (nm in c(names(features), paste0("primary:", PRIORITY))) {
    hit <- if (startsWith(nm, "primary:")) dmr$primary == sub("primary:", "", nm) else fm[, nm]
    hy <- dmr$label == "HypoMethylated"
    tab <- matrix(c(sum(hit & hy), sum(!hit & hy), sum(hit & !hy), sum(!hit & !hy)), 2)
    ft <- fisher.test(tab)
    dir_out[[length(dir_out) + 1]] <- data.table(comparison = cmp, feature = nm,
      hypo_hit = tab[1, 1], hypo_n = sum(hy), hyper_hit = tab[1, 2], hyper_n = sum(!hy),
      odds_ratio = unname(ft$estimate), p_value = ft$p.value)
  }

  # enrichment vs null: one draw per permutation, scored on every feature and primary category
  for (subname in c("all", "HypoMethylated", "HyperMethylated")) {
    reg <- if (subname == "all") dmr else dmr[label == subname]
    widths <- reg$end - reg$start
    ofm <- flag_matrix(reg); oprim <- primary_of(ofm)
    obs <- c(colMeans(ofm), setNames(sapply(PRIORITY, function(p) mean(oprim == p)), paste0("primary:", PRIORITY)))
    for (nullname in c("callable", "callable_cpg")) {
      drawfun <- if (nullname == "callable") draw_callable else draw_callable_cpg
      null <- t(replicate(N_PERM, {
        d <- drawfun(widths, ci); f <- flag_matrix(d); pr <- primary_of(f)
        c(colMeans(f), sapply(PRIORITY, function(p) mean(pr == p)))
      }))
      colnames(null) <- names(obs)
      for (nm in names(obs)) {
        nv <- null[, nm]
        enr_out[[length(enr_out) + 1]] <- data.table(comparison = cmp, subset = subname, feature = nm, null_model = nullname,
          n_regions = nrow(reg), n_hit = round(obs[[nm]] * nrow(reg)), pct_obs = 100 * obs[[nm]], pct_null = 100 * mean(nv),
          fold = obs[[nm]] / mean(nv),
          p_enriched = (sum(nv >= obs[[nm]]) + 1) / (N_PERM + 1),
          p_depleted = (sum(nv <= obs[[nm]]) + 1) / (N_PERM + 1))
      }
    }
  }
}

# per-comparison output
tag <- paste(COMPS, collapse = "+")
fwrite(rbindlist(assign_out), file.path(PERDIR, paste0(tag, "_assignment.tsv")), sep = "\t")
fwrite(rbindlist(enr_out),    file.path(PERDIR, paste0(tag, "_enrichment.tsv")), sep = "\t")
fwrite(rbindlist(dir_out),    file.path(PERDIR, paste0(tag, "_direction.tsv")), sep = "\t")
quit(save = "no")
}

# MODE=combine: pool comparisons, apply FDR, summarize
if (MODE == "combine") {
rd <- function(suffix) {
  f <- list.files(PERDIR, pattern = paste0("_", suffix, "\\.tsv$"), full.names = TRUE)
  rbindlist(lapply(f, fread))
}
assign_dt <- rd("assignment"); enr <- rd("enrichment"); dirs <- rd("direction")
stopifnot(uniqueN(enr$comparison) == 3, uniqueN(dirs$comparison) == 3)
enr[, p_two_sided := pmin(1, 2 * pmin(p_enriched, p_depleted))]
enr[, fdr := p.adjust(p_two_sided, method = "BH"), by = .(null_model, subset)]
dirs[, fdr := p.adjust(p_value, method = "BH"), by = comparison]
fwrite(assign_dt, file.path(OUTDIR, "dmr_gene_feature_assignment.tsv"), sep = "\t")
fwrite(enr,  file.path(OUTDIR, "dmr_gene_feature_enrichment.tsv"), sep = "\t")
fwrite(dirs, file.path(OUTDIR, "dmr_gene_feature_direction.tsv"), sep = "\t")

}
