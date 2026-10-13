# Methylation discordance between haplotypes at homologous TE copies lifted to inland coordinates.
# Copies of the same family overlapping across haplotypes are clustered into loci. For each pair of
# haplotypes, over shared copies with >= 20 CpG observations in both:
#   mad = mean |difference in % methylation|
#   z2  = mean squared coverage-standardized difference, (p1 - p2) / sqrt(pbar (1 - pbar) (1/n1 + 1/n2))
# A hybrid's two haplotypes are compared with pairs of coastal and inland haplotypes from different pure
# birds, which differ through inherited (cis) sequence alone.

suppressPackageStartupMessages(library(data.table))
O  <- "methylation/hap_discordance"
TE <- "methylation/te_copies"
MIN_COV <- as.integer(Sys.getenv("MIN_COV", "20"))
COLS <- c("sample_hap", "status", "chrom", "start", "end", "divergence", "class_family", "sum_cov", "sum_mod", "n_cpg")
TAXON <- c(AH14K04 = "hybrid", MF12K02 = "hybrid",
           AH19K01 = "coastal", AH19K05 = "coastal", JF08K03 = "coastal", JF08K05 = "coastal",
           AH19K17 = "inland", AH21K03 = "inland", KF14K01 = "inland", KF15K02 = "inland")
HAPS <- as.vector(t(outer(names(TAXON), c("hap1", "hap2"), paste, sep = "_")))

# lifted copies joined back to their methylation records
all <- rbindlist(lapply(HAPS, function(h) {
  lf <- file.path(O, "lifted", paste0(h, ".te_inland.bed"))
  l <- fread(lf, header = FALSE, select = 1:4, col.names = c("ichrom", "istart", "iend", "name"))
  l[, `:=`(chrom = sub("_[0-9]+_[0-9]+$", "", name),
           start = as.integer(sub("^.*_([0-9]+)_[0-9]+$", "\\1", name)),
           end   = as.integer(sub("^.*_([0-9]+)$", "\\1", name)))]
  t <- rbind(fread(file.path(TE, paste0(h, ".fulllength.tsv")), header = FALSE, col.names = COLS),
             fread(file.path(TE, paste0(h, ".fragmented.tsv")), header = FALSE, col.names = COLS))
  m <- merge(l, t, by = c("chrom", "start", "end"))
  m <- unique(m, by = c("chrom", "start", "end"))
  m[, .(hap = h, ichrom, istart, iend, len_own = end - start, len_inl = iend - istart,
        status, divergence, class_family, sum_cov, sum_mod)]
}))
all <- all[len_inl > 0 & len_inl / len_own > 0.5 & len_inl / len_own < 2]  
all[, bird := sub("_hap[12]$", "", hap)][, taxon := TAXON[bird]]
all[, p := sum_mod / sum_cov]
all[, age := cut(divergence, c(-Inf, 5, 15, Inf), labels = c("young", "mid", "old"), right = FALSE)]

# cluster homologous copies into loci
setorder(all, class_family, ichrom, istart)
all[, locus := {
  brk <- c(TRUE, istart[-1] >= cummax(iend)[-.N]); cumsum(brk)
}, by = .(class_family, ichrom)]
all[, locus := .GRP, by = .(class_family, ichrom, locus)]
all <- all[order(-len_own)][!duplicated(all[order(-len_own), .(hap, locus)])]

# pairwise discordance
W <- all[sum_cov >= MIN_COV, .(hap, locus, p, n = sum_cov, age, class = sub("/.*$", "", class_family))]
setkey(W, hap)
pair_stats <- function(h1, h2, filt = NULL) {
  a <- W[h1]; b <- W[h2]
  m <- merge(a, b, by = "locus", suffixes = c("1", "2"))
  m[, `:=`(age = age1, class = class1)]
  if (!is.null(filt)) m <- m[eval(filt)]
  if (nrow(m) < 200) return(NULL)
  pb <- pmin(pmax((m$p1 * m$n1 + m$p2 * m$n2) / (m$n1 + m$n2), 0.01), 0.99)
  z2 <- (m$p1 - m$p2)^2 / (pb * (1 - pb) * (1 / m$n1 + 1 / m$n2))
  data.table(n_loci = nrow(m), mad = 100 * mean(abs(m$p1 - m$p2)), z2 = mean(z2), cor = cor(m$p1, m$p2))
}
b1 <- function(h) sub("_hap[12]$", "", h)
pairs <- CJ(h1 = HAPS, h2 = HAPS)[h1 < h2]
pairs[, `:=`(t1 = TAXON[b1(h1)], t2 = TAXON[b1(h2)], same_bird = b1(h1) == b1(h2))]
pairs[, type := fifelse(same_bird, paste0("within_", t1),
                fifelse(t1 == "hybrid" | t2 == "hybrid", "involves_hybrid",
                fifelse(t1 == t2, paste0("cross_", t1), "cross_subspecies")))]
pairs <- pairs[type != "involves_hybrid"]
strata <- list(all = NULL, young = quote(age == "young"), old = quote(age == "old"),
               LTR = quote(class == "LTR"), DNA = quote(class == "DNA"), MITE = quote(class == "MITE"), LINE = quote(class == "LINE"))
res <- rbindlist(lapply(names(strata), function(s) {
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    r <- pair_stats(pairs$h1[i], pairs$h2[i], strata[[s]])
    if (is.null(r)) NULL else cbind(pairs[i, .(h1, h2, type)], stratum = s, r)
  }))
}))
fwrite(res, file.path(O, "hap_discordance_pairs.tsv"), sep = "\t")
