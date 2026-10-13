# Do hybrids lose methylation at young TE copies? Copies are split by divergence from consensus
# (< 5%, 5-15%, > 15%), class and family. Per individual, methylation is the mean over copies.
# With 2 hybrids among 10 birds, p is an exact one-sided rank test over the 45 ways of labelling two birds
# as hybrids (minimum 1/45), for (a) methylation in the stratum and (b) the within-bird contrast with that
# bird's old (> 15%) copies of the same class.

suppressPackageStartupMessages(library(data.table))
D <- "methylation/te_copies"
COLS <- c("sample_hap", "status", "chrom", "start", "end", "divergence", "class_family", "sum_cov", "sum_mod", "n_cpg")
MIN_COPIES <- 30L

files <- list.files(D, pattern = "\\.(fulllength|fragmented)\\.tsv$", full.names = TRUE)
dt <- rbindlist(lapply(files, fread, header = FALSE, col.names = COLS))
TAXON <- c(AH14K04 = "hybrid", MF12K02 = "hybrid",
           AH19K01 = "coastal", AH19K05 = "coastal", JF08K03 = "coastal", JF08K05 = "coastal",
           AH19K17 = "inland", AH21K03 = "inland", KF14K01 = "inland", KF15K02 = "inland")
dt[, bird := sub("_hap[12]$", "", sample_hap)]
dt[, taxon := TAXON[bird]]
stopifnot(!anyNA(dt$taxon))
dt[, pct := 100 * sum_mod / sum_cov]
dt[, class := sub("/.*$", "", class_family)]
dt[, age := cut(divergence, c(-Inf, 5, 15, Inf), labels = c("young_lt5", "mid_5to15", "old_gt15"), right = FALSE)]

birds <- names(TAXON); is_h <- TAXON[birds] == "hybrid"
combos <- combn(length(birds), 2)
exact_p <- function(v) {
  v <- v[birds]; if (anyNA(v)) return(NA_real_)
  obs <- mean(v[is_h]); null <- apply(combos, 2, function(i) mean(v[i]))
  mean(null <= obs + 1e-12)
}
summarise_stratum <- function(sub, ref, label) {
  b <- sub[, .(m = mean(pct), n = .N), by = bird]
  if (nrow(b) < 10 || min(b$n) < MIN_COPIES) return(NULL)
  r <- ref[, .(m_ref = mean(pct)), by = bird]
  b <- merge(b, r, by = "bird", all.x = TRUE); b[, contrast := m - m_ref]
  v <- setNames(b$m, b$bird); cv <- setNames(b$contrast, b$bird)
  tm <- function(x, t) mean(x[TAXON[names(x)] == t])
  data.table(stratum = label, min_copies_per_bird = min(b$n),
    coastal = tm(v, "coastal"), inland = tm(v, "inland"), hybrid = tm(v, "hybrid"),
    hyb_minus_midparent = tm(v, "hybrid") - (tm(v, "coastal") + tm(v, "inland")) / 2,
    p_level = exact_p(v),
    contrast_hyb_minus_midparent = tm(cv, "hybrid") - (tm(cv, "coastal") + tm(cv, "inland")) / 2,
    p_contrast = exact_p(cv),
    hybrids_rank_level = paste(sort(rank(v)[names(v) %in% names(TAXON)[TAXON == "hybrid"]]), collapse = ","))
}

res <- list()
# all TEs by age and status
for (a in levels(dt$age)) for (s in c("fulllength", "fragmented", "any")) {
  sub <- dt[age == a & (s == "any" | status == s)]
  res[[length(res) + 1]] <- summarise_stratum(sub, dt[age == "old_gt15"], paste("ALL", a, s))
}
# by class and age
for (cl in c("LTR", "LINE", "DNA", "MITE", "SINE", "RC")) for (a in levels(dt$age)) {
  res[[length(res) + 1]] <- summarise_stratum(dt[class == cl & age == a], dt[class == cl & age == "old_gt15"],
                                              paste("class", cl, a))
}
# young copies by family, for families with enough young copies in every bird
fams <- dt[age == "young_lt5", .N, by = .(bird, class_family)][, .(minN = min(N), nb = .N), by = class_family][nb == 10 & minN >= MIN_COPIES, class_family]
for (f in fams) {
  cl <- sub("/.*$", "", f)
  res[[length(res) + 1]] <- summarise_stratum(dt[class_family == f & age == "young_lt5"],
                                              dt[class == cl & age == "old_gt15"], paste("family", f, "young_lt5"))
  res[[length(res) + 1]] <- summarise_stratum(dt[class_family == f & age == "young_lt5" & status == "fulllength"],
                                              dt[class == cl & age == "old_gt15"], paste("family", f, "young_lt5 fulllength"))
}
res <- rbindlist(res)
res[, `:=`(fdr_level = p.adjust(p_level, "BH"), fdr_contrast = p.adjust(p_contrast, "BH"))]
fwrite(res, file.path(D, "young_te_methylation.tsv"), sep = "\t")
