# Do hybrids show more read-level methylation disorder at TE copies? Per bird (hap1 + hap2 pooled), over copies
# with >= 5 informative reads: PDR = discordant reads / reads, and the share of intermediate reads. Each stratum is
# also contrasted with the same bird's old (> 15%) copies of the same class. Exact one-sided rank test with
# 2 hybrids among 10 birds (hybrids higher; minimum p = 1/45).

suppressPackageStartupMessages(library(data.table))
D <- "methylation/read_disorder"
MIN_READS <- 5L
TAXON <- c(AH14K04 = "hybrid", MF12K02 = "hybrid",
           AH19K01 = "coastal", AH19K05 = "coastal", JF08K03 = "coastal", JF08K05 = "coastal",
           AH19K17 = "inland", AH21K03 = "inland", KF14K01 = "inland", KF15K02 = "inland")
f <- list.files(file.path(D, "per_copy"), pattern = "\\.disorder\\.tsv$", full.names = TRUE)
dt <- rbindlist(lapply(f, fread))
dt <- dt[n_reads >= MIN_READS]
dt[, bird := sub("_hap[12]$", "", sample_hap)][, taxon := TAXON[bird]]
stopifnot(!anyNA(dt$taxon))
dt[, class := sub("/.*$", "", class_family)]
dt[, age := cut(divergence, c(-Inf, 5, 15, Inf), labels = c("young", "mid", "old"), right = FALSE)]

birds <- names(TAXON); is_h <- TAXON[birds] == "hybrid"; combos <- combn(10, 2)
exact_hi <- function(v) { v <- v[birds]; if (anyNA(v)) return(NA_real_)
  obs <- mean(v[is_h]); mean(apply(combos, 2, function(i) mean(v[i])) >= obs - 1e-12) }
bird_stats <- function(x) x[, .(pdr = sum(discordant) / sum(n_reads), inter = sum(intermediate) / sum(n_reads),
                                meth = sum(sum_frac) / sum(n_reads), copies = .N), by = bird]
summ <- function(sub, ref, label) {
  b <- bird_stats(sub); if (nrow(b) < 10 || min(b$copies) < 30) return(NULL)
  r <- bird_stats(ref)[, .(bird, pdr_ref = pdr, inter_ref = inter)]
  b <- merge(b, r, by = "bird")
  tm <- function(col, t) mean(b[TAXON[bird] == t][[col]])
  mid <- function(col) tm(col, "hybrid") - (tm(col, "coastal") + tm(col, "inland")) / 2
  v <- function(col) setNames(b[[col]], b$bird)
  b[, `:=`(c_pdr = pdr - pdr_ref, c_inter = inter - inter_ref)]
  data.table(stratum = label, min_copies = min(b$copies),
    pdr_coastal = tm("pdr", "coastal"), pdr_inland = tm("pdr", "inland"), pdr_hybrid = tm("pdr", "hybrid"),
    p_pdr = exact_hi(v("pdr")), inter_hyb_vs_mid = mid("inter"), p_inter = exact_hi(v("inter")),
    contrast_pdr_hyb_vs_mid = mid("c_pdr"), p_contrast_pdr = exact_hi(v("c_pdr")),
    contrast_inter_hyb_vs_mid = mid("c_inter"), p_contrast_inter = exact_hi(v("c_inter")),
    hybrid_ranks_pdr = paste(sort(rank(v("pdr"))[c("AH14K04", "MF12K02")]), collapse = ","))
}
res <- list()
for (a in levels(dt$age)) for (s in c("any", "fulllength")) {
  sub <- dt[age == a & (s == "any" | status == s)]
  res[[length(res) + 1]] <- summ(sub, dt[age == "old"], paste("ALL", a, s))
}
for (cl in c("LTR", "LINE", "DNA", "MITE")) for (a in levels(dt$age))
  res[[length(res) + 1]] <- summ(dt[class == cl & age == a], dt[class == cl & age == "old"], paste("class", cl, a))
fams <- dt[age == "young", .N, by = .(bird, class_family)][, .(m = min(N), nb = .N), by = class_family][nb == 10 & m >= 30, class_family]
for (fm in fams) res[[length(res) + 1]] <- summ(dt[class_family == fm & age == "young"],
                                                dt[class == sub("/.*$", "", fm) & age == "old"], paste("family", fm, "young"))
res <- rbindlist(res)
fwrite(res, file.path(D, "read_disorder_results.tsv"), sep = "\t")
