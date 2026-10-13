# Do the chromosome 1 and 5 blocks depart from what their recombination rate predicts?
# Window statistics are averaged to 1 Mb tiles. The background expectation is a GAM fitted to the 500
# background tiles (log y ~ s(log10 r); sensitivity model adds s(log10 SNP density)), and block excess is
# the mean residual of its tiles. Uncertainty: bootstrap of background tiles (refitting the GAM) with a
# moving-block bootstrap of block tiles (runs of 4), 2,000 replicates.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
P <- "arg/post"
set.seed(20260921)
NBOOT <- 2000; BLOCK_LEN <- 4

w <- fread(file.path(P, "window_summary.tsv"))
tl <- w[, .(Tb = mean(T_between), Tw = mean(T_within), mono = mean(mono_recip), r = r[1],
            dens = mean(sites_per_kb), region = names(sort(table(region), decreasing = TRUE))[1],
            chrom = chrom[1], start = min(win_start), nwin = .N), by = tile]
tl[, `:=`(lr = log10(r), ld = log10(dens))]
bg <- tl[region == "background"]
blocks <- list(chr1 = tl[region == "core_01"][order(start)], chr5 = tl[region == "core_05"][order(start)])
stopifnot(all(is.finite(tl$lr)))
rng_bg <- range(bg$lr)
for (b in names(blocks)) stopifnot(all(blocks[[b]]$lr >= rng_bg[1] & blocks[[b]]$lr <= rng_bg[2]))

resp <- list(
  T_between = list(col = "Tb",   log = TRUE,  fam = gaussian()),
  T_within  = list(col = "Tw",   log = TRUE,  fam = gaussian()),
  monophyly = list(col = "mono", log = FALSE, fam = quasibinomial()))
fit <- function(d, rs, sens = FALSE) {
  d <- copy(d); d[, y := if (rs$log) log(get(rs$col)) else get(rs$col)]
  f <- if (sens) y ~ s(lr, k = 6) + s(ld, k = 6) else y ~ s(lr, k = 6)
  gam(f, data = d, family = rs$fam, method = "REML")
}
excess <- function(m, blk, rs) {
  pred <- predict(m, newdata = blk, type = "response")
  obs  <- blk[[rs$col]]
  if (rs$log) mean(log(obs) - pred) else mean(obs - pred)
}
mbb <- function(n, L) {
  starts <- sample.int(n - L + 1, ceiling(n / L), replace = TRUE)
  idx <- unlist(lapply(starts, function(s) s:(s + L - 1)))
  idx[seq_len(n)]
}
fmt <- function(rs, x) if (rs$log) sprintf("%+.1f%%", 100 * (exp(x) - 1)) else sprintf("%+.3f", x)

out <- file(file.path(P, "arg_test_results.txt"), "w")
say <- function(...) { s <- sprintf(...); cat(s, "\n"); writeLines(s, out) }
say("Tiles: %d background, %d chr1 block, %d chr5 block. Bootstrap: %d replicates, block length %d tiles.",
    nrow(bg), nrow(blocks$chr1), nrow(blocks$chr5), NBOOT, BLOCK_LEN)

tile_out <- list(); curve_out <- list()
for (nm in names(resp)) {
  rs <- resp[[nm]]
  m  <- fit(bg, rs); ms <- fit(bg, rs, sens = TRUE)
  say("\n==== %s ====  (primary model: deviance explained %.1f%%)", nm, 100 * summary(m)$dev.expl)
  for (b in names(blocks)) {
    blk <- blocks[[b]]
    res <- if (rs$log) log(blk[[rs$col]]) - predict(m, blk, type = "response") else blk[[rs$col]] - predict(m, blk, type = "response")
    ac  <- if (length(res) > 2) cor(res[-1], res[-length(res)]) else NA
    est <- excess(m, blk, rs); est_s <- excess(ms, blk, rs)
    boot <- replicate(NBOOT, {
      mb <- fit(bg[sample.int(nrow(bg), replace = TRUE)], rs)
      excess(mb, blk[mbb(nrow(blk), BLOCK_LEN)], rs)
    })
    ci <- quantile(boot, c(.025, .975)); p <- (sum(boot <= 0) + 1) / (NBOOT + 1)
    say("  %s block: excess %s  [95%% CI %s to %s]  one-sided p = %.4f | lag-1 residual autocorrelation %.2f | with density covariate: %s",
        b, fmt(rs, est), fmt(rs, ci[1]), fmt(rs, ci[2]), p, ac, fmt(rs, est_s))
    tile_out[[paste(nm, b)]] <- data.table(stat = nm, block = b, tile = blk$tile, residual = res)
  }
  # model-free comparison with background tiles matched on r and density
  core <- rbind(blocks$chr1, blocks$chr5)
  mt <- bg[r >= quantile(core$r, .05) & r <= quantile(core$r, .95) & dens >= quantile(core$dens, .05) & dens <= quantile(core$dens, .95)]
  say("  model-free check: %d matched background tiles, median %s = %.4g; chr1 block median %.4g; chr5 block median %.4g",
      nrow(mt), nm, median(mt[[rs$col]]), median(blocks$chr1[[rs$col]]), median(blocks$chr5[[rs$col]]))
  # background curve with 95% band
  g <- data.table(lr = seq(rng_bg[1], rng_bg[2], length.out = 200))
  pr <- predict(m, g, type = "link", se.fit = TRUE)
  inv <- if (rs$log) exp else function(x) m$family$linkinv(x)
  curve_out[[nm]] <- data.table(stat = nm, r = 10^g$lr, fit = inv(pr$fit), lo = inv(pr$fit - 1.96 * pr$se.fit), hi = inv(pr$fit + 1.96 * pr$se.fit))
}
close(out)
fwrite(rbindlist(tile_out), file.path(P, "arg_test_tiles.tsv"), sep = "\t")
fwrite(rbindlist(curve_out), file.path(P, "arg_test_curves.tsv"), sep = "\t")
cat("\nwrote arg_test_results.txt, arg_test_tiles.tsv, arg_test_curves.tsv\n")
