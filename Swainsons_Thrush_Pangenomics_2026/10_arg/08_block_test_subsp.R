# As 07_block_test.R, testing within-subspecies coalescence for coastal and inland separately and adding
# bootstrap intervals for the SNP-density sensitivity model.

suppressPackageStartupMessages({ library(data.table); library(mgcv) })
P <- "arg/post"
set.seed(20260921)
NBOOT <- 2000; BLOCK_LEN <- 4

w <- fread(file.path(P, "window_summary.tsv"))
tl <- w[, .(Tb = mean(T_between), Tw = mean(T_within), Tco = mean(T_within_co), Tin = mean(T_within_in),
            mono = mean(mono_recip), r = r[1], dens = mean(sites_per_kb),
            region = names(sort(table(region), decreasing = TRUE))[1],
            chrom = chrom[1], start = min(win_start), nwin = .N), by = tile]
tl[, `:=`(lr = log10(r), ld = log10(dens))]
bg <- tl[region == "background"]
blocks <- list(chr1 = tl[region == "core_01"][order(start)], chr5 = tl[region == "core_05"][order(start)])

resp <- list(
  T_between      = list(col = "Tb",  log = TRUE, fam = gaussian()),
  T_within       = list(col = "Tw",  log = TRUE, fam = gaussian()),
  T_within_coast = list(col = "Tco", log = TRUE, fam = gaussian()),
  T_within_inl   = list(col = "Tin", log = TRUE, fam = gaussian()))
fit <- function(d, rs, sens = FALSE) {
  d <- copy(d); d[, y := log(get(rs$col))]
  f <- if (sens) y ~ s(lr, k = 6) + s(ld, k = 6) else y ~ s(lr, k = 6)
  gam(f, data = d, family = rs$fam, method = "REML")
}
excess <- function(m, blk, rs) mean(log(blk[[rs$col]]) - predict(m, newdata = blk, type = "response"))
mbb <- function(n, L) {
  starts <- sample.int(n - L + 1, ceiling(n / L), replace = TRUE)
  unlist(lapply(starts, function(s) s:(s + L - 1)))[seq_len(n)]
}
pct <- function(x) sprintf("%+.1f%%", 100 * (exp(x) - 1))

out <- file(file.path(P, "arg_test_results_subsp.txt"), "w")
say <- function(...) { s <- sprintf(...); cat(s, "\n"); writeLines(s, out) }
say("Tiles: %d background, %d chr1 block, %d chr5 block. Bootstrap: %d replicates, block length %d tiles.",
    nrow(bg), nrow(blocks$chr1), nrow(blocks$chr5), NBOOT, BLOCK_LEN)
say("p = one-sided share of bootstrap replicates on the other side of zero from the estimate.")

tile_out <- list(); curve_out <- list()
for (nm in names(resp)) {
  rs <- resp[[nm]]
  for (sens in c(FALSE, TRUE)) {
    m <- fit(bg, rs, sens)
    say("\n==== %s | %s ====  (deviance explained %.1f%%)", nm,
        if (sens) "sensitivity: + SNP density" else "primary: recombination only", 100 * summary(m)$dev.expl)
    for (b in names(blocks)) {
      blk <- blocks[[b]]
      est <- excess(m, blk, rs)
      boot <- replicate(NBOOT, {
        mb <- fit(bg[sample.int(nrow(bg), replace = TRUE)], rs, sens)
        excess(mb, blk[mbb(nrow(blk), BLOCK_LEN)], rs)
      })
      ci <- quantile(boot, c(.025, .975))
      p <- if (est >= 0) (sum(boot <= 0) + 1) / (NBOOT + 1) else (sum(boot >= 0) + 1) / (NBOOT + 1)
      say("  %s block: excess %s  [95%% CI %s to %s]  one-sided p = %.4f", b, pct(est), pct(ci[1]), pct(ci[2]), p)
      if (!sens) tile_out[[paste(nm, b)]] <- data.table(stat = nm, block = b, tile = blk$tile,
                                                         residual = log(blk[[rs$col]]) - predict(m, blk, type = "response"))
    }
    if (!sens) {
      core <- rbind(blocks$chr1, blocks$chr5)
      mt <- bg[r >= quantile(core$r, .05) & r <= quantile(core$r, .95) & dens >= quantile(core$dens, .05) & dens <= quantile(core$dens, .95)]
      say("  model-free check: %d matched background tiles, median %.4g; chr1 block median %.4g; chr5 block median %.4g",
          nrow(mt), median(mt[[rs$col]]), median(blocks$chr1[[rs$col]]), median(blocks$chr5[[rs$col]]))
      g <- data.table(lr = seq(min(bg$lr), max(bg$lr), length.out = 200))
      pr <- predict(m, g, type = "link", se.fit = TRUE)
      curve_out[[nm]] <- data.table(stat = nm, r = 10^g$lr, fit = exp(pr$fit), lo = exp(pr$fit - 1.96 * pr$se.fit),
                                    hi = exp(pr$fit + 1.96 * pr$se.fit))
    }
  }
}
close(out)
fwrite(rbindlist(tile_out), file.path(P, "arg_test_tiles_subsp.tsv"), sep = "\t")
fwrite(rbindlist(curve_out), file.path(P, "arg_test_curves_subsp.tsv"), sep = "\t")
cat("\nwrote arg_test_results_subsp.txt, arg_test_tiles_subsp.tsv, arg_test_curves_subsp.tsv\n")
