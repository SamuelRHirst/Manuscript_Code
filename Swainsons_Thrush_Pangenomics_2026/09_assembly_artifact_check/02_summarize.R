# Summaries of the window scan:
#   A. assembly quality of the core blocks and flanks vs the autosomal background
#   B. differentiated-SV enrichment after removing SVs in suspect sequence
#   C. within the background, Spearman correlation of differentiated SVs per window with assembly quality

library(data.table)

dir <- "artifact_check"
win <- fread(file.path(dir, "windows.tsv"))
sv  <- fread(file.path(dir, "sv_sites.tsv"))

region_of <- function(chrom, pos) {
    mb <- pos / 1e6
    fifelse(chrom == "Chromosome_01" & mb >= 46.0 & mb < 59.5, "chr1_core",
    fifelse(chrom == "Chromosome_01" & mb >= 25.0 & mb < 100,  "chr1_flank",
    fifelse(chrom == "Chromosome_05" & mb >= 21   & mb < 39,   "chr5_core",
    fifelse(chrom == "Chromosome_05" & mb >= 21   & mb < 60,   "chr5_flank", "background"))))
}
regions <- c("chr1_core", "chr1_flank", "chr5_core", "chr5_flank", "background")

win[, region := region_of(chrom, (start + end) / 2)]
win[, len := end - start]
sv[, region := region_of(chrom, start)]
sv[, diff := fixed_class != "none"]

win[, `:=`(gap_frac = gap_bp / len, sat_frac = sat_bp / len, rm_frac = rm_bp / len,
           edta_frac = edta_bp / len, ends_per_Mb = contig_ends / (len / 1e6))]
svw <- sv[, .(sv_sites = .N,
              sv_ineligible_AN = sum(eligible == 0 & subtype != "SV_Complex"),
              n_diff = sum(diff)),
          by = .(chrom, wstart = as.integer(floor(start / 1e5) * 1e5))]
win[, start := as.integer(start)]
win <- merge(win, svw, by.x = c("chrom", "start"), by.y = c("chrom", "wstart"), all.x = TRUE)
for (v in c("sv_sites", "sv_ineligible_AN", "n_diff")) set(win, which(is.na(win[[v]])), v, 0L)
win[, blind_frac := fifelse(sv_sites > 0, sv_ineligible_AN / sv_sites, NA_real_)]

metrics <- win[, .(windows = .N,
                   haps_covered50_median = as.numeric(median(haps_covered50)),
                   pct_windows_lt16_haps = 100 * mean(haps_covered50 < 16),
                   mean_hap_cov = mean(mean_hap_cov),
                   contig_ends_per_Mb = sum(contig_ends) / (sum(len) / 1e6),
                   pct_gap = 100 * sum(gap_bp) / sum(len),
                   pct_satellite = 100 * sum(sat_bp) / sum(len),
                   pct_RM = 100 * sum(rm_bp) / sum(len),
                   pct_EDTA = 100 * sum(edta_bp) / sum(len),
                   odgi_depth_median = median(odgi_depth, na.rm = TRUE),
                   pct_SV_sites_failing_AN = 100 * sum(sv_ineligible_AN) / sum(sv_sites)),
               by = region][match(regions, region)][!is.na(region)]
fwrite(metrics, file.path(dir, "A_assembly_metrics_by_region.tsv"), sep = "\t")

sv[, wstart := as.integer(floor(start / 1e5) * 1e5)]
sv[, in_weak_window := paste(chrom, wstart) %in% win[haps_covered50 < 16, paste(chrom, start)]]
filters <- list(
    none               = rep(TRUE, nrow(sv)),
    no_contig_end_10kb = sv$contig_ends_10kb == 0,
    lt3_ends_50kb      = sv$contig_ends_50kb < 3,
    no_sat_gap_1kb     = sv$sat_1kb == 0 & sv$gap_1kb == 0,
    ge16_haps_window   = !sv$in_weak_window,
    all_combined       = sv$contig_ends_10kb == 0 & sv$sat_1kb == 0 & sv$gap_1kb == 0 & !sv$in_weak_window)

enrich <- rbindlist(lapply(names(filters), function(f) {
    Mb <- if (f %in% c("ge16_haps_window", "all_combined"))
        win[haps_covered50 >= 16, .(Mb = sum(len) / 1e6), by = region]
    else win[, .(Mb = sum(len) / 1e6), by = region]
    r <- merge(sv[filters[[f]], .(diff = sum(diff), eligible = sum(eligible == 1)), by = region], Mb, by = "region")
    r[, `:=`(filter = f, diff_per_Mb = diff / Mb, diff_share_pct = 100 * diff / eligible)]
    bg <- r[region == "background"]
    r[, `:=`(fold_per_Mb = diff_per_Mb / bg$diff_per_Mb, fold_per_SV = diff_share_pct / bg$diff_share_pct)]
}))
enrich <- merge(enrich, enrich[filter == "none", .(region, diff0 = diff)], by = "region")
enrich[, pct_diff_retained := 100 * diff / diff0][, diff0 := NULL]
setcolorder(enrich, c("filter", "region"))
enrich <- enrich[order(match(filter, names(filters)), match(region, regions))]
fwrite(enrich, file.path(dir, "B_enrichment_under_filters.tsv"), sep = "\t")

bgw <- win[region == "background" & sv_sites > 0]
vars <- c("haps_covered50", "mean_hap_cov", "ends_per_Mb", "gap_frac", "sat_frac",
          "rm_frac", "edta_frac", "odgi_depth", "blind_frac")
cors <- rbindlist(lapply(vars, function(v) {
    x <- bgw[[v]]; ok <- is.finite(x)
    if (sum(ok) < 10 || sd(x[ok]) == 0) return(data.table(metric = v, rho = NA_real_, p = NA_real_, n = sum(ok)))
    ct <- suppressWarnings(cor.test(x[ok], bgw$n_diff[ok], method = "spearman"))
    data.table(metric = v, rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}))
fwrite(cors, file.path(dir, "C_background_correlations.tsv"), sep = "\t")
