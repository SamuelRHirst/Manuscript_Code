# Coastal ancestry of the two hybrids (four haplotypes) at differentiated autosomal SVs,
# by region on chromosomes 1 and 5 versus the autosomal background. Significance from a
# circular-shift permutation that rotates ancestry values along the genome-ordered SVs,
# preserving ancestry tract structure.

library(data.table)

N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000"))

svs <- fread("fixed_svs/diff_svs.tsv")
svs[, chr := sub(".*#", "", chrom)]
chr_sizes <- fread("fixed_svs/chr_sizes.txt")
svs <- svs[chr %in% chr_sizes[!grepl("^Chromosome_(W|Z)", chr), chr]]

svs <- svs[hybrid_AN == 4]
svs[, hybrid_DAC := round(hybrid_freq * 4)]
svs[, has_hybrid_derived := hybrid_DAC >= 1]
svs[, coastal_ancestry := ifelse(direction == "coastal_fixed", hybrid_DAC / 4, (4 - hybrid_DAC) / 4)]

CORE <- data.table(chr = c("Chromosome_01", "Chromosome_05"),
                   start = c(46.0e6, 21e6), end = c(59.5e6, 39e6))
PERIPHERY <- data.table(chr = c(rep("Chromosome_01", 4), "Chromosome_05"),
                        start = c(25.0e6, 59.5e6, 73.0e6, 86.5e6, 39e6),
                        end = c(46.0e6, 73.0e6, 86.5e6, 100.0e6, 60e6))

assign_region <- function(chr, pos) {
    region <- rep("background", length(chr))
    for (i in seq_len(nrow(PERIPHERY)))
        region[chr == PERIPHERY$chr[i] & pos >= PERIPHERY$start[i] & pos < PERIPHERY$end[i]] <- "periphery"
    for (i in seq_len(nrow(CORE)))
        region[chr == CORE$chr[i] & pos >= CORE$start[i] & pos < CORE$end[i]] <- "core"
    region
}
svs[, region := assign_region(chr, bedStart)]
bg <- svs[region == "background"]

summarise_tier <- function(sub, label) {
    data.table(label = label, n = nrow(sub),
               pct_hybrid_derived = 100 * mean(sub$has_hybrid_derived),
               p_derived_vs_bg = prop.test(c(sum(sub$has_hybrid_derived), sum(bg$has_hybrid_derived)),
                                           c(nrow(sub), nrow(bg)))$p.value,
               mean_coastal_ancestry = mean(sub$coastal_ancestry),
               pct_majority_coastal = 100 * mean(sub$coastal_ancestry > 0.5),
               pct_majority_inland = 100 * mean(sub$coastal_ancestry < 0.5),
               p_wilcox_vs_bg = wilcox.test(sub$coastal_ancestry, bg$coastal_ancestry)$p.value)
}

# Chromosomes 1 and 5 are kept separate: their ancestry biases run in opposite directions
make_tiers <- function(d) {
    tiers <- list()
    for (ch in c("Chromosome_01", "Chromosome_05")) {
        lbl <- sub("Chromosome_0", "Chr", ch)
        tiers[[paste(lbl, "core")]] <- d$chr == ch & d$region == "core"
        tiers[[paste(lbl, "periphery only")]] <- d$chr == ch & d$region == "periphery"
        tiers[[paste(lbl, "extension")]] <- d$chr == ch & d$region %in% c("core", "periphery")
    }
    tiers
}
tiers <- make_tiers(svs)

fwrite(rbindlist(Map(function(sel, lbl) summarise_tier(svs[sel], lbl), tiers, names(tiers))),
       "fixed_svs/hybrid_ancestry_stratified.tsv", sep = "\t")

# Circular-shift permutation
set.seed(42)
setorder(svs, chr, bedStart)
ca <- svs$coastal_ancestry
der <- as.numeric(svs$has_hybrid_derived)
bg_sel <- svs$region == "background"
rotate <- function(v, k) v[((seq_along(v) + k - 1L) %% length(v)) + 1L]

tiers <- make_tiers(svs)

perm_tier <- function(sel, label) {
    if (sum(sel) < 20) return(NULL)
    obs_ca <- mean(ca[sel]) - mean(ca[bg_sel])
    obs_der <- mean(der[sel]) - mean(der[bg_sel])
    ks <- sample.int(length(ca), N_PERM, replace = TRUE)
    d_ca <- vapply(ks, function(k) { r <- rotate(ca, k); mean(r[sel]) - mean(r[bg_sel]) }, numeric(1))
    d_der <- vapply(ks, function(k) { r <- rotate(der, k); mean(r[sel]) - mean(r[bg_sel]) }, numeric(1))
    data.table(label = label, n = sum(sel),
               delta_coastal_ancestry = obs_ca,
               p_perm_coastal_ancestry = (1 + sum(abs(d_ca) >= abs(obs_ca))) / (1 + N_PERM),
               null_lo_ca = quantile(d_ca, 0.025), null_hi_ca = quantile(d_ca, 0.975),
               delta_pct_derived = obs_der,
               p_perm_derived = (1 + sum(abs(d_der) >= abs(obs_der))) / (1 + N_PERM),
               null_lo_der = quantile(d_der, 0.025), null_hi_der = quantile(d_der, 0.975))
}

fwrite(rbindlist(Map(perm_tier, tiers, names(tiers))),
       "fixed_svs/hybrid_ancestry_permutation.tsv", sep = "\t")
