# Density of differentiated SVs in the chromosome 1 and 5 regions versus the autosomal background.
#   core:      ancestry blocks of Blain et al. 2025
#   periphery: rest of the Justen et al. 2025 zone around each core (flanking region)
#   extension: core + periphery (full region)
# Significance from a circular-shift permutation that rotates the observed SV positions
# around the concatenated autosomes, preserving their clustering.

library(data.table)

N_PERM <- as.integer(Sys.getenv("N_PERM", unset = "10000"))
set.seed(42)

svs <- fread("fixed_svs/diff_svs.tsv")
svs[, chr := sub(".*#", "", chrom)]
chr_sizes <- fread("fixed_svs/chr_sizes.txt")
AUTOSOMES <- chr_sizes[!grepl("^Chromosome_(W|Z)", chr), chr]
svs <- svs[chr %in% AUTOSOMES]

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

core_len <- sum(CORE$end - CORE$start)
peri_len <- sum(PERIPHERY$end - PERIPHERY$start)
auto_len <- sum(chr_sizes[chr %in% AUTOSOMES, length])
lens <- c(core = core_len, periphery = peri_len, extension = core_len + peri_len,
          background = auto_len - core_len - peri_len)

density_block <- function(sub, label) {
    n <- c(core = sum(sub$region == "core"), periphery = sum(sub$region == "periphery"),
           background = sum(sub$region == "background"))
    n["extension"] <- n["core"] + n["periphery"]
    dens <- n / (lens[names(n)] / 1e6)
    pt <- function(a, b) poisson.test(c(n[a], n[b]), c(lens[a], lens[b]))$p.value
    data.table(comparison = label,
               n_core = n["core"], n_periphery = n["periphery"],
               n_extension = n["extension"], n_background = n["background"],
               dens_core = dens["core"], dens_periphery = dens["periphery"],
               dens_extension = dens["extension"], dens_background = dens["background"],
               fold_core_vs_bg = dens["core"] / dens["background"],
               fold_ext_vs_bg = dens["extension"] / dens["background"],
               fold_core_vs_periphery = dens["core"] / dens["periphery"],
               p_core_vs_bg = pt("core", "background"),
               p_ext_vs_bg = pt("extension", "background"),
               p_core_vs_periphery = pt("core", "periphery"))
}

subsets <- c(list(Overall = svs), split(svs, by = "fixed_class"), split(svs, by = "direction"))
for (fc in unique(svs$fixed_class))
    for (dr in unique(svs$direction))
        subsets[[paste(fc, dr, sep = "_x_")]] <- svs[fixed_class == fc & direction == dr]
fwrite(rbindlist(Map(density_block, subsets, names(subsets))),
       "fixed_svs/sv_density_stratified.tsv", sep = "\t")

bg_dens <- sum(svs$region == "background") / (lens["background"] / 1e6)
windows <- rbind(CORE[, tier := "core"], PERIPHERY[, tier := "periphery"])
per_window <- windows[, {
    n <- svs[chr == .BY$chr & bedStart >= start & bedStart < end, .N]
    len <- end - start
    .(start_mb = start / 1e6, end_mb = end / 1e6, length_mb = len / 1e6, n_svs = n,
      svs_per_mb = n / (len / 1e6), fold_vs_background = n / (len / 1e6) / bg_dens)
}, by = .(tier, chr, start, end)][, !c("start", "end")]
fwrite(per_window, "fixed_svs/sv_density_per_window.tsv", sep = "\t")

# Circular-shift permutation on the concatenated autosomes
csa <- chr_sizes[chr %in% AUTOSOMES][order(chr)]
csa[, offset := cumsum(c(0, head(length, -1)))]
G <- sum(csa$length)
off <- setNames(csa$offset, csa$chr)
svs[, gpos := bedStart + off[chr]]

brk <- sort(unique(c(0, CORE$start + off[CORE$chr], CORE$end + off[CORE$chr],
                     PERIPHERY$start + off[PERIPHERY$chr], PERIPHERY$end + off[PERIPHERY$chr], G)))
mids <- (head(brk, -1) + tail(brk, -1)) / 2
seg_region <- vapply(mids, function(m) {
    if (any(m >= CORE$start + off[CORE$chr] & m < CORE$end + off[CORE$chr])) return("core")
    if (any(m >= PERIPHERY$start + off[PERIPHERY$chr] & m < PERIPHERY$end + off[PERIPHERY$chr])) return("periphery")
    "background"
}, character(1))
classify <- function(gp) seg_region[findInterval(gp, brk, rightmost.closed = TRUE)]
stopifnot(identical(classify(svs$gpos), svs$region))

folds <- function(gp) {
    tb <- table(factor(classify(gp), levels = c("core", "periphery", "background")))
    d <- unname(as.numeric(tb) / (lens[c("core", "periphery", "background")] / 1e6))
    ext <- unname((tb[["core"]] + tb[["periphery"]]) / (lens["extension"] / 1e6))
    c(core_vs_bg = d[1] / d[3], ext_vs_bg = ext / d[3], core_vs_pe = d[1] / d[2])
}

perm_block <- function(sub, label) {
    if (nrow(sub) < 20) return(NULL)
    obs <- folds(sub$gpos)
    null <- t(vapply(runif(N_PERM, 0, G), function(s) folds((sub$gpos + s) %% G), obs))
    res <- data.table(comparison = label, n = nrow(sub))
    for (k in names(obs)) {
        nd <- null[, k][is.finite(null[, k])]
        res[[paste0("fold_", k)]] <- obs[[k]]
        res[[paste0("nullmed_", k)]] <- median(nd)
        res[[paste0("nullhi_", k)]] <- quantile(nd, 0.975)
        res[[paste0("pperm_", k)]] <- (1 + sum(nd >= obs[[k]])) / (1 + length(nd))
    }
    res
}

perm_subsets <- c(list(Overall = svs), split(svs, by = "fixed_class"), split(svs, by = "direction"))
fwrite(rbindlist(Map(perm_block, perm_subsets, names(perm_subsets))),
       "fixed_svs/sv_density_permutation.tsv", sep = "\t")
