# Could a large inversion in the chromosome 1 or 5 blocks be hidden by contigs breaking at its breakpoints?
# From the RagTag minimap2 alignments of each pseudo-haplotype (primary, MAPQ >= 5, blocks >= 10 kb):
#   inverted segments: blocks on the strand opposite to the rest of their contig, merged within 100 kb
#   detectability: per 100 kb, pseudo-haplotypes whose contig runs collinearly >= 250 kb on both sides
# Positive control: the 4.8 Mb chromosome 4 inversion.

library(data.table)

outdir <- "artifact_check/orientation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
RAGTAG <- "inversions/ragtag"

MIN_BLOCK <- 10000
MIN_MAPQ  <- 5
MERGE_GAP <- 100000
FLANK     <- 250000
MIN_SEG   <- 50000

regions <- data.table(region = c("chr1_full", "chr5_full", "chr4_control"),
                      chrom  = c("Chromosome_01", "Chromosome_05", "Chromosome_04"),
                      start  = c(25e6, 21e6, 70e6), end = c(100e6, 60e6, 80e6))
cores <- data.table(chrom = c("Chromosome_01", "Chromosome_05", "Chromosome_04"),
                    cstart = c(46e6, 21e6, 72.19e6), cend = c(59.5e6, 39e6, 77.03e6),
                    label = c("chr1 core", "chr5 core", "chr4 SyRI inversion"))
birds <- data.table(bird = c("AH19K17", "AH21K03", "KF14K01", "KF15K02", "AH19K01", "AH19K05", "JF08K03", "JF08K05",
                             "AH14K04", "MF12K02"),
                    group = rep(c("inland", "coastal", "hybrid"), c(4, 4, 2)))
haps <- birds[, .(hap = paste0(bird, "_hap", 1:2)), by = .(bird, group)]

paf <- rbindlist(lapply(haps$hap, function(h) {
  f <- file.path(RAGTAG, h, "ragtag.scaffold.asm.paf")
  x <- fread(f, header = FALSE, sep = "\t", fill = TRUE, select = 1:13,
             col.names = c("qname", "qlen", "qs", "qe", "strand", "chrom", "tlen", "ts", "te", "nmatch", "alen",
                           "mapq", "tp"))
  x[, hap := h]
}))
paf <- paf[tp == "tp:A:P" & mapq >= MIN_MAPQ & chrom %in% regions$chrom]
paf <- merge(paf, haps, by = "hap")
paf[, contig := paste(hap, qname)]

# majority strand of each contig on each chromosome, by aligned target bp (blocks >= MIN_BLOCK)
big <- paf[te - ts >= MIN_BLOCK]
maj <- big[, .(bp = sum(te - ts)), by = .(contig, chrom, strand)][order(-bp)][, .SD[1], by = .(contig, chrom)]
big <- merge(big, maj[, .(contig, chrom, maj_strand = strand)], by = c("contig", "chrom"))
big[, orient := fifelse(strand == maj_strand, "collinear", "inverted")]

# inverted segments
segs <- big[orient == "inverted"][order(contig, chrom, ts)]
segs[, grp := cumsum(c(1, diff(ts) > 0 & (ts[-1] - cummax(te)[-.N]) > MERGE_GAP)), by = .(contig, chrom)]
segs <- segs[, .(hap = hap[1], bird = bird[1], group = group[1], start = min(ts), end = max(te),
                 aligned_bp = sum(te - ts), n_blocks = .N), by = .(contig, chrom, grp)]
col <- big[orient == "collinear", .(contig, chrom, cs = ts, ce = te)]
segs[, `:=`(left_flank = mapply(function(cn, ch, s) any(col$contig == cn & col$chrom == ch & col$ce <= s + 1e4), contig, chrom, start),
            right_flank = mapply(function(cn, ch, e) any(col$contig == cn & col$chrom == ch & col$cs >= e - 1e4), contig, chrom, end))]
segs[, type := fifelse(left_flank & right_flank, "internal", "terminal")]
segs[, len := end - start]
segs <- segs[len >= MIN_SEG]
segs <- merge(segs, regions[, .(region, chrom, rstart = start, rend = end)], by = "chrom",
              allow.cartesian = TRUE)[start < rend & end > rstart]
fwrite(segs[order(region, start), .(region, chrom, hap, bird, group, start, end, len, aligned_bp, n_blocks, type)],
       file.path(outdir, "inverted_segments.tsv"), sep = "\t")

# cluster segments across haplotypes (single-linkage on >= 50% reciprocal overlap) and count carriers
clus <- segs[order(chrom, start)]
if (nrow(clus)) {
  clus[, id := .I]; parent <- clus$id
  find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
  for (i in seq_len(nrow(clus))) for (j in seq_len(nrow(clus))) if (i < j && clus$chrom[i] == clus$chrom[j]) {
    ov <- min(clus$end[i], clus$end[j]) - max(clus$start[i], clus$start[j])
    if (ov > 0 && ov / clus$len[i] >= 0.5 && ov / clus$len[j] >= 0.5) parent[find(j)] <- find(i)
  }
  clus[, cluster := sapply(id, find)]
  clsum <- clus[, .(start = as.numeric(median(start)), end = as.numeric(median(end)), len_median = as.numeric(median(len)),
                    n_haps = uniqueN(hap), n_internal = sum(type == "internal"),
                    inland_haps = uniqueN(hap[group == "inland"]), coastal_haps = uniqueN(hap[group == "coastal"]),
                    hybrid_haps = uniqueN(hap[group == "hybrid"])), by = .(region, chrom, cluster)][order(region, start)]
  clsum[, cluster := NULL]
} else clsum <- data.table()
fwrite(clsum, file.path(outdir, "inversion_clusters.tsv"), sep = "\t")

# breakpoint detectability
spans <- big[orient == "collinear", .(s = min(ts), e = max(te)), by = .(contig, hap, chrom)]
grid <- regions[, .(pos = seq(start, end - 1e5, by = 1e5) + 5e4), by = .(region, chrom)]
grid[, haps_spanning := mapply(function(ch, p) uniqueN(spans[chrom == ch & s <= p - FLANK & e >= p + FLANK, hap]), chrom, pos)]
grid <- merge(grid, cores[, .(chrom, cstart, cend)], by = "chrom")
grid[, in_core := pos >= cstart & pos < cend]
fwrite(grid[, .(region, chrom, pos, haps_spanning, in_core)], file.path(outdir, "breakpoint_detectability.tsv"), sep = "\t")
pw <- grid[, .(positions = .N, median_haps = as.numeric(median(haps_spanning)),
               pct_ge16 = 100 * mean(haps_spanning >= 16), pct_ge10 = 100 * mean(haps_spanning >= 10),
               pct_lt4 = 100 * mean(haps_spanning < 4)), by = .(region, in_core)][order(region, -in_core)]
fwrite(pw, file.path(outdir, "breakpoint_detectability_summary.tsv"), sep = "\t")
