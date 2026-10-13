# Swainson's Thrush Pangenomics

Code for the Swainson's thrush (*Catharus ustulatus*) structural variant pangenome.

Scripts are numbered in the order they were run. Each is run from the project root, which is expected to contain:

```
reads/       <sample>.hifi_reads.bam and .fastq.gz
refs/        inland and coastal reference assemblies and annotation, C. minimus assembly, LDhat map
```

All other directories (`assembly/`, `pangenome/`, `repeats/`, `inversions/`, ...) are created by the scripts.
Per-sample steps take the sample (and haplotype) as arguments and were run in parallel on an HPC cluster.

## Samples

| Sample | Taxon | Data |
|---|---|---|
| Inland_ref | inland | chromosome-level reference (graph backbone) |
| Coastal_ref | coastal | chromosome-level reference |
| AH19K01, AH19K05, JF08K03, JF08K05 | coastal | PacBio HiFi, phased diploid |
| AH19K17, AH21K03, KF14K01, KF15K02 | inland | PacBio HiFi, phased diploid |
| AH14K04, MF12K02 | hybrid | PacBio HiFi, phased diploid |
| bCatMin1 | *C. minimus* (outgroup) | GCA_052095055.2 |

## 01_assembly

| Script | Step |
|---|---|
| `01_hifiasm.sh` | hifiasm assembly, hap1/hap2 pseudo-haplotypes |
| `02_genomescope.sh` | k-mer genome size and heterozygosity (KMC + GenomeScope2) |
| `03_purge_dups.sh` | Remove retained haplotypic duplication |
| `04_repurge_hap1.sh` | Re-purge three hap1 assemblies with the matching hap2 cutoffs |
| `05_busco.sh` | BUSCO (aves_odb10) on purged haplotypes |
| `06_assembly_stats.sh` | Size, contiguity and GC before and after purging |

## 02_pangenome

Chromosome-partitioned PGGB graphs against the inland reference.

| Script | Step |
|---|---|
| `01_stage_fastas.sh` | PanSN naming, bgzip and index for all 23 assemblies |
| `02_wfmash_map.sh` | Map assemblies to the inland reference |
| `03_remap_unaligned.sh` | Rescue unmapped contigs with relaxed parameters |
| `04_collect_rescues.sh` | Best rescue hit per contig |
| `05_partition.sh` | Split sequences by reference chromosome (`chrom_list.txt`) |
| `06_merge_Z.sh` | Merge the three Z scaffolds into one partition |
| `07_pggb.sh` | PGGB per chromosome |
| `08_graph_stats.sh` | odgi graph statistics |
| `09_combine_graphs.sh` | Whole-genome GFA |

The W chromosome is built but excluded from all downstream analyses.

## 03_variant_calling

Variants called from the graphs relative to the inland reference, polarized with *C. minimus*.

| Script | Step |
|---|---|
| `01_deconstruct.sh` | vg deconstruct per chromosome |
| `02_normalize.sh` | vcfbub + vcfwave |
| `03_pad_ploidy.sh` | Uniform sample columns, diploid genotypes |
| `04_merge_norm.sh` | Concatenate chromosomes, merge multiallelics, left-align |
| `05_ancestral_allele.sh` | Ancestral allele from the outgroup (`ancestral_allele.py`), outgroup removed |
| `06_genomic_regions.sh` | Exon / intron / promoter / intergenic BED from the inland annotation |
| `07_annotate_classify.sh` | Add genomic region, classify variants (`classify_variants.py`) |
| `08_biallelic.sh` | Split multiallelic sites for population genetics |
| `09_variant_table.py` | Variant counts by class, region and population |

`classify_variants.py` is adapted from the [scrub-jay pangenome](https://github.com/harvardinformatics/scrub-jay-genomics) variant classification.

## 04_repeats

Repeat annotation of all 22 thrush assemblies (`genomes.tsv`) and TE annotation of SVs. `refs/Catharus_TE_library.fasta` is the curated *Catharus* repeat library of Termignoni-Garcia et al. (2026).

| Script | Step |
|---|---|
| `01_repeatmasker.sh` | RepeatMasker with the *Catharus* library |
| `02_edta.sh` | EDTA, including structurally intact elements |
| `03_edta_summary.sh` | Masked and intact TE content, recomputed from the EDTA GFFs |
| `04_graffite.sh` | GraffiTE on SVs, annotations transferred to the biallelic VCF |
| `05_repeat_class_tables.sh` | Repeat content by class per genome and for TE-derived SVs (`class_map.awk`) |

## 05_inversions

Inversions from three approaches, merged into loci by complete-linkage clustering at 50% reciprocal overlap.

| Script | Step |
|---|---|
| `01_ragtag.sh` | Scaffold each assembly onto the inland reference |
| `02_syri_prep.sh` | Chromosome-named queries and matching references for SyRI |
| `03_syri.sh` | minimap2 + SyRI per haplotype |
| `04_svim_asm.sh` | minimap2 + SVIM-asm per bird |
| `05_jasmine.sh` | Merge SyRI and SVIM-asm calls across birds |
| `06_inversion_tables.sh` | Genotype tables per caller (`inversion_table.py`), within-caller merge (`merge_within_caller.py`), cross-caller loci (`combine_inversions.py`) |

## 06_sfs_dfe

| Script | Step |
|---|---|
| `01_annotate_features.sh` | CDS / UTR / intron / intergenic label per variant (`parse_gff.py`) |
| `02_extract_sfs.R` | Unfolded SFS per feature and variant class |
| `03_fastdfe.py` | DFE with fastDFE, ancestral misidentification estimated |

## 07_fst

`01_fst.sh`: Weir-Cockerham Fst between four inland and four coastal individuals for SNPs, INDELs and SVs, genome-wide and in 500 kb windows.

## 08_fixed_svs

| Script | Step |
|---|---|
| `01_filter_diff_svs.sh` | Strictly and nearly fixed SVs between subspecies |
| `02_sv_density.R` | Density in the chromosome 1 and 5 core blocks and full regions vs background, circular-shift permutation |
| `03_hybrid_ancestry.R` | Coastal ancestry of the hybrids at differentiated SVs, circular-shift permutation |
| `04_sv_normalization.sh` | Enrichment per SV and against windows of high SNP Fst |

## 09_assembly_artifact_check

| Script | Step |
|---|---|
| `01_window_scan.sh` | Assembly quality per 100 kb window and per SV site |
| `02_summarize.R` | Quality of the blocks vs background, enrichment after removing SVs in suspect sequence |
| `03_sniffles.sh` | Read-based SV calls with Sniffles2 |
| `04_sniffles_concordance.sh` | Joint genotyping and comparison with pangenome SVs (`concordance.py`) |
| `05_orientation.R` | Inverted segments and breakpoint detectability from RagTag alignments |

## 10_arg

Ancestral recombination graphs with SINGER v0.1.9 for the eight pure individuals, in 1 Mb tiles across the chromosome 1 (25-100 Mb) and 5 (21-60 Mb) regions and 500 random background tiles. `refs/rhomap_res_win_100kb_chrnames.csv` is the inland LDhat map of Justen et al. (2025).

| Script | Step |
|---|---|
| `01_prep_vcf.sh` | Biallelic SNPs per autosome |
| `02_window_table.sh` | SNPs and diversity per 500 kb |
| `03_recomb_maps.sh` | Recombination maps from the LDhat map |
| `04_tiles.sh` | Tile table with per-tile recombination rate |
| `05_singer.sh` | SINGER per tile, posterior samples 80-99 as tree sequences |
| `06_tile_stats.sh` | Coordinates (`tile_offset.py`), per-window statistics (`tile_stats.py`), convergence (`convergence.py`), window summary (`window_summary.py`) |
| `07_block_test.R` | Block departure from the recombination expectation (GAM, moving-block bootstrap) |
| `08_block_test_subsp.R` | As above with within-subspecies coalescence for each subspecies |

## 11_methylation

### haplotype

CpG methylation of each pseudo-haplotype in its own assembly coordinates.

| Script | Step |
|---|---|
| `01_align_cpg.sh` | pbmm2 to the bird's own assembly, pb-CpG-tools |
| `02_te_methylation.sh` | CpG counts and methylation by TE annotation and TE class |
| `03_aggregate_te_methylation.R` | Per-individual methylation, taxon and paired class tests |
| `04_te_copies.sh` | Per-copy TE methylation, full-length and fragmented |
| `05_te_copy_methylation.R` | Taxon tests on per-copy methylation |
| `06_young_te_methylation.R` | Methylation by TE age, class and family, hybrid rank tests |
| `07_hap_liftover.sh` | TE copies lifted to inland coordinates |
| `08_hap_discordance.R` | Methylation discordance between haplotypes at shared copies |
| `09_read_disorder.py` | Read-level methylation disorder per TE copy |
| `10_read_disorder_tests.R` | Hybrid rank tests on read-level disorder |

### dmr

Differentially methylated regions from reads aligned to the inland reference, with MethBat v1.1.0 (`comparisons.tsv`, `samples.tsv`).

| Script | Step |
|---|---|
| `01_align_pileup.sh` | pbmm2 to the inland reference, methbat pileup |
| `02_signature.sh` | methbat signature per comparison |
| `03_graffite_intervals.sh` | Reference intervals of TE-annotated SVs |
| `04_te_enrichment_permutation.R` | DMR-TE overlap against uniform and callable-sequence nulls |
| `05_coverage_bias.sh` | Read coverage at DMRs relative to each individual's genome-wide mean |
| `06_te_class_enrichment.R` | TE class and family overlap, position and CpG-matched nulls |
| `07_direction_fisher.R` | Hypo- vs hypermethylated DMR TE overlap, one-sided Fisher |
| `08_direction_interaction.R` | Direction x comparison type interaction, likelihood ratio test |
| `09_transgressive.R` | Hybrid DMRs classified against both parental means |
| `10_gene_features.R` | Gene features and CpG islands of DMRs, CpG-matched null |
| `11_dmr_te_age.R` | Age of the TE copies under DMRs |
