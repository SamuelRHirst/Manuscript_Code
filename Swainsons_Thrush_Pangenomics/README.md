# Swainson's Thrush Pangenomics

Code for the Swainson's thrush (*Catharus ustulatus*) structural variant pangenome.

Scripts are numbered in the order they were run. Each is run from the project root, which is expected to contain:

```
reads/       <sample>.hifi_reads.fastq.gz
refs/        inland and coastal reference assemblies, C. minimus assembly
assembly/    created by 01_assembly
pangenome/   created by 02_pangenome
```

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
