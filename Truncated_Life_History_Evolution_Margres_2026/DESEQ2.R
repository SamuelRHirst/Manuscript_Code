# Libraries
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(EnhancedVolcano)

# 1. Data Loading and Global Preparation ----------------------------------

# Load counts and metadata
countData_all <- as.matrix(read.csv("gene_count_matrix.csv", row.names = "gene_id"))
colData_all   <- read.csv("Cadam_Meta.csv", row.names = "SampleID")

# Factor levels to ensure consistent comparisons (Control/Reference first)
colData_all$Maturity <- factor(colData_all$Maturity, levels = c("Juvenile", "Adult"))
colData_all$Locality <- factor(colData_all$Locality, levels = c("Island", "Mainland"))

# Ensure sample order matches
countData_all <- countData_all[, rownames(colData_all)]

# Define gene sets of interest
venom_pattern <- "V-"
tf_pattern    <- "TF-"

# 2. Define Core Analysis Function ----------------------------------------

run_venom_deseq <- function(counts, metadata, design_formula, contrast_vec, 
                            output_name, fc_cutoff = 1, p_cutoff = 0.05, 
                            x_lim = c(-5, 5), y_lim = c(0, 5)) {
  
  # Initialize DESeq2
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = design_formula)
  dds <- DESeq(dds)
  res <- results(dds, contrast = contrast_vec)
  
  # Prepare results dataframe
  all_genes <- data.frame(
    gene = rownames(res),
    FDR = res$padj,
    log2FoldChange = res$log2FoldChange
  ) %>% filter(!is.na(FDR))
  
  # Identify significant Venom and TF genes for highlighting
  sig_venom <- all_genes %>% 
    filter(grepl(venom_pattern, gene) & abs(log2FoldChange) >= fc_cutoff & FDR < p_cutoff)
  
  sig_tf <- all_genes %>% 
    filter(grepl(tf_pattern, gene) & abs(log2FoldChange) >= fc_cutoff & FDR < p_cutoff)
  
  # Save results
  write.csv(rbind(sig_venom, sig_tf), paste0("DESEQ2_", output_name, ".csv"))
  
  # Define custom colors for Volcano plot
  keyvals <- ifelse(res$log2FoldChange < -fc_cutoff, '#AA000E', # Comparison Group 2 Biased
             ifelse(res$log2FoldChange > fc_cutoff, '#00008B',  # Comparison Group 1 Biased
             'black'))
  
  # Override colors for specific categories
  keyvals[rownames(res) %in% sig_venom$gene] <- '#CCFF23' # Venom
  keyvals[rownames(res) %in% sig_tf$gene]    <- '#5D148F' # TF
  keyvals[is.na(keyvals)] <- 'black'
  
  names(keyvals)[keyvals == '#00008B'] <- paste(contrast_vec[2], "Biased")
  names(keyvals)[keyvals == '#AA000E'] <- paste(contrast_vec[3], "Biased")
  names(keyvals)[keyvals == '#CCFF23'] <- 'Biased Venom Gene'
  names(keyvals)[keyvals == '#5D148F'] <- 'TF'

  # Generate Plot
  plot <- EnhancedVolcano(res,
                  lab = rownames(res),
                  x = 'log2FoldChange', y = 'padj',
                  selectLab = sig_venom$gene,
                  pCutoff = p_cutoff, FCcutoff = fc_cutoff,
                  pointSize = 2.0, labSize = 3,
                  drawConnectors = TRUE,
                  title = paste("Comparison:", output_name),
                  subtitle = paste("FC Cutoff:", fc_cutoff, "| p-value:", p_cutoff),
                  xlim = x_lim, ylim = y_lim,
                  colCustom = keyvals)
  
  return(plot)
}

# 3. Execute Comparisons -------------------------------------------------

# Example 1: Adult Mainland vs Island
adult_meta <- colData_all %>% filter(Maturity == "Adult")
adult_counts <- countData_all[, rownames(adult_meta)]

run_venom_deseq(adult_counts, adult_meta, ~Locality, 
                c("Locality", "Mainland", "Island"), 
                "Adult_Mainland_vs_Island", fc_cutoff = 1)

# Example 2: Island Adult vs Juvenile
island_meta <- colData_all %>% filter(Locality == "Island")
island_counts <- countData_all[, rownames(island_meta)]

run_venom_deseq(island_counts, island_meta, ~Maturity, 
                c("Maturity", "Adult", "Juvenile"), 
                "Island_Adult_vs_Juvenile", fc_cutoff = 1)

# Example 3: Mainland Adult vs Juvenile
mainland_meta <- colData_all %>% filter(Locality == "Mainland")
mainland_counts <- countData_all[, rownames(mainland_meta)]

run_venom_deseq(mainland_counts, mainland_meta, ~Maturity, 
                c("Maturity", "Adult", "Juvenile"), 
                "Mainland_Adult_vs_Juvenile", fc_cutoff = 2, y_lim = c(0, 45))
