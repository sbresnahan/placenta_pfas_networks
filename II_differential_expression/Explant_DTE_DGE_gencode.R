# =============================================================================
# PFBS Differential Expression Analysis in Human Placental Explants
# =============================================================================
# 
# Purpose: Analyze differential gene and transcript expression associated with 
#          PFBS (perfluorobutane sulfonic acid) exposure in human placental 
#          explant cultures
#
# Study Design: In vitro experimental study examining dose-response effects of 
#               PFBS exposure on placental gene expression patterns at both 
#               transcript and gene levels
#
# Experimental Design:
#   - Human placental explants from multiple donors (3 third trimester, 
#     3 second trimester placentas, plus HTR-8/SVneo cell line)
#   - PFBS exposure concentrations: Control (0 μM), 5 μM, 20 μM
#   - Technical replicates for each donor-concentration combination
#   - Total samples: 24 (18 tissue + 6 cell line samples)
#
# Data: - Transcript expression quantified using Salmon against GENCODE v45
#       - Gene-level aggregation from transcript quantifications
#       - Overdispersion correction applied to account for read-to-transcript 
#         mapping ambiguity (Baldoni et al. method)
#       - Sample metadata including donor ID, trimester, treatment concentration
#
# Methods: - RUVSeq (RUVr) for unwanted variation removal
#          - DESeq2 for differential expression analysis
#          - ACAT (Aggregated Cauchy Association Test) for gene-level aggregation
#          - isotwas for isoform confirmation testing in transcript-level analysis
#          - Dose-response modeling with scaled drug concentration as continuous variable
# =============================================================================

# Load packages ---------------------------------------------------------------#

library(dplyr)
library(plyr)
library(tidyr)
library(rtracklayer)
library(ggplot2)
library(ggExtra)
library(edgeR)
library(tximport)
library(DESeq2)
library(fishpond)
library(tximeta)
library(readxl)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(BiocParallel)
library(gridExtra)
library(RUVSeq)
library(EDASeq)
library(biomaRt)
library(reshape2)
library(ACAT)
library(isotwas)

# -----------------------------------------------------------------------------#

# =============================================================================
# COMPLETED: DATA PREPROCESSING AND TRANSCRIPT ANNOTATION
# =============================================================================

# SQANTI3 isoform classification table
# Contains structural annotations for all transcripts identified by long-read sequencing
# and processed through ESPRESSO correction and SQANTI3 classification pipeline
classif <- read.table("/rsrch5/home/epi/bhattacharya_lab/data/Placenta_LRRNAseq/STB/SQANTI3/ESPRESSO_corrected_SC_filtered/ESPRESSO_corrected_SC_filtered_classification.txt",header=T)

# Reorder and rename structural categories for better interpretation
classif$structural_category <- factor(as.character(classif$structural_category),
                                      levels=c("full-splice_match",
                                               "incomplete-splice_match",
                                               "novel_in_catalog",
                                               "novel_not_in_catalog",
                                               "genic",
                                               "antisense",
                                               "fusion",
                                               "intergenic",
                                               "genic_intron"))

# Apply more intuitive naming for structural categories
classif$structural_category <- revalue(classif$structural_category, 
                                       c("full-splice_match"="FSM",
                                         "incomplete-splice_match"="ISM",
                                         "novel_in_catalog"="NIC",
                                         "novel_not_in_catalog"="NNC",
                                         "genic"="Genic Genomic",
                                         "antisense"="Antisense",
                                         "fusion"="Fusion",
                                         "intergenic"="Intergenic",
                                         "genic_intron"="Genic Intron"))
names(classif)[1] <- "transcript_id"

# Create transcript-to-gene mapping table for aggregation analyses
tx2g <- classif[,c("transcript_id","associated_gene")]
names(tx2g)[2] <- "gene_id"

# =============================================================================
# EXPERIMENTAL METADATA AND SAMPLE INFORMATION
# =============================================================================

# Load and process experimental metadata
## Sample covariates including donor information, treatment conditions, and technical details
metadata <- read.table("/rsrch5/home/epi/bhattacharya_lab/data/Elkin2025_Placenta/Drelichman_PFBS_RNA_Tracking.tsv",header=T,sep="\t")
names(metadata) <- c("sample_id","sample_type","drug_concentration","replicate","trimester")
metadata$replicate <- as.character(metadata$replicate)

# Standardize sample identifiers for consistency with quantification files
metadata$sample_id <- c("PFBS_Placenta1_Control","PFBS_Placenta1_5um","PFBS_Placenta1_20um","PFBS_Placenta2_Control",
                        "PFBS_Placenta2_5um","PFBS_Placenta2_20um","PFBS_Placenta3_Control","PFBS_Placenta3_5um",
                        "PFBS_Placenta3_20um","PFBS_2nd_placenta_N1_Control","PFBS_2nd_placenta_N1_5um","PFBS_2nd_placenta_N1_20um",
                        "PFBS_2nd_placenta_N2_Control","PFBS_2nd_placenta_N2_5um","PFBS_2nd_placenta_N2_20um","PFBS_2nd_placenta_N3_Control",
                        "PFBS_2nd_placenta_N3_5um","PFBS_2nd_placenta_N3_20um","PFBS_HTR_N1_Control","PFBS_HTR_N1_5um",
                        "PFBS_HTR_N2_Control","PFBS_HTR_N2_5um","PFBS_HTR_N3_Control","PFBS_HTR_N3_5um")

# =============================================================================
# TRANSCRIPT QUANTIFICATION AND EXPRESSION DATA LOADING
# =============================================================================

## Read in Salmon quantifications
# Salmon was used to quantify transcript expression against the long-read defined transcript annotation
dir_quant <- "/rsrch5/home/epi/bhattacharya_lab/data/Elkin2025_Placenta/salmon_quants/Placenta_LR_assembly"
dirs.salmon <- paste0(dir_quant,"/",metadata$sample_id)
files.salmon <- paste0(dirs.salmon,"/quant.sf")

# Import transcript-level quantifications with overdispersion correction
s <- catchSalmon(dirs.salmon) # transcript-level expression with overdispersion estimates

# Import gene-level aggregated expression
txi.g <- tximport(files.salmon, type = "salmon", tx2gene = tx2g,
                  geneIdCol="gene_id",txIdCol="transcript_id",dropInfReps=T,txOut=F) # gene-level expression

# Import transcript-level expression and apply overdispersion correction
## Correct tx expression for overdispersion from mapping ambiguity (Baldoni et al. 2019)
## This accounts for uncertainty in read-to-transcript assignment in transcript quantification
txi <- tximport(files.salmon, type = "salmon", tx2gene = tx2g,
                geneIdCol="gene_id",txIdCol="transcript_id",dropInfReps=T,txOut=T)

# Apply overdispersion correction to transcript counts
txi$counts <- txi$counts/s$annotation$Overdispersion
txi$abundance <- cpm(txi$counts)

# =============================================================================
# DATA ORGANIZATION AND EXPORT FOR DOWNSTREAM ANALYSES
# =============================================================================

## Save files for downstream analyses
metadata$pfas <- "PFBS"

# Ensure consistent sample naming across all data objects
dimnames(txi$abundance)[[2]] <- metadata$sample_id
dimnames(txi$counts)[[2]] <- metadata$sample_id
dimnames(txi$length)[[2]] <- metadata$sample_id

dimnames(txi.g$abundance)[[2]] <- metadata$sample_id
dimnames(txi.g$counts)[[2]] <- metadata$sample_id
dimnames(txi.g$length)[[2]] <- metadata$sample_id

# Export processed data for analysis pipeline
save(list=c("classif","metadata","tx2g","txi","txi.g"),
     file="/rsrch5/home/epi/bhattacharya_lab/data/Elkin2025_Placenta/R/Elkin2025_Placenta.RData")

# -----------------------------------------------------------------------------#


# =============================================================================
# BEGIN HERE: DIFFERENTIAL EXPRESSION ANALYSIS PIPELINE
# =============================================================================

setwd("/Users/stbresnahan/Desktop/PFAS/Elkin2025_Placenta")
a <- load("Elkin2025_Placenta_gencode.RData")

##### Loaded objects:
##### classif = SQANTI3 isoform classification table with structural annotations
##### metadata = sample covariates including treatment conditions and donor information
##### tx2g = transcript-to-gene mapping table for aggregation analyses
##### txi = transcript-level expression corrected for read-to-transcript assignment (RTA) overdispersion
##### txi.g = gene-level expression aggregated from transcript quantifications

# Add RIN scores
metadata$RIN <- as.numeric(scale(c(
  7.5, 6.6, 4.6, 5.9, 6.4, 5.7,
  8.6, 9.0, 9.3, 6.2, 8.0, 5.9,
  5.0, 5.1, 6.7, 4.0, 6.8, 7.1,
  9.9, 9.8, 9.9, 9.8, 9.9, 9.8
)))

# =============================================================================
# QUALITY CONTROL AND SAMPLE RELATIONSHIP ASSESSMENT
# =============================================================================

## Check sample relationships using Principal Component Analysis
# Create DESeq2 object for variance stabilizing transformation
dds <- DESeqDataSetFromTximport(txi = txi, colData = metadata, design=~1)
vst <- assay(varianceStabilizingTransformation(dds))

# Perform PCA to assess sample clustering and potential batch effects
pca <- prcomp(t(vst))
pca.df <- data.frame(pca$x[,1:2])
pca.df$sample_id <- metadata$sample_id
pca.df <- left_join(pca.df,metadata,by="sample_id")

# Calculate variance explained by first two principal components
var_explained <- pca$sdev^2
pve <- var_explained / sum(var_explained)
pve[1:2] * 100

# Visualize sample relationships
fig1 <- ggplot(pca.df,aes(x=PC1,y=PC2,color=sample_type)) +
  geom_point() +
  xlab("PC1: 38.79% variance") +
  ylab("PC2: 8.23% variance") +
  ggtitle("Elkin2025 RNA-seq PCA")
fig1

# =============================================================================
# SETUP PARALLEL PROCESSING AND BIOMART FOR GENE ANNOTATION
# =============================================================================

# Configure parallel processing for computational efficiency
parallel_workers <- 4
register(MulticoreParam(parallel_workers))
message(paste("Using", parallel_workers, "cores for parallel processing"))

# Setup biomaRt for gene ID conversion and annotation
message("Setting up biomaRt connection...")
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror="useast")

# =============================================================================
# GENE SYMBOL ANNOTATION FUNCTION
# =============================================================================

# Function to convert Ensembl IDs to HGNC gene symbols
get_gene_symbols <- function(ensembl_ids) {
  # Query biomaRt for gene information
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Create mapping dataframe
  df = data.frame(ensembl_gene_id = ensembl_ids)
  
  # Handle duplicate mappings by random selection (rare edge case)
  gene_info_unique = gene_info %>%
    group_by(ensembl_gene_id) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Left join and fill missing symbols with Ensembl IDs as fallback
  df = df %>%
    left_join(gene_info_unique, by = "ensembl_gene_id")
  df$hgnc_symbol = ifelse(is.na(df$hgnc_symbol),
                          df$ensembl_gene_id,
                          df$hgnc_symbol)
  
  return(df$hgnc_symbol)
}

# =============================================================================
# MAIN DIFFERENTIAL EXPRESSION ANALYSIS FUNCTION
# =============================================================================

# Comprehensive RUVr + DESeq2 pipeline for PFBS differential expression analysis
# Supports both transcript-level (with isoform aggregation/confirmation) and gene-level analyses
run_ruvr_deseq2 <- function(se, treatment_var, transcript_to_gene, 
                            covariates = c(), k = 0, 
                            just_hist = FALSE,
                            fdr_threshold = .1,
                            tx.aggregate = F) {
  # Differential expression analysis pipeline combining RUVr normalization with DESeq2
  # 
  # Arguments:
  #   se: SummarizedExperiment object containing count data and sample metadata
  #   treatment_var: character string specifying the treatment variable column name in colData(se)
  #   transcript_to_gene: data.frame with transcript_id and gene_id columns for aggregation
  #   covariates: character vector of covariate column names to include in the model (default: empty)
  #   k: integer specifying number of unwanted variation factors for RUVr correction (default: 0, no correction)
  #   just_hist: logical, if TRUE only generates p-value histogram for QC (default: FALSE)
  #   fdr_threshold: numeric, FDR threshold for significance testing in isoform confirmation (default: 0.1)
  #   tx.aggregate: logical, if TRUE performs transcript-to-gene aggregation with isoform confirmation (default: FALSE)
  # 
  # Returns:
  #   List containing:
  #     - dds: DESeqDataSet object with fitted model
  #     - results: DESeq2 results object with differential expression statistics
  #     - W: RUVr factors matrix (if k > 0)
  #     - plot: p-value histogram for quality control
  #     - results_aggregate: aggregated results with isoform confirmation (if tx.aggregate = TRUE)
  # 
  # Methods:
  #   - RUVr: Removes unwanted variation using residuals from initial model fit
  #   - DESeq2: Negative binomial modeling for differential expression
  #   - ACAT: Aggregates transcript p-values to gene level (if tx.aggregate = TRUE)
  #   - isotwas: Performs isoform confirmation testing (if tx.aggregate = TRUE)
  
  message(paste("\nRunning analysis for", treatment_var, "with k =", k))
  
  # Filter samples with complete covariate data to ensure balanced design
  vars <- c(treatment_var, covariates)
  sample_mask <- complete.cases(as.data.frame(colData(se)[, vars]))
  se <- se[, sample_mask]
  
  # Prepare count matrix and sample information for analysis
  counts_matrix <- round(assays(se)$counts)
  sample_info <- as.data.frame(colData(se))
  
  # Apply RUVr correction if k > 0 to remove unwanted variation
  if (k > 0) {
    # Create design matrix for RUVr
    design <- model.matrix(as.formula(paste("~", paste(vars, collapse = "+"))), data = sample_info)
    
    # Prepare DGEList object for RUVr
    y <- DGEList(counts = counts_matrix)
    keep <- rowSums(cpm(y) > 1) >= 3  # Filter low-expressed genes/transcripts
    y <- y[keep, ]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    
    # Fit model and extract residuals for RUVr
    fit <- glmFit(y, design)
    res <- residuals(fit, type = "deviance")
    
    # Apply RUVr to estimate unwanted variation factors
    ruvr_out <- RUVr(x = y$counts, k = k, residuals = res)
    W <- ruvr_out$W
    rownames(W) <- colData(se)$sample_id
    
    # Align W matrix with SummarizedExperiment samples
    sample_ids <- colData(se)$sample_id
    if (!all(sample_ids %in% rownames(W))) {
      stop("Mismatch between sample IDs in RUVr output and SE object")
    }
    W <- W[sample_ids, , drop = FALSE]
    colnames(W) <- paste0("W_", 1:ncol(W))
    
    # Add RUVr factors to colData for inclusion in DESeq2 model
    se_coldata <- as.data.frame(colData(se))
    rownames(se_coldata) = colData(se)$sample_id
    W_df <- as.data.frame(W)
    if (!all(rownames(se_coldata) == rownames(W_df))) {
      stop("Sample order mismatch between colData(se) and W matrix")
    }
    colData(se) <- DataFrame(cbind(se_coldata, W_df))
    vars <- c(vars, colnames(W))
  }
  
  # Run DESeq2 differential expression analysis
  design_formula <- as.formula(paste("~", paste(vars, collapse = "+")))
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix, 
                                colData = colData(se), 
                                design = design_formula)
  
  # Filter lowly expressed transcripts/genes
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  
  # Run DESeq2 differential expression pipeline
  dds <- DESeq(dds)
  res <- results(dds, name = treatment_var)
  
  # Generate p-value histogram for quality control assessment
  pvals <- na.omit(res$pvalue)
  p <- ggplot(data.frame(pvalue = pvals), aes(x = pvalue)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = paste("P-value histogram for", treatment_var, 
                       ifelse(k == 0, "(no RUVr)", paste("(RUVr k=", k, ")"))),
         x = "P-value", y = "Count")
  print(p)
  
  if (just_hist) return(invisible(NULL))
  
  # Transcript-level analysis with gene-level aggregation and isoform confirmation
  if(tx.aggregate==T){
    # ==========================================================================
    # GENE-LEVEL AGGREGATION AND ISOFORM CONFIRMATION TESTING
    # ==========================================================================
    
    cat('Aggregating isoform associations for FDR/FWER control\n')
    
    # Prepare results for aggregation
    res$transcript_id = rownames(res)
    res = merge(as.data.frame(res), transcript_to_gene, by='transcript_id')
    
    # Aggregate transcript p-values to gene level using ACAT
    # ACAT is particularly suited for combining correlated p-values from isoforms of the same gene
    res_aggregate <- res %>%
      filter(!is.na(pvalue)) %>%
      group_by(gene_id) %>%
      group_split() %>%
      lapply(function(df) {
        data.frame(
          gene_id = df$gene_id[1],
          P = ACAT::ACAT(df$pvalue)
        )
      }) %>%
      bind_rows()
    
    # Apply FDR correction at gene level for screening
    res_aggregate$FDR = p.adjust(res_aggregate$P,'fdr')
    colnames(res_aggregate)[2:3] = c('Screen_P','Screen_FDR')
    
    # Merge back with transcript results
    res = merge(res, res_aggregate, by='gene_id')
    
    # Apply isoform confirmation testing using isotwas
    # Tests whether individual isoforms within significant genes show detectable effects
    isoform_confirm <- res %>%
      group_by(gene_id) %>%
      group_modify(~ {
        gene_data <- .x
        # Only test confirmation for genes passing screening threshold
        if (gene_data$Screen_FDR[1] < fdr_threshold) {
          confirm_p <- isotwas::p_confirm(gene_data$pvalue, alpha = fdr_threshold)
          gene_data <- mutate(gene_data, Confirmation.P = confirm_p)
        } else {
          gene_data <- mutate(gene_data, Confirmation.P = 1)
        }
        gene_data
      }) %>%
      ungroup()
    
    # Return comprehensive results including both screening and confirmation
    return(list(
      dds = dds,
      results = res,
      results_aggregate = isoform_confirm,
      W = if (k > 0) W else NULL,
      plot = p
    ))
  }else{
    # Return standard results for gene-level analysis
    return(list(
      dds = dds,
      results = res,
      W = if (k > 0) W else NULL,
      plot = p
    ))
  }
}

# =============================================================================
# VOLCANO PLOT VISUALIZATION FUNCTION
# =============================================================================

# Create publication-ready volcano plots with gene symbol annotations
create_volcano_plot <- function(deseq2_results, treatment_var, output_dir = ".", 
                                fc_threshold = 0.5, fdr_threshold = 0.1, 
                                top_genes = 15, k = NULL,
                                type=c("DTE","DGE")) {
  # Generate volcano plots for differential expression results with gene symbol labeling
  # 
  # Arguments:
  #   deseq2_results: either DESeq2 results object or list from run_ruvr_deseq2() containing results
  #   treatment_var: character string specifying treatment variable name for plot title and filename
  #   output_dir: character string specifying directory path for saving plots (default: current directory)
  #   fc_threshold: numeric, absolute log2 fold change threshold for significance coloring (default: 0.5)
  #   fdr_threshold: numeric, FDR threshold for significance coloring (default: 0.1)
  #   top_genes: integer, maximum number of significant genes to label on plot (default: 15, set to 0 for no labels)
  #   k: integer, number of RUVr factors used (for plot title, auto-detected if NULL)
  #   type: character, analysis type for plot title - either "DTE" (differential transcript expression) or "DGE" (differential gene expression)
  # 
  # Returns:
  #   List containing:
  #     - plot: ggplot object of volcano plot
  #     - data: data.frame with plot data including gene symbols and significance categories
  #     - sig_genes: data.frame subset containing only significantly differentially expressed genes/transcripts
  # 
  # Color Scheme:
  #   - grey: non-significant (FDR >= threshold)
  #   - forestgreen: significant but small fold change (|log2FC| <= threshold)
  #   - maroon: significantly upregulated (log2FC > threshold, FDR < threshold)
  #   - navy: significantly downregulated (log2FC < -threshold, FDR < threshold)
  # 
  # Output:
  #   Saves PDF plot to output_dir with filename pattern: {treatment_var}_volcanoPlot_{type}.pdf
  
  # Extract results from analysis output
  if (is.list(deseq2_results) && "results" %in% names(deseq2_results)) {
    res <- deseq2_results$results
    k_value <- if(is.null(k)) {
      ifelse(is.null(deseq2_results$W), 0, ncol(deseq2_results$W))
    } else {
      k
    }
  } else {
    res <- deseq2_results
    k_value <- if(is.null(k)) 0 else k
  }
  
  # Convert to dataframe and prepare gene annotations
  res_df <- as.data.frame(res)
  
  # Extract clean Ensembl IDs (remove version numbers if present)
  if(type=="DGE"){res_df$gene_id <- row.names(res_df)}else{
    res_df$transcript_id <- row.names(res_df)
    res_df <- left_join(res_df,tx2g)
    res_df$transcript_id <- NULL
  }
  ensembl_ids <- sapply(strsplit(res_df$gene_id, '[.]'), function(x) x[1])
  res_df$ensembl <- ensembl_ids
  
  # Retrieve gene symbols for annotation
  tryCatch({
    symbols <- get_gene_symbols(ensembl_ids)
    res_df$Gene <- symbols
  }, error = function(e) {
    message("Error retrieving gene symbols: ", e$message, "\nUsing Ensembl IDs as labels.")
    res_df$Gene <- res_df$ensembl
  })
  
  # Calculate FDR and assign colors for visualization
  res_df$FDR <- p.adjust(res_df$pvalue, 'fdr')
  
  # Color coding scheme:
  # grey: non-significant
  # forestgreen: significant but small fold change
  # maroon: significantly upregulated (large fold change)
  # navy: significantly downregulated (large fold change)
  res_df$Color <- ifelse(res_df$FDR > fdr_threshold, 'grey',
                         ifelse(res_df$FDR < fdr_threshold & res_df$log2FoldChange > fc_threshold, 'maroon',
                                ifelse(res_df$FDR < fdr_threshold & res_df$log2FoldChange < -fc_threshold, 'navy',
                                       'forestgreen')))
  
  # Determine genes to label on plot
  if (top_genes > 0) {
    sig_genes <- res_df[res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold,]
    
    if (nrow(sig_genes) == 0) {
      message("No significant genes found with FDR < ", fdr_threshold, 
              " and |log2FC| > ", fc_threshold)
      res_df$Label <- ''
    } else if (nrow(sig_genes) <= top_genes) {
      res_df$Label <- ifelse(res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold, 
                             res_df$Gene, '')
    } else {
      # Select top genes by significance
      sig_genes <- sig_genes[order(sig_genes$FDR),][1:top_genes,]
      res_df$Label <- ifelse(res_df$ensembl %in% sig_genes$ensembl, res_df$Gene, '')
    }
  } else {
    res_df$Label <- ''
  }
  
  res_df$Color <- factor(res_df$Color, levels = c('forestgreen', 'grey', 'maroon', 'navy'))
  
  # Set appropriate plot title based on analysis type
  if(type=="DTE"){
    plot.title=paste0("Gene names of differentially-expressed transcripts\nRUVr k = ", k_value)
  }else{
    plot.title=paste0("Differentially expressed genes\nRUVr k = ", k_value)
  }
  
  # Create volcano plot
  volcano_plot <- ggplot(aes(x = log2FoldChange,
                             y = -log10(FDR),
                             color = Color),
                         data = res_df) +
    geom_vline(xintercept = c(fc_threshold, -fc_threshold),
               color = 'black', linetype = 2) +
    geom_hline(yintercept = -log10(fdr_threshold),
               color = 'black', linetype = 2) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c('forestgreen', 'grey', 'maroon', 'navy'), drop=F) +
    ggrepel::geom_label_repel(aes(label = Label), 
                              size = 3, max.overlaps = 100, force = 10) +
    theme_minimal() +
    xlab(expression(log[2]~"FC")) +
    ylab(expression(-log[10]~"adjusted P-value")) +
    guides(color = FALSE) +
    ggtitle(plot.title)
  
  # Save and display plot
  if(type=="DTE"){
    filename <- paste0(output_dir, "/", treatment_var, "_volcanoPlot_DTE.pdf")
  }else{
    filename <- paste0(output_dir, "/", treatment_var, "_volcanoPlot_DGE.pdf")
  }
  
  ggsave(plot = volcano_plot, filename = filename, height = 6, width = 8)
  print(volcano_plot)
  
  return(list(
    plot = volcano_plot,
    data = res_df,
    sig_genes = res_df[res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold,]
  ))
}

# =============================================================================
# CSV OUTPUT FUNCTION
# =============================================================================

# Function to export results to standardized TSV format
csv_out = function(result, treatment_var, output_dir, type=c("DTE","DGE")){
  # Export differential expression results to tab-separated value (TSV) files
  # 
  # Arguments:
  #   result: list object returned from run_ruvr_deseq2() containing DESeq2 results
  #   treatment_var: character string specifying treatment variable name for filename and metadata
  #   output_dir: character string specifying directory path for saving output files
  #   type: character, analysis type for filename - either "DTE" (differential transcript expression) or "DGE" (differential gene expression)
  # 
  # Returns:
  #   No return value (side effect: writes TSV file to disk)
  # 
  # Output Format:
  #   Tab-separated file with columns:
  #     - Ensembl: Ensembl gene/transcript ID
  #     - Gene: HGNC gene symbol (or Ensembl ID if symbol unavailable)
  #     - Transcript: transcript ID (for DTE analysis, may be NA for DGE)
  #     - Outcome: treatment variable name
  #     - log2FC: log2 fold change estimate
  #     - SE: standard error of log2 fold change
  #     - P: nominal p-value
  #     - FDR: Benjamini-Hochberg adjusted p-value
  # 
  # Output Filename:
  #   Pattern: {output_dir}/{treatment_var}_{type}.tsv
  #   Example: results/drug_concentration_numeric_DTE.tsv
  
  # Extract and annotate results
  rrr = data.frame(result$results)
  if(type=="DTE"){
    rrr$transcript_id = row.names(rrr)
    rrr <- left_join(rrr,tx2g)
    rrr$ensembl = rrr$gene_id
  }else{
    rrr$transcript_id = NA
    rrr$gene_id = row.names(rrr)
    rrr$ensembl = rrr$gene_id
  }
  
  # Add gene symbols for better interpretation
  rrr$Gene = get_gene_symbols(sapply(strsplit(rrr$ensembl,'[.]'),
                                     function(x) x[1]))
  
  # Create standardized output format
  out_rrr = data.frame(
    Ensembl = rrr$ensembl,
    Gene = rrr$Gene,
    Transcript = rrr$transcript_id,
    Outcome = treatment_var,
    log2FC = rrr$log2FoldChange,
    SE = rrr$lfcSE,
    P = rrr$pvalue,
    FDR = p.adjust(rrr$pvalue)
  )
  
  # Sort by significance for easier interpretation
  out_rrr = out_rrr[order(out_rrr$FDR, decreasing = F),]
  
  # Export to TSV with appropriate filename
  if(type=="DTE"){fname=paste0(output_dir, "/", treatment_var, "_DTE.tsv")}else{
    fname=paste0(output_dir, "/", treatment_var, "_DGE.tsv")
  }
  
  data.table::fwrite(out_rrr,
                     fname,
                     sep='\t',
                     col.names=T,
                     row.names=F,
                     quote=F)
}

# =============================================================================
# DIFFERENTIAL TRANSCRIPT EXPRESSION ANALYSES
# =============================================================================

##### Primary Analysis: Tissue samples (exclude cell line), DTE with PFBS dose-response
##### Analysis Design:
##### - Treatment variable: drug_concentration_numeric (scaled/standardized)
##### - Interpretation: DTE = significant change in transcript expression with 1 SD change in PFBS exposure
##### - Covariates: trimester (to account for developmental differences)
##### - RUVr correction: k=8 factors to remove unwanted variation

# Prepare dose-response variable
metadata$drug_concentration_numeric <- as.numeric(sub("^\\s*(\\d+).*", "\\1", metadata$drug_concentration))

## Tissue-only metadata
metadata.t <- metadata[!metadata$trimester=="cell line",]

## All samples
metadata$drug_concentration_numeric <- as.numeric(scale(metadata$drug_concentration_numeric))
metadata.t$drug_concentration_numeric <- as.numeric(scale(metadata.t$drug_concentration_numeric))

# Utility function to filter tximport objects by sample exclusion list
filter_tximport <- function(x,filter.list){
  # Filter samples from tximport objects (abundance, counts, length matrices)
  # 
  # Arguments:
  #   x: tximport object containing abundance, counts, and length matrices
  #   filter.list: character vector of sample IDs to exclude from analysis
  # 
  # Returns:
  #   tximport object with specified samples removed from all matrices
  # 
  # Usage:
  #   Commonly used to exclude specific sample types (e.g., cell lines) 
  #   or failed samples from downstream differential expression analysis
  x1 <- data.frame(x[["abundance"]])
  x[["abundance"]] <- as.matrix(x1[,!names(x1)%in%filter.list])
  x2 <- data.frame(x[["counts"]])
  x[["counts"]] <- as.matrix(x2[,!names(x2)%in%filter.list])
  x3 <- data.frame(x[["length"]])
  x[["length"]] <- as.matrix(x3[,!names(x3)%in%filter.list])
  return(x)
}

# Filter out cell line samples for tissue-specific analysis if required
txi.tissue <- filter_tximport(txi,metadata[metadata$trimester=="cell line","sample_id"])
dds.tissue <- DESeqDataSetFromTximport(txi = txi.tissue, 
                                       colData = metadata.t, 
                                       design=~drug_concentration_numeric)

dds <- DESeqDataSetFromTximport(txi = txi, 
                                colData = metadata, 
                                design=~drug_concentration_numeric)

dds.tissue <- DESeqDataSetFromTximport(txi = txi.tissue, 
                                       colData = metadata.t, 
                                       design=~drug_concentration_numeric)

# Run differential transcript expression analysis
result <- run_ruvr_deseq2(se = dds.tissue, 
                          treatment_var = "drug_concentration_numeric", 
                          covariates = c("trimester"), k = 8,
                          just_hist = F, fdr_threshold = .1,
                          transcript_to_gene = tx2g,
                          tx.aggregate = F)

# Generate volcano plot for DTE results
volcano <- create_volcano_plot(result, treatment_var = "drug_concentration_numeric", 
                               output_dir = 'results_gencode',type="DTE")

# Export DTE results to TSV
csv_out(result, treatment_var = "drug_concentration_numeric", 
        output_dir = 'results_gencode', type="DTE")


# Analysis of just term samples
## Filter out cell line samples for tissue-specific analysis if required
metadata.term <- metadata[metadata$trimester=="term",]
txi.term <- filter_tximport(txi,metadata[!metadata$trimester=="term","sample_id"])
dds.term <- DESeqDataSetFromTximport(txi = txi.term, 
                                     colData = metadata.term, 
                                     design=~drug_concentration_numeric)

vst <- assay(varianceStabilizingTransformation(dds.term))
pca <- prcomp(t(vst))
pca.df <- data.frame(pca$x[,1:2])
pca.df$sample_id <- metadata.term$sample_id
pca.df <- left_join(pca.df,metadata.term,by="sample_id")
var_explained <- pca$sdev^2
pve <- var_explained / sum(var_explained)
pve[1:2] * 100
fig2 <- ggplot(pca.df,aes(x=PC1,y=PC2,color=sample_type)) +
  geom_point() +
  xlab("PC1: 22.76% variance") +
  ylab("PC2: 16.59% variance") +
  ggtitle("Elkin2025 RNA-seq Term Explant PCA")
fig2

## Run differential transcript expression analysis
result <- run_ruvr_deseq2(se = dds.term, 
                          treatment_var = "drug_concentration_numeric", 
                          covariates = c("replicate"), k = 2,
                          just_hist = F, fdr_threshold = .1,
                          transcript_to_gene = tx2g,
                          tx.aggregate = F)

## Generate volcano plot for DTE results
volcano <- create_volcano_plot(result, treatment_var = "drug_concentration_numeric", 
                               output_dir = 'results_gencode',type="DTE")

# Export DTE results to TSV
csv_out(result, treatment_var = "drug_concentration_numeric", 
        output_dir = 'results_gencode', type="DTE")

# =============================================================================
# DIFFERENTIAL GENE EXPRESSION ANALYSES
# =============================================================================

##### Secondary Analysis: Gene-level analysis with same experimental design
##### Analysis Design:
##### - Treatment variable: drug_concentration_numeric (scaled/standardized)
##### - Interpretation: DGE = significant change in gene expression with 1 SD change in PFBS exposure
##### - Covariates: trimester (to account for developmental differences)
##### - RUVr correction: k=8 factors to remove unwanted variation
##### - No isoform aggregation needed (already at gene level)

dds.g <- DESeqDataSetFromTximport(txi = txi.g, colData = metadata, design=~1)

# Filter out cell line samples for tissue-specific gene-level analysis if required
txi.g.tissue <- filter_tximport(txi.g,metadata[metadata$trimester=="cell line","sample_id"])
dds.g.tissue <- DESeqDataSetFromTximport(txi = txi.g.tissue, 
                                         colData = metadata[!metadata$trimester=="cell line",], 
                                         design=~drug_concentration_numeric)

# Run differential gene expression analysis
result.g <- run_ruvr_deseq2(se = dds.g, 
                            treatment_var = "drug_concentration_numeric", 
                            covariates = c("trimester"), k = 10,
                            just_hist = F, fdr_threshold = .1,
                            transcript_to_gene = tx2g)

# Generate volcano plot for DGE results
volcano <- create_volcano_plot(result.g, treatment_var = "drug_concentration_numeric", 
                               output_dir = 'results_gencode', type="DGE")

# Export DGE results to TSV
csv_out(result.g, treatment_var = "drug_concentration_numeric", 
        output_dir = 'results_gencode', type="DGE")

# Analysis of just term samples
## Filter out cell line samples for tissue-specific analysis if required
txi.term.g <- filter_tximport(txi.g,metadata[!metadata$trimester=="term","sample_id"])
dds.term.g <- DESeqDataSetFromTximport(txi = txi.term.g, 
                                       colData = metadata.term, 
                                       design=~drug_concentration_numeric)

result.g.term <- run_ruvr_deseq2(se = dds.term.g, 
                                 treatment_var = "drug_concentration_numeric", 
                                 covariates = c("replicate"), k = 3,
                                 just_hist = F, fdr_threshold = .1,
                                 transcript_to_gene = tx2g)

# Generate volcano plot for DGE results
volcano <- create_volcano_plot(result.g.term, treatment_var = "drug_concentration_numeric", 
                               output_dir = 'results_gencode', type="DGE")

# Export DGE results to TSV
csv_out(result.g.term, treatment_var = "drug_concentration_numeric", 
        output_dir = 'results_gencode', type="DGE")