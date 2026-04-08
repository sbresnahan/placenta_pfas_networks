require(readxl)
require(SummarizedExperiment)
require(limma)
require(edgeR)
require(DESeq2)

setwd("/Users/stbresnahan/Desktop/PFAS")

# Load transcript and gene expression estimates
## Quantified using Salmon against our long-read defined placenta-specific transcript annotation
## Corrected for overdispersion from read-to-transcript mapping ambiguity following methods by Baldoni et al.

a = load('GUSTO_GENCODE_quants.RData')
b = load('PFAS_in.RData')

cols.all <- colnames(txi.gencode.g$abundance)
cols.keep <- gse_gusto.tx@colData@rownames
cols.filter <- setdiff(cols.all,cols.keep)

filter_tximport <- function(x,filter.list){
  x1 <- data.frame(x[["abundance"]])
  x[["abundance"]] <- as.matrix(x1[,!names(x1)%in%filter.list])
  x2 <- data.frame(x[["counts"]])
  x[["counts"]] <- as.matrix(x2[,!names(x2)%in%filter.list])
  x3 <- data.frame(x[["length"]])
  x[["length"]] <- as.matrix(x3[,!names(x3)%in%filter.list])
  return(x)
}

txi.gencode.g <- filter_tximport(txi.gencode.g,cols.filter)
txi.gencode <- filter_tximport(txi.gencode,cols.filter)

gse_gusto.tx <- SummarizedExperiment(
  assays = list(
    counts    = txi.gencode.g$counts,
    abundance = txi.gencode.g$abundance,
    length    = txi.gencode.g$length
  ),
  colData = colData(gse_gusto.tx)
)

se_gusto.tx <- SummarizedExperiment(
  assays = list(
    counts    = txi.gencode$counts,
    abundance = txi.gencode$abundance,
    length    = txi.gencode$length
  ),
  colData = colData(se_gusto.tx)
)

# Load covariates matrices
covs = read_excel('20220823-Full_200_RNAseq_covars_v2.xlsx')

pfas_mat = read_excel('mat_pfas_imputed.xlsx')
colnames(pfas_mat)[-1] = paste0(colnames(pfas_mat)[-1],'_mat')
pfas_child = read_excel('cord_pfas_imputed.xlsx')
colnames(pfas_child)[-1] = paste0(colnames(pfas_child)[-1],'_child')

# Load transcript-to-gene mapping tables
b = load("tx2g_ESPRESSO_assembly_SC_filtered.RData")


# First, merge the PFAS maternal and child data
# Convert both to data frames (in case they're tibbles)
pfas_mat_df <- as.data.frame(pfas_mat)
pfas_child_df <- as.data.frame(pfas_child)

# Check the name of the identifier column (assuming it's the first column)
id_column_name <- colnames(pfas_mat_df)[1]
# Confirm it's the same in both datasets
stopifnot(id_column_name == colnames(pfas_child_df)[1])

# Merge maternal and child PFAS data
pfas_combined <- merge(pfas_mat_df, pfas_child_df, by = id_column_name, all = TRUE)

# Now merge with colData from SummarizedExperiment objects
# For se_gusto.tx (transcript-level)
se_coldata <- colData(se_gusto.tx)
se_coldata_df <- as.data.frame(se_coldata)
covs = subset(covs,SubjectID %in% se_coldata_df$SubjectID)
covs = covs[match(se_coldata_df$SubjectID,covs$SubjectID),]
se_coldata_df$sex = covs$sex

# Merge PFAS data with transcript-level colData
merged_coldata_tx <- merge(se_coldata_df, pfas_combined, 
                           by.x = "SubjectID", by.y = 'SubjectID', 
                           all.x = TRUE)

# Update the colData in se_gusto.tx
colData(se_gusto.tx) <- DataFrame(merged_coldata_tx)

# For gse_gusto.tx (gene-level)
gse_coldata <- colData(gse_gusto.tx)
gse_coldata_df <- as.data.frame(gse_coldata)

# Merge PFAS data with gene-level colData
merged_coldata_gene <- merge(gse_coldata_df, pfas_combined, 
                             by.x = "SubjectID", by.y = 'SubjectID', 
                             all.x = TRUE)

# Update the colData in gse_gusto.tx
colData(gse_gusto.tx) <- DataFrame(merged_coldata_gene)



# Set up for differential expression analysis with RUVSeq and DESeq2
library(ggplot2)
library(DESeq2)
library(EnhancedVolcano)
library(BiocParallel)
library(gridExtra)
library(RUVSeq)
library(EDASeq)
library(edgeR)  # Make sure edgeR is loaded
library(biomaRt) # For gene ID conversion
library(reshape2) # For melting data frames
library(ACAT)
library(isotwas)

# Set up parallel processing with 4 cores
parallel_workers <- 4
register(MulticoreParam(parallel_workers))
message(paste("Using", parallel_workers, "cores for parallel processing"))

# Set up biomaRt for gene ID conversion
message("Setting up biomaRt connection...")
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",mirror="useast")


# Function to get gene symbols from Ensembl IDs
get_gene_symbols <- function(ensembl_ids) {
  gene_info <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  
  # Create a named vector for easy mapping
  df = data.frame(ensembl_gene_id = ensembl_ids)
  
  # Deduplicate gene_info by randomly selecting one row per ensembl_gene_id
  gene_info_unique = gene_info %>%
    group_by(ensembl_gene_id) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Left join with df
  df = df %>%
    left_join(gene_info_unique, by = "ensembl_gene_id")
  df$hgnc_symbol = ifelse(is.na(df$hgnc_symbol),
                          df$ensembl_gene_id,
                          df$hgnc_symbol)
  
  return(df$hgnc_symbol)
}
# Simplified and safer RUVr + DESeq2 pipeline for PFAS analysis

library(DESeq2)
library(RUVSeq)
library(edgeR)
library(ggplot2)
library(SummarizedExperiment)


run_ruvr_deseq2 <- function(se, pfas_var, transcript_to_gene, 
                            covariates = c(), k = 0, 
                            output_dir = ".", just_hist = FALSE,
                            fdr_threshold = .1,aggregate=F) {
  message(paste("\nRunning analysis for", pfas_var, "with k =", k))
  
  # Subset samples with complete covariates
  vars <- c(pfas_var, covariates)
  sample_mask <- complete.cases(as.data.frame(colData(se)[, vars]))
  se <- se[, sample_mask]
  
  # Round counts and prepare sample info
  counts_matrix <- round(assays(se)$counts)
  sample_info <- as.data.frame(colData(se))
  
  if (k > 0) {
    design <- model.matrix(as.formula(paste("~", paste(vars, collapse = "+"))), data = sample_info)
    y <- DGEList(counts = counts_matrix)
    keep <- rowSums(cpm(y) > 1) >= 3
    y <- y[keep, ]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    res <- residuals(fit, type = "deviance")
    
    ruvr_out <- RUVr(x = y$counts, k = k, residuals = res)
    W <- ruvr_out$W
    rownames(W) <- colData(se)$SubjectID  # Set sample names explicitly
    
    # Subset W to match se sample order
    sample_ids <- colData(se)$SubjectID
    if (!all(sample_ids %in% rownames(W))) {
      stop("Mismatch between sample IDs in RUVr output and SE object")
    }
    W <- W[sample_ids, , drop = FALSE]
    colnames(W) <- paste0("W_", 1:ncol(W))
    
    # Add W to colData safely
    se_coldata <- as.data.frame(colData(se))
    rownames(se_coldata) = colData(se)$SubjectID
    W_df <- as.data.frame(W)
    if (!all(rownames(se_coldata) == rownames(W_df))) {
      stop("Sample order mismatch between colData(se) and W matrix")
    }
    colData(se) <- DataFrame(cbind(se_coldata, W_df))
    vars <- c(vars, colnames(W))
  }
  
  # DESeq2 analysis
  design_formula <- as.formula(paste("~", paste(vars, collapse = "+")))
  dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = colData(se), design = design_formula)
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds)
  res <- results(dds, name = pfas_var)
  
  # Histogram of p-values
  pvals <- na.omit(res$pvalue)
  p <- ggplot(data.frame(pvalue = pvals), aes(x = pvalue)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black") +
    theme_minimal() +
    labs(title = paste("P-value histogram for", pfas_var, ifelse(k == 0, "(no RUVr)", paste("(RUVr k=", k, ")"))),
         x = "P-value", y = "Count")
  print(p)
  
  if (just_hist) return(invisible(NULL))
  
  if(aggregate){
    cat('Aggregating isoform associations for FDR/FWER control')
    require(tidyverse)
    require(ACAT)
    res$transcript_id = rownames(res)
    res = merge(as.data.frame(res),
                transcript_to_gene,
                by='transcript_id')
    res_aggregate = res %>%
      group_by(gene_id) %>%
      summarize(P = ACAT(pvalue))
    res_aggregate$FDR = p.adjust(res_aggregate$P,'fdr')
    colnames(res_aggregate)[2:3] = c('Screen_P','Screen_FDR')
    
    require(isotwas)
    res = merge(res,res_aggregate,by='gene_id')
    isoform_confirm <- res %>%
      group_by(gene_id) %>%
      group_modify(~ {
        gene_data <- .x
        if (gene_data$Screen_FDR[1] < fdr_threshold) {
          confirm_p <- isotwas::p_confirm(gene_data$pvalue, alpha = fdr_threshold)
          gene_data <- mutate(gene_data, Confirmation.P = confirm_p)
        } else {
          gene_data <- mutate(gene_data, Confirmation.P = 1)
        }
        gene_data
      }) %>%
      ungroup()
  }else{
    isoform_confirm="NA"
  }
  
  # Return results
  return(list(
    dds = dds,
    results = res,
    results_aggregate = isoform_confirm,
    W = if (k > 0) W else NULL,
    plot = p
  ))
}

# Function to create a volcano plot from DESeq2 results
create_volcano_plot <- function(deseq2_results, pfas_var, output_dir = ".", 
                                fc_threshold = 0.5, fdr_threshold = 0.1, 
                                top_genes = 15, k = NULL) {
  # Extract results from run_ruvr_deseq2 output
  if (is.list(deseq2_results) && "results_gencode" %in% names(deseq2_results)) {
    # If input is the complete output from run_ruvr_deseq2
    res <- deseq2_results$results
    k_value <- if(is.null(k)) {
      ifelse(is.null(deseq2_results$W), 0, ncol(deseq2_results$W))
    } else {
      k
    }
  } else {
    # If input is just the DESeq2 results object
    res <- deseq2_results
    k_value <- if(is.null(k)) 0 else k
  }
  
  # Convert results to a data frame
  res_df <- as.data.frame(res)
  
  # Extract Ensembl IDs (remove version numbers) from rownames
  ensembl_ids <- sapply(strsplit(res_df$gene_id, '[.]'), function(x) x[1])
  res_df$ensembl <- ensembl_ids
  
  # Get gene symbols using biomaRt - providing clean Ensembl IDs
  tryCatch({
    # First try using existing get_gene_symbols function
    symbols <- get_gene_symbols(ensembl_ids)
    res_df$Gene <- symbols
  }, error = function(e) {
    message("Error retrieving gene symbols: ", e$message, "\nUsing Ensembl IDs as labels.")
    res_df$Gene <- res_df$ensembl
  })
  
  # Calculate FDR
  res_df$FDR <- p.adjust(res_df$pvalue, 'fdr')
  
  # Add color coding for the plot
  res_df$Color <- ifelse(res_df$FDR > fdr_threshold, 'grey',
                         ifelse(res_df$FDR < fdr_threshold & res_df$log2FoldChange > fc_threshold, 'maroon',
                                ifelse(res_df$FDR < fdr_threshold & res_df$log2FoldChange < -fc_threshold, 'navy',
                                       'forestgreen')))
  
  # Determine which genes to label
  if (top_genes > 0) {
    # Find significant genes (based on FDR and fold change)
    sig_genes <- res_df[res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold,]
    
    # If fewer than top_genes are significant, just use all of them
    if (nrow(sig_genes) == 0) {
      message("No significant genes found with FDR < ", fdr_threshold, 
              " and |log2FC| > ", fc_threshold)
      res_df$Label <- ''
    } else if (nrow(sig_genes) <= top_genes) {
      res_df$Label <- ifelse(res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold, 
                             res_df$Gene, '')
    } else {
      # Otherwise, sort by significance and take the top ones
      sig_genes <- sig_genes[order(sig_genes$FDR),][1:top_genes,]
      res_df$Label <- ifelse(res_df$ensembl %in% sig_genes$ensembl, res_df$Gene, '')
    }
  } else {
    # If top_genes is 0 or negative, don't label any genes
    res_df$Label <- ''
  }
  
  res_df$Color <- factor(res_df$Color, levels = c('forestgreen', 'grey', 'maroon', 'navy'))
  
  # Create volcano plot
  volcano_plot <- ggplot(aes(x = log2FoldChange,
                             y = -log10(FDR),
                             color = Color),
                         data = res_df) +
    geom_vline(xintercept = c(fc_threshold, -fc_threshold),
               color = 'black',
               linetype = 2) +
    geom_hline(yintercept = -log10(fdr_threshold),
               color = 'black',
               linetype = 2) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c('forestgreen', 'grey', 'maroon', 'navy'),drop=F) +
    ggrepel::geom_label_repel(aes(label = Label), 
                              size = 3,
                              max.overlaps = 100,
                              force = 10) +
    theme_minimal() +
    xlab(expression(log[2]~"FC")) +
    ylab(expression(-log[10]~"adjusted P-value")) +
    guides(color = FALSE) +
    ggtitle(paste0("Gene names of differentially-expressed transcripts (",
                   ifelse(grepl('_mat', pfas_var), 'Maternal', 'Fetal'),
                   ' ', strsplit(pfas_var, '_')[[1]][1], ')\nRUVr k = ', k_value))
  
  # Save the plot
  filename <- paste0(output_dir, "/", pfas_var, "_volcanoPlot_DTE.pdf")
  ggsave(plot = volcano_plot, filename = filename, height = 6, width = 8)
  
  # Also display the plot
  print(volcano_plot)
  
  # Return the plot and processed data
  return(list(
    plot = volcano_plot,
    data = res_df,
    sig_genes = res_df[res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold,]
  ))
}

# Function to create a heatmap of significant genes using ggplot2
create_heatmap <- function(deseq2_results, pfas_var, output_dir = ".", 
                           fdr_threshold = 0.1, fc_threshold = 0.5, 
                           max_genes = 50, k = NULL) {
  # Load required packages
  library(ggplot2)
  library(scales)
  library(reshape2)
  
  # Extract results from run_ruvr_deseq2 output
  if (is.list(deseq2_results) && "results_gencode" %in% names(deseq2_results)) {
    # If input is the complete output from run_ruvr_deseq2
    res <- deseq2_results$results
    dds <- deseq2_results$dds
    k_value <- if(is.null(k)) {
      ifelse(is.null(deseq2_results$W), 0, ncol(deseq2_results$W))
    } else {
      k
    }
  } else {
    # If just the DESeq2 results object was provided, we need the dds object
    stop("Need both results and dds objects. Please provide the full output from run_ruvr_deseq2().")
  }
  
  # Convert results to a data frame
  res_df <- as.data.frame(res)
  
  # Extract clean Ensembl IDs (remove version numbers) from rownames
  ensembl_ids <- sapply(strsplit(rownames(res_df), '[.]'), function(x) x[1])
  res_df$ensembl <- ensembl_ids
  
  # Get gene symbols using biomaRt - providing clean Ensembl IDs
  tryCatch({
    # First try using existing get_gene_symbols function
    symbols <- get_gene_symbols(ensembl_ids)
    res_df$Gene <- symbols
  }, error = function(e) {
    message("Error retrieving gene symbols: ", e$message, "\nUsing Ensembl IDs as labels.")
    res_df$Gene <- res_df$ensembl
  })
  
  # Calculate FDR
  res_df$FDR <- p.adjust(res_df$pvalue, 'fdr')
  
  # Filter for significant genes
  sig_genes <- res_df[res_df$FDR < fdr_threshold & abs(res_df$log2FoldChange) > fc_threshold,]
  
  # Check if we have any significant genes
  if (nrow(sig_genes) == 0) {
    message("No significant genes found for ", pfas_var, " using FDR < ", fdr_threshold, 
            " and |log2FC| > ", fc_threshold)
    return(NULL)
  }
  
  # If too many significant genes, select the top ones by significance
  if (nrow(sig_genes) > max_genes) {
    message("Found ", nrow(sig_genes), " significant genes. Limiting heatmap to top ", max_genes, 
            " genes by significance.")
    sig_genes <- sig_genes[order(sig_genes$FDR),][1:max_genes,]
  }
  
  message("Creating heatmap with ", nrow(sig_genes), " significant genes.")
  
  # Get normalized counts for significant genes
  # Use variance stabilizing transformation for better visualization
  vst_dds <- vst(dds)
  
  # Prepare the count matrix with proper IDs for matching
  vst_counts <- assay(vst_dds)
  
  # Extract the base Ensembl IDs (without version) from rownames for matching
  vst_ensembl_ids <- sapply(strsplit(rownames(vst_counts), '[.]'), function(x) x[1])
  
  # Create a matching table to map between rownames and clean Ensembl IDs
  matching_table <- data.frame(
    original_id = rownames(vst_counts),
    ensembl_id = vst_ensembl_ids,
    stringsAsFactors = FALSE
  )
  
  # Find the significant genes in the vst count matrix
  sig_gene_indices <- match(sig_genes$ensembl, matching_table$ensembl_id)
  
  # Check if all significant genes were found
  if (all(is.na(sig_gene_indices))) {
    message("None of the significant genes found in the count matrix. Check gene ID matching.")
    return(NULL)
  }
  
  # Filter out NAs (genes not found in the count matrix)
  sig_gene_indices <- sig_gene_indices[!is.na(sig_gene_indices)]
  sig_genes_found <- sig_genes[sig_genes$ensembl %in% matching_table$ensembl_id[sig_gene_indices],]
  
  if (nrow(sig_genes_found) < nrow(sig_genes)) {
    message("Only ", nrow(sig_genes_found), " out of ", nrow(sig_genes), 
            " significant genes found in count matrix.")
  }
  
  # Get the original rownames for these genes
  original_ids <- matching_table$original_id[sig_gene_indices]
  
  # Extract counts for these genes
  sig_counts <- vst_counts[original_ids,]
  
  # Create gene display names for the plot
  gene_display_names <- rep("", nrow(sig_counts))
  ensembl_ids_clean <- sapply(strsplit(rownames(sig_counts), '[.]'), function(x) x[1])
  
  # Try to use gene symbols if available
  if ("Gene" %in% colnames(sig_genes_found)) {
    # Create mapping from Ensembl ID to gene symbol
    gene_symbol_map <- setNames(sig_genes_found$Gene, sig_genes_found$ensembl)
    
    # Set display names
    for (i in 1:length(ensembl_ids_clean)) {
      if (ensembl_ids_clean[i] %in% names(gene_symbol_map)) {
        gene_display_names[i] <- gene_symbol_map[ensembl_ids_clean[i]]
      } else {
        gene_display_names[i] <- ensembl_ids_clean[i]
      }
    }
  } else {
    gene_display_names <- ensembl_ids_clean
  }
  
  # Create a data frame for sample annotations
  sample_anno <- data.frame(Sample = colnames(sig_counts), stringsAsFactors = FALSE)
  rownames(sample_anno) <- sample_anno$Sample
  
  # Add sex annotation if available
  if ("sex" %in% colnames(colData(dds))) {
    sample_anno$Sex <- colData(dds)$sex
  }
  
  # Add PFAS variable value
  if (pfas_var %in% colnames(colData(dds))) {
    sample_anno$PFAS <- colData(dds)[[pfas_var]]
  }
  
  # Add GA annotation if available
  if ("GA" %in% colnames(colData(dds))) {
    sample_anno$GA <- colData(dds)$GA
  }
  
  # Scale the data by row (z-score)
  sig_counts_scaled <- t(scale(t(sig_counts)))
  
  # Reshape data for ggplot2
  melted_data <- reshape2::melt(sig_counts_scaled, varnames = c("Gene", "Sample"), value.name = "Expression")
  
  # Replace gene identifiers with display names
  melted_data$Gene <- factor(melted_data$Gene, 
                             levels = rownames(sig_counts), 
                             labels = gene_display_names)
  
  # Add sample annotations
  melted_data$Sex <- sample_anno$Sex[match(melted_data$Sample, rownames(sample_anno))]
  if ("PFAS" %in% colnames(sample_anno)) {
    melted_data$PFAS <- sample_anno$PFAS[match(melted_data$Sample, rownames(sample_anno))]
  }
  if ("GA" %in% colnames(sample_anno)) {
    melted_data$GA <- sample_anno$GA[match(melted_data$Sample, rownames(sample_anno))]
  }
  
  # Determine sample clustering (if desired)
  # This uses hierarchical clustering to order the samples
  sample_hc <- hclust(dist(t(sig_counts_scaled)), method = "complete")
  sample_order <- sample_hc$labels[sample_hc$order]
  
  # Determine gene clustering
  gene_hc <- hclust(dist(sig_counts_scaled), method = "complete")
  gene_order <- gene_hc$labels[gene_hc$order]
  
  # Apply the clustering order to factors
  melted_data$Sample <- factor(melted_data$Sample, levels = sample_order)
  melted_data$Gene <- factor(melted_data$Gene, 
                             levels = gene_display_names[match(gene_order, rownames(sig_counts))])
  
  # Create a title for the plot
  plot_title <- paste("Significant genes for", pfas_var, "(RUVr k =", k_value, ")")
  
  # Create the heatmap using ggplot2
  heatmap_plot <- ggplot(melted_data, aes(x = Sample, y = Gene, fill = Expression)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "navy",
      mid = "white",
      high = "firebrick3",
      midpoint = 0,
      name = "Z-score"
    ) +
    theme_minimal() +
    labs(
      title = plot_title,
      x = "Sample",
      y = "Gene"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank()
    )
  
  # If we have sample annotations, add them as colored bars above the heatmap
  annotation_plots <- list()
  
  # Sex annotation
  if ("Sex" %in% colnames(melted_data)) {
    # Get unique samples in the correct order
    samples_ordered <- levels(melted_data$Sample)
    sex_data <- data.frame(
      Sample = samples_ordered,
      Sex = sample_anno$Sex[match(samples_ordered, rownames(sample_anno))],
      stringsAsFactors = FALSE
    )
    
    sex_plot <- ggplot(sex_data, aes(x = Sample, y = 1, fill = Sex)) +
      geom_tile() +
      scale_fill_manual(values = c("pink", "lightblue")) +
      theme_void() +
      theme(
        legend.position = "top",
        axis.text.x = element_blank()
      ) +
      labs(title = "Sex")
    
    annotation_plots$sex <- sex_plot
  }
  
  # PFAS annotation
  if ("PFAS" %in% colnames(melted_data)) {
    samples_ordered <- levels(melted_data$Sample)
    pfas_data <- data.frame(
      Sample = samples_ordered,
      PFAS = sample_anno$PFAS[match(samples_ordered, rownames(sample_anno))],
      stringsAsFactors = FALSE
    )
    
    pfas_plot <- ggplot(pfas_data, aes(x = Sample, y = 1, fill = PFAS)) +
      geom_tile() +
      scale_fill_gradient(low = "lightblue", high = "darkblue") +
      theme_void() +
      theme(
        legend.position = "top",
        axis.text.x = element_blank()
      ) +
      labs(title = pfas_var)
    
    annotation_plots$pfas <- pfas_plot
  }
  
  # GA annotation
  if ("GA" %in% colnames(melted_data)) {
    samples_ordered <- levels(melted_data$Sample)
    ga_data <- data.frame(
      Sample = samples_ordered,
      GA = sample_anno$GA[match(samples_ordered, rownames(sample_anno))],
      stringsAsFactors = FALSE
    )
    
    ga_plot <- ggplot(ga_data, aes(x = Sample, y = 1, fill = GA)) +
      geom_tile() +
      scale_fill_gradient(low = "lightyellow", high = "orange") +
      theme_void() +
      theme(
        legend.position = "top",
        axis.text.x = element_blank()
      ) +
      labs(title = "GA")
    
    annotation_plots$ga <- ga_plot
  }
  
  # Combine all plots using cowplot or patchwork if annotations exist
  final_plot <- heatmap_plot
  
  # Use patchwork or cowplot to combine plots
  if (length(annotation_plots) > 0) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      # Use patchwork if available
      library(patchwork)
      
      # Start with an empty plot list
      combined_plot <- NULL
      
      # Add each annotation plot
      for (anno_name in names(annotation_plots)) {
        if (is.null(combined_plot)) {
          combined_plot <- annotation_plots[[anno_name]]
        } else {
          combined_plot <- combined_plot / annotation_plots[[anno_name]]
        }
      }
      
      # Add the heatmap at the bottom with more space
      final_plot <- combined_plot / heatmap_plot + 
        plot_layout(heights = c(rep(0.1, length(annotation_plots)), 1))
      
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      # Use cowplot if available
      library(cowplot)
      
      # Calculate annotation plot heights
      anno_height <- 0.1
      plot_heights <- c(rep(anno_height, length(annotation_plots)), 1)
      
      # Combine all plots
      plots_to_combine <- c(annotation_plots, list(heatmap = heatmap_plot))
      final_plot <- cowplot::plot_grid(
        plotlist = plots_to_combine,
        ncol = 1,
        align = 'v',
        rel_heights = plot_heights
      )
    } else {
      message("For annotation bars, please install either 'patchwork' or 'cowplot' package")
      final_plot <- heatmap_plot
    }
  }
  
  # Save the plot
  ggsave(
    filename = paste0(output_dir, "/", pfas_var, "_significant_genes_heatmap_ggplot_DTE.pdf"),
    plot = final_plot,
    width = 10,
    height = min(12, 5 + nrow(sig_counts_scaled) * 0.15),
    limitsize = FALSE
  )
  
  # Also display the plot
  print(final_plot)
  
  # Return the data used for the heatmap
  return(list(
    plot = final_plot,
    data = melted_data,
    sig_counts = sig_counts,
    sig_genes = sig_genes_found,
    sample_anno = sample_anno
  ))
}

csv_out = function(result,pfas_var,output_dir,tx2g){
  
  rrr = data.frame(result$results)
  rrr$transcript_id <- row.names(rrr)
  rrr <- left_join(rrr,tx2g.gencode,by="transcript_id")
  rrr$ensembl = rrr$gene_id
  
  rrr$Gene = get_gene_symbols(sapply(strsplit(rrr$ensembl,'[.]'),
                                     function(x) x[1]))
  rrr = as.data.frame(rrr)
  out_rrr = data.frame(Ensembl = rrr$ensembl,
                       Gene = rrr$Gene,
                       Transcript = rrr$transcript_id,
                       Outcome = pfas_var,
                       log2FC = rrr$log2FoldChange,
                       SE = rrr$lfcSE,
                       P = rrr$pvalue,
                       FDR = p.adjust(rrr$pvalue))
  out_rrr = out_rrr[order(out_rrr$FDR,decreasing = F),]
  # out_rrr = out_rrr[!duplicated(out_rrr$Ensembl),]
  
  data.table::fwrite(out_rrr,
                     paste0(output_dir, "/", pfas_var, "_DTE.tsv"),
                     sep='\t',
                     col.names=T,
                     row.names=F,
                     quote=F)
  
}

# Visualize correlation in effect sizes 
vis_efcor <- function(result_mat,result_child,exposure,output_dir){
  
  df_mat <- as.data.frame(result_mat$results)
  df_child <- as.data.frame(result_child$results)

  df_merged <- dplyr::inner_join(
    df_mat %>% dplyr::select(transcript_id, gene_id, log2FoldChange, pvalue) %>% rename(lfc_maternal = log2FoldChange,
                                                                  pvalue_maternal = pvalue),
    df_child %>% dplyr::select(transcript_id, log2FoldChange, pvalue) %>% rename(lfc_fetal = log2FoldChange,
                                                                    pvalue_fetal = pvalue),
    by = "transcript_id"
  )
  
  # Get concordant transcripts
  concordant <- df_merged %>%
    dplyr::filter(
      pvalue_maternal < 0.05,
      pvalue_fetal < 0.05,
      abs(lfc_maternal) > 0.5,
      abs(lfc_fetal) > 0.5,
      sign(lfc_maternal) == sign(lfc_fetal)
    ) %>%
    dplyr::pull(transcript_id)
  
  if (length(concordant) == 0) {
    df_merged$symbol <- ""
  }else{
    concordant <- data.frame(transcript_id=concordant)
    concordant <- left_join(concordant,df_merged[,c("transcript_id","gene_id")],by="transcript_id")
    concordant$ensembl <- sapply(strsplit(concordant$gene_id, '[.]'), function(x) x[1])
    symbols <- get_gene_symbols(concordant$ensembl)
    concordant$symbol <- symbols
    concordant <- concordant[,c(2,4)]
    concordant <- concordant[!duplicated(concordant),]
    df_merged <- left_join(df_merged,concordant)
    df_merged[is.na(df_merged$symbol),"symbol"] <- ""
  }
  
  df_merged <- df_merged %>%
    dplyr::mutate(symbol = ifelse(
      abs(lfc_maternal) < 0.5 |
        abs(lfc_fetal) < 0.5 |
        pvalue_maternal > 0.05 |
        pvalue_fetal > 0.05,
      "", symbol
    ))
  
  # Plot correlation of effect sizes
  
  exposure <- strsplit(exposure,"_")[[1]][1]
  
  cor_val <- cor.test(df_merged$lfc_maternal, df_merged$lfc_fetal, method = "pearson")
  
  g1 <- ggplot(df_merged, aes(x = lfc_maternal, y = lfc_fetal)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      title = paste0(paste(paste0("Correlation of DTE Effect Sizes between","\nMaternal and Fetal"),exposure,sep=" "),
                     "\n(r = ",round(cor_val[["estimate"]][["cor"]],2),
                     ", p = ",signif(cor_val[["p.value"]], digits = 3),")"),
      x = paste("Maternal",exposure,"log2FC",sep=" "),
      y = paste("Fetal",exposure,"log2FC",sep=" ")
    ) +
    theme_minimal() +
    ggrepel::geom_label_repel(aes(label = symbol), 
                              size = 3,
                              max.overlaps = 100,
                              force = 10)
  
  # Plot correlation of exposures
  
  cov_mat <- data.frame(result_mat[["dds"]]@colData@listData)
  cov_child <- data.frame(result_child[["dds"]]@colData@listData)
  
  cov_merged <- dplyr::inner_join(
    cov_mat %>%
      dplyr::select(SampleID, exposure_value = paste0(exposure,"_mat")) %>%
      dplyr::rename(PFBA_mat = exposure_value),
    
    cov_child %>%
      dplyr::select(SampleID, exposure_value = paste0(exposure,"_child")) %>%
      dplyr::rename(PFBA_child = exposure_value),
    
    by = "SampleID"
  )
  
  names(cov_merged)[2:3] <- c("cov_mat","cov_child")
  
  cor_val2 <- cor.test(cov_merged$cov_mat, cov_merged$cov_child, method = "pearson")
  
  g2 <- ggplot(cov_merged, aes(x = cov_mat, y = cov_child)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    labs(
      title = paste0(paste("Correlation of Maternal and Fetal Blood",exposure,sep=" "),
                     "\n(r = ",round(cor_val2[["estimate"]][["cor"]],2),
                     ", p = ",signif(cor_val2[["p.value"]], digits = 3),")"),
      x = paste("Maternal Blood",exposure,sep=" "),
      y = paste("Fetal Blood",exposure,sep=" ")
    ) +
    theme_minimal()
  
  g3 <- cowplot::plot_grid(g1,g2)
  
  ggsave(
    filename = paste0(output_dir, "/", exposure, "_MatChildCorr_DTE.pdf"),
    plot = g3,
    width = 8.8,
    height = 4.22,
    units="in"
  )
  
  return(g3)
}


all_pfas = colnames(pfas_combined)[c(2:9,11:18)]


## Maternal
result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[1], 
                          covariates = c("sex", "GA"), 
                          k = 6, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[1], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[1],output_dir = "results_gencode")
result_PFBA_mat <- result


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[2], 
                          covariates = c("sex", "GA"), 
                          k = 0, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[2], output_dir = "results_gencode")
csv_out(result, pfas_var = all_pfas[2],output_dir = "results_gencode")
result_PFOA_mat <- result

result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[3], 
                          covariates = c("sex", "GA"), 
                          k = 2, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[3], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[3],output_dir = "results_gencode")
result_PFNA_mat <- result


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[4], 
                          covariates = c("sex", "GA"), 
                          k = 1, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[4], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[4],output_dir = "results_gencode")
result_PFDA_mat <- result

result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[5], 
                          covariates = c("sex", "GA"), 
                          k = 1, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[5], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[5],output_dir = "results_gencode")
result_PFUnDA_mat <- result

result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[6], 
                          covariates = c("sex", "GA"), 
                          k = 4, output_dir = "results_gencode", 
                          just_hist = F, fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[6], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[6],output_dir = "results_gencode")
result_PFBS_mat <- result



result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[7], 
                          covariates = c("sex", "GA"), 
                          k = 1, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[7], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[7],output_dir = "results_gencode")
result_PFHxS_mat <- result


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[8], 
                          covariates = c("sex", "GA"), 
                          k = 1, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[8], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[8],output_dir = "results_gencode")
result_PFOS_mat <- result


## Child (with mat-child correlations)
result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[9], 
                          covariates = c("sex", "GA"), 
                          k = 3, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[9], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[9],output_dir = "results_gencode")

result_PFBA_child <- result
vis_efcor(result_PFBA_mat,result_PFBA_child,all_pfas[9],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[10], 
                          covariates = c("sex", "GA"), 
                          k = 1, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[10], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[10],output_dir = "results_gencode")

result_PFOA_child <- result
vis_efcor(result_PFOA_mat,result_PFOA_child,all_pfas[10],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[11], 
                          covariates = c("sex", "GA"), 
                          k = 5, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[11], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[11],output_dir = "results_gencode")

result_PFNA_child <- result
vis_efcor(result_PFNA_mat,result_PFNA_child,all_pfas[11],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[12], 
                          covariates = c("sex", "GA"), 
                          k = 3, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[12], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[12],output_dir = "results_gencode")

result_PFDA_child <- result
vis_efcor(result_PFDA_mat,result_PFDA_child,all_pfas[12],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[13], 
                          covariates = c("sex", "GA"), 
                          k = 5, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[13], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[13],output_dir = "results_gencode")

result_PFUnDA_child <- result
vis_efcor(result_PFUnDA_mat,result_PFUnDA_child,all_pfas[13],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[14], 
                          covariates = c("sex", "GA"), 
                          k = 4, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[14], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[14],output_dir = "results_gencode")

result_PFBS_child <- result
vis_efcor(result_PFBS_mat,result_PFBS_child,all_pfas[14],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[15], 
                          covariates = c("sex", "GA"), 
                          k = 10, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[15], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[15],output_dir = "results_gencode")

result_PFHxS_child <- result
vis_efcor(result_PFHxS_mat,result_PFHxS_child,all_pfas[15],"results_gencode")


result <- run_ruvr_deseq2(se = se_gusto.tx, 
                          pfas_var = all_pfas[16], 
                          covariates = c("sex", "GA"), 
                          k = 10, output_dir = "results_gencode", 
                          just_hist = FALSE,fdr_threshold = .1,
                          transcript_to_gene = tx2g.assembly)
volcano <- create_volcano_plot(result, pfas_var = all_pfas[16], output_dir = "results_gencode")
csv_out(result,pfas_var = all_pfas[16],output_dir = "results_gencode")

result_PFOS_child <- result
vis_efcor(result_PFOS_mat,result_PFOS_child,all_pfas[16],"results_gencode")
