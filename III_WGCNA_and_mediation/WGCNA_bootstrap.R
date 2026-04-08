#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript WGCNA_bootstrap_iteration.R <iteration_number>")
}

iteration_number <- as.numeric(args[1])

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
library(WGCNA)

setwd("/rsrch5/scratch/epi/stbresnahan/PFAS")

load("input_for_WGCNA.RData")
load("WGCNA_parameter_grid.RData")


correlate_modules_traits <- function(datExpr = NULL, 
                                     moduleColors = NULL, 
                                     datTraits, 
                                     cols = NULL,
                                     MEs = NULL,
                                     plot_heatmap = TRUE,
                                     plot_individual = FALSE,
                                     cor_method = "pearson",
                                     p_adjust_method = "BH",
                                     min_module_size = 10,
                                     output_dir = NULL) {
  
  # Input validation - check if MEs provided or need to calculate
  if (is.null(MEs)) {
    if (is.null(datExpr) || is.null(moduleColors)) {
      stop("Either provide MEs directly, or provide both datExpr and moduleColors")
    }
    
    if (!is.matrix(datExpr) && !is.data.frame(datExpr)) {
      stop("datExpr must be a matrix or data frame")
    }
    
    if (length(moduleColors) != ncol(datExpr)) {
      stop("Length of moduleColors must equal number of genes in datExpr")
    }
    
    if (nrow(datExpr) != nrow(datTraits)) {
      stop("Number of samples in datExpr must equal number of rows in datTraits")
    }
    
    calculate_MEs <- TRUE
  } else {
    if (nrow(MEs) != nrow(datTraits)) {
      stop("Number of samples in MEs must equal number of rows in datTraits")
    }
    calculate_MEs <- FALSE
    cat("Using pre-calculated module eigengenes\n")
  }
  
  # Handle column selection
  if (is.null(cols)) {
    cols <- colnames(datTraits)
    cat("Analyzing all", length(cols), "traits in datTraits\n")
  } else {
    missing_cols <- setdiff(cols, colnames(datTraits))
    if (length(missing_cols) > 0) {
      stop("Columns not found in datTraits: ", paste(missing_cols, collapse = ", "))
    }
    cat("Analyzing", length(cols), "specified traits:", paste(cols, collapse = ", "), "\n")
  }
  
  # Select trait columns and handle different data types
  traits_selected <- datTraits[, cols, drop = FALSE]
  
  # Convert factors to numeric if possible
  for (col in colnames(traits_selected)) {
    if (is.factor(traits_selected[[col]])) {
      if (nlevels(traits_selected[[col]]) == 2) {
        # Binary factor - convert to 0/1
        traits_selected[[col]] <- as.numeric(traits_selected[[col]]) - 1
        cat("Converted binary factor", col, "to numeric (0/1)\n")
      } else if (nlevels(traits_selected[[col]]) > 2) {
        # Multi-level factor - create dummy variables
        cat("Multi-level factor detected for", col, "- creating dummy variables\n")
        dummy_vars <- model.matrix(~ . - 1, data = traits_selected[col])
        colnames(dummy_vars) <- paste0(col, "_", colnames(dummy_vars))
        
        # Remove original column and add dummy variables
        traits_selected <- traits_selected[, !colnames(traits_selected) %in% col, drop = FALSE]
        traits_selected <- cbind(traits_selected, dummy_vars)
      }
    } else if (is.character(traits_selected[[col]])) {
      # Convert character to factor then numeric
      traits_selected[[col]] <- as.factor(traits_selected[[col]])
      traits_selected[[col]] <- as.numeric(traits_selected[[col]]) - 1
      cat("Converted character", col, "to numeric\n")
    }
  }
  
  # Basic statistics
  if (calculate_MEs) {
    nGenes <- ncol(datExpr)
    nSamples <- nrow(datExpr)
  } else {
    nSamples <- nrow(MEs)
    nGenes <- "Not calculated (using pre-calculated MEs)"
  }
  nTraits <- ncol(traits_selected)
  
  cat("Dataset info:\n")
  cat("- Genes:", nGenes, "\n")
  cat("- Samples:", nSamples, "\n") 
  cat("- Traits to analyze:", nTraits, "\n")
  
  # Calculate or use provided module eigengenes
  if (calculate_MEs) {
    cat("Calculating module eigengenes...\n")
    MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs <- orderMEs(MEs0)
    
    # Filter small modules if requested
    module_sizes <- table(moduleColors)
    large_modules <- names(module_sizes)[module_sizes >= min_module_size]
    
    if (min_module_size > 1) {
      # Keep only eigengenes for large modules (plus grey if present)
      keep_MEs <- paste0("ME", large_modules)
      if ("MEgrey" %in% names(MEs)) keep_MEs <- c(keep_MEs, "MEgrey")
      MEs <- MEs[, names(MEs) %in% keep_MEs, drop = FALSE]
      cat("Filtered to", ncol(MEs), "modules with >=", min_module_size, "genes\n")
    }
  } else {
    cat("Using provided module eigengenes\n")
    # Ensure MEs are properly ordered
    MEs <- orderMEs(MEs)
  }
  
  nModules <- ncol(MEs)
  modNames <- substring(names(MEs), 3)
  
  cat("Analyzing", nModules, "modules:", paste(modNames, collapse = ", "), "\n")
  
  # Module-trait correlations
  cat("Computing module-trait correlations...\n")
  moduleTraitCor <- cor(MEs, traits_selected, use = "pairwise.complete.obs", method = cor_method)
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Adjust p-values for multiple testing
  moduleTraitPvalue_adj <- matrix(p.adjust(as.vector(moduleTraitPvalue), method = p_adjust_method),
                                  nrow = nrow(moduleTraitPvalue),
                                  ncol = ncol(moduleTraitPvalue))
  rownames(moduleTraitPvalue_adj) <- rownames(moduleTraitPvalue)
  colnames(moduleTraitPvalue_adj) <- colnames(moduleTraitPvalue)
  
  # Gene-module membership (MM) - only calculate if we have expression data
  if (calculate_MEs) {
    cat("Computing gene-module membership...\n")
    geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs", method = cor_method))
    MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
    names(geneModuleMembership) <- paste("MM", modNames, sep = "")
    names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  } else {
    cat("Skipping gene-module membership (no expression data provided)\n")
    geneModuleMembership <- NULL
    MMPvalue <- NULL
  }
  
  # Gene-trait significance (GS) for each trait - only if we have expression data
  if (calculate_MEs) {
    cat("Computing gene-trait significance...\n")
    geneTraitSignificance <- list()
    GSPvalue <- list()
    
    for (trait in colnames(traits_selected)) {
      trait_clean <- make.names(trait)  # Clean trait name for R naming conventions
      
      gs_cor <- cor(datExpr, traits_selected[, trait, drop = FALSE], 
                    use = "pairwise.complete.obs", method = cor_method)
      gs_pval <- corPvalueStudent(as.matrix(gs_cor), nSamples)
      
      geneTraitSignificance[[trait_clean]] <- as.data.frame(gs_cor)
      GSPvalue[[trait_clean]] <- as.data.frame(gs_pval)
      
      names(geneTraitSignificance[[trait_clean]]) <- paste("GS", trait_clean, sep = ".")
      names(GSPvalue[[trait_clean]]) <- paste("p.GS", trait_clean, sep = ".")
    }
  } else {
    cat("Skipping gene-trait significance (no expression data provided)\n")
    geneTraitSignificance <- NULL
    GSPvalue <- NULL
  }
  
  # Create summary statistics
  significant_correlations <- sum(moduleTraitPvalue_adj < 0.05)
  strong_correlations <- sum(abs(moduleTraitCor) > 0.5 & moduleTraitPvalue_adj < 0.05)
  
  cat("Summary:\n")
  cat("- Significant correlations (adj. p < 0.05):", significant_correlations, "out of", nModules * nTraits, "\n")
  cat("- Strong significant correlations (|r| > 0.5, adj. p < 0.05):", strong_correlations, "\n")
  
  # Create plots
  if (plot_heatmap) {
    cat("Creating correlation heatmap ordered to separate positive/negative correlations...\n")
    
    # Strategy: Order modules and traits to group positive and negative correlations
    # Calculate mean correlation (not absolute) to distinguish positive vs negative patterns
    mean_cors_modules <- rowMeans(moduleTraitCor)
    mean_cors_traits <- colMeans(moduleTraitCor)
    
    # Order modules: strongest positive first, then strongest negative
    positive_modules <- which(mean_cors_modules >= 0)
    negative_modules <- which(mean_cors_modules < 0)
    
    # Within positive modules, order by decreasing mean correlation
    pos_order <- positive_modules[order(mean_cors_modules[positive_modules], decreasing = TRUE)]
    # Within negative modules, order by increasing mean correlation (most negative first)
    neg_order <- negative_modules[order(mean_cors_modules[negative_modules], decreasing = FALSE)]
    
    module_order <- c(pos_order, neg_order)
    
    # Same strategy for traits
    positive_traits <- which(mean_cors_traits >= 0)
    negative_traits <- which(mean_cors_traits < 0)
    
    pos_trait_order <- positive_traits[order(mean_cors_traits[positive_traits], decreasing = TRUE)]
    neg_trait_order <- negative_traits[order(mean_cors_traits[negative_traits], decreasing = FALSE)]
    
    trait_order <- c(pos_trait_order, neg_trait_order)
    
    # Reorder correlation matrix based on correlation strength
    moduleTraitCor_ordered <- moduleTraitCor[module_order, trait_order, drop = FALSE]
    moduleTraitPvalue_adj_ordered <- moduleTraitPvalue_adj[module_order, trait_order, drop = FALSE]
    
    # Prepare text matrix for ordered heatmap
    textMatrix_ordered <- paste(signif(moduleTraitCor_ordered, 2), "\n(",
                                signif(moduleTraitPvalue_adj_ordered, 1), ")", sep = "")
    dim(textMatrix_ordered) <- dim(moduleTraitCor_ordered)
    
    # Calculate optimal plot dimensions
    heatmap_width <- max(8, min(20, 2 + nTraits * 1.2))
    heatmap_height <- max(6, min(16, 2 + nModules * 0.4))
    
    # Create directory if specified
    if (!is.null(output_dir) && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Create ordered heatmap
    if (!is.null(output_dir)) {
      pdf(file.path(output_dir, "module_trait_heatmap_ordered.pdf"), 
          width = heatmap_width, height = heatmap_height)
    } else {
      sizeGrWindow(heatmap_width, heatmap_height)
    }
    
    par(mar = c(8, 10, 4, 4))
    
    # Create ordered labels
    ordered_module_labels <- rownames(moduleTraitCor)[module_order]
    ordered_trait_labels <- colnames(moduleTraitCor)[trait_order]
    
    # Create ordered heatmap
    labeledHeatmap(Matrix = moduleTraitCor_ordered,
                   xLabels = ordered_trait_labels,
                   yLabels = ordered_module_labels,
                   ySymbols = ordered_module_labels,
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix_ordered,
                   setStdMargins = FALSE,
                   cex.text = max(0.4, min(0.8, 15/max(nModules, nTraits))),
                   cex.lab = 1.0,
                   zlim = c(-1, 1),
                   legendLabel = paste("Correlation (", cor_method, ")", sep = ""),
                   main = paste("Module-Trait Relationships (Adjusted p-values,", p_adjust_method, "method)"))
    
    if (!is.null(output_dir)) {
      dev.off()
      cat("Ordered heatmap saved to:", file.path(output_dir, "module_trait_heatmap_ordered.pdf"), "\n")
    }
  }
  
  # Individual trait plots (optional)
  if (plot_individual) {
    cat("Creating individual trait correlation plots...\n")
    
    for (i in seq_along(cols)) {
      trait_name <- cols[i]
      trait_clean <- make.names(trait_name)
      
      if (!is.null(output_dir)) {
        pdf(file.path(output_dir, paste0("trait_", trait_clean, "_correlations.pdf")), 
            width = 10, height = 6)
      } else {
        sizeGrWindow(10, 6)
      }
      
      # Plot module correlations for this trait
      par(mfrow = c(1, 2))
      
      # Module correlations
      trait_cors <- moduleTraitCor[, trait_name]
      trait_pvals <- moduleTraitPvalue_adj[, trait_name]
      
      barplot(trait_cors, 
              names.arg = modNames,
              main = paste("Module correlations with", trait_name),
              ylab = paste("Correlation (", cor_method, ")", sep = ""),
              col = ifelse(trait_pvals < 0.05, "red", "grey"),
              las = 2)
      abline(h = c(-0.5, 0.5), lty = 2, col = "blue")
      legend("topright", legend = c("p.adj < 0.05", "p.adj >= 0.05"), 
             fill = c("red", "grey"))
      
      # P-value plot
      barplot(-log10(trait_pvals), 
              names.arg = modNames,
              main = paste("-log10(adj. p-value) for", trait_name),
              ylab = "-log10(adj. p-value)",
              col = ifelse(trait_pvals < 0.05, "red", "grey"),
              las = 2)
      abline(h = -log10(0.05), lty = 2, col = "blue")
      
      if (!is.null(output_dir)) {
        dev.off()
      }
    }
  }
  
  # Return comprehensive results
  results <- list(
    # Main results
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    moduleTraitPvalue_adj = moduleTraitPvalue_adj,
    
    # Module eigengenes
    MEs = MEs,
    
    # Gene-level results
    geneModuleMembership = geneModuleMembership,
    MMPvalue = MMPvalue,
    geneTraitSignificance = geneTraitSignificance,
    GSPvalue = GSPvalue,
    
    # Metadata
    moduleColors = if(calculate_MEs) moduleColors else NULL,
    traits_analyzed = colnames(traits_selected),
    cor_method = cor_method,
    p_adjust_method = p_adjust_method,
    nSamples = nSamples,
    nGenes = nGenes,
    nModules = nModules,
    nTraits = nTraits,
    
    # Summary statistics
    summary_stats = list(
      significant_correlations = significant_correlations,
      strong_correlations = strong_correlations,
      module_sizes = if(calculate_MEs) table(moduleColors) else "Not available"
    )
  )
  
  cat("Analysis complete!\n")
  return(results)
}


create_numeric_module_mapping <- function(wgcna_result) {
  
  # Extract colors and MEs
  moduleColors <- wgcna_result$colors
  MEs <- wgcna_result$MEs
  
  # Get unique color labels and their frequencies
  color_info <- table(moduleColors)
  unique_colors <- names(color_info)
  
  # Get ME labels (extract numeric part)
  ME_labels <- colnames(MEs)
  numeric_labels <- as.numeric(gsub("ME", "", ME_labels))
  
  # Create the mapping using WGCNA's labels2colors function
  # The key insight: WGCNA assigns colors based on module size order
  # Module 0 = grey (unassigned), then largest module = turquoise, etc.
  
  # Method 1: Direct mapping using the color vector
  # Find which numeric label corresponds to each color
  
  # Get all numeric labels present in the data
  if ("numericLabels" %in% names(wgcna_result)) {
    # If numeric labels were saved
    numeric_assignments <- wgcna_result$numericLabels
  } else {
    # Convert colors back to numeric using WGCNA standard
    library(WGCNA)
    # This recreates the numeric assignments
    color_to_numeric <- c("grey" = 0)
    standard_colors <- standardColors()
    
    # Find which colors are actually used (excluding grey)
    used_colors <- unique_colors[unique_colors != "grey"]
    
    # Sort by module size (WGCNA assigns numbers by size)
    color_sizes <- color_info[used_colors]
    sorted_colors <- names(sort(color_sizes, decreasing = TRUE))
    
    # Assign numbers 1, 2, 3, ... to colors by size
    for (i in seq_along(sorted_colors)) {
      color_to_numeric[sorted_colors[i]] <- i
    }
    
    numeric_assignments <- color_to_numeric[moduleColors]
  }
  
  # Create comprehensive mapping
  mapping_list <- list()
  
  for (me_label in ME_labels) {
    numeric_label <- as.numeric(gsub("ME", "", me_label))
    
    # Find corresponding color
    corresponding_color <- names(color_info)[sapply(names(color_info), function(color) {
      any(numeric_assignments == numeric_label & moduleColors == color)
    })]
    
    if (length(corresponding_color) == 0) {
      # Handle edge case where no genes assigned to this module
      corresponding_color <- "unknown"
      module_size <- 0
    } else {
      corresponding_color <- corresponding_color[1]  # Take first match
      module_size <- color_info[corresponding_color]
    }
    
    mapping_list[[me_label]] <- data.frame(
      ME_label = me_label,
      numeric_label = numeric_label,
      module_color = corresponding_color,
      module_size = as.numeric(module_size),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine into single dataframe
  mapping_df <- do.call(rbind, mapping_list)
  rownames(mapping_df) <- NULL
  
  # Sort by numeric label for clarity
  mapping_df <- mapping_df[order(mapping_df$numeric_label), ]
  
  # Add validation
  mapping_df$genes_available <- !is.na(mapping_df$module_size) & mapping_df$module_size > 0
  
  cat("=== MODULE MAPPING CREATED ===\n")
  cat("Total modules:", nrow(mapping_df), "\n")
  cat("Total genes:", length(moduleColors), "\n\n")
  
  # Show the mapping
  print(mapping_df)
  
  return(mapping_df)
}

extract_significant_associations <- function(correlation_results, 
                                             p_threshold = 0.05, 
                                             cor_threshold = 0.3) {
  
  cor_matrix <- correlation_results$moduleTraitCor
  pval_matrix <- correlation_results$moduleTraitPvalue
  
  # Find significant associations
  significant_idx <- which(abs(cor_matrix) >= cor_threshold & 
                             pval_matrix <= p_threshold, arr.ind = TRUE)
  
  if (nrow(significant_idx) == 0) {
    cat("No significant associations found with current thresholds\n")
    return(data.frame())
  }
  
  # Create results data frame
  results_df <- data.frame(
    Module = rownames(cor_matrix)[significant_idx[, 1]],
    Trait = colnames(cor_matrix)[significant_idx[, 2]],
    Correlation = cor_matrix[significant_idx],
    P_value = correlation_results$moduleTraitPvalue[significant_idx],
    Adj_P_value = pval_matrix[significant_idx],
    stringsAsFactors = FALSE
  )
  
  # Sort by absolute correlation
  results_df <- results_df[order(abs(results_df$Correlation), decreasing = TRUE), ]
  
  cat("Found", nrow(results_df), "significant module-trait associations\n")
  return(results_df)
}

get_genes_by_module_id <- function(wgcna_result, module_id, gene_names = NULL, mapping=NULL) {
  
  moduleColors <- wgcna_result$colors
  
  # Determine what type of identifier was provided
  if (grepl("^ME", module_id)) {
    # ME label provided (e.g., "ME5")
    target_numeric <- as.numeric(gsub("ME", "", module_id))
    # Need to find corresponding color
    if(is.null(mapping)){
      mapping <- create_numeric_module_mapping(wgcna_result)
    }
    target_color <- mapping$module_color[mapping$numeric_label == target_numeric]
  } else if (is.numeric(module_id)) {
    # Numeric label provided
    if(is.null(mapping)){
      mapping <- create_numeric_module_mapping(wgcna_result)
    }
    target_color <- mapping$module_color[mapping$numeric_label == module_id]
  } else {
    # Color provided
    target_color <- module_id
  }
  
  if (length(target_color) == 0 || !target_color %in% moduleColors) {
    stop("Module identifier '", module_id, "' not found")
  }
  
  # Get gene indices
  gene_indices <- which(moduleColors == target_color)
  
  if (is.null(gene_names)) {
    return(gene_indices)
  } else {
    return(gene_names[gene_indices])
  }
}

run_optimized_wgcna <- function(datExpr, 
                                power = NULL,
                                minModuleSize = 30,
                                deepSplit = 3,
                                mergeCutHeight = 0.25,
                                detectCutHeight = NULL,
                                maxBlockSize = 5000) {
  
  cat("=== Running Optimized WGCNA ===\n")
  
  # Auto-select power if not provided
  if (is.null(power)) {
    sft <- pickSoftThreshold(datExpr, networkType = "signed", verbose = 5)
    power <- ifelse(is.na(sft$powerEstimate), 10, sft$powerEstimate)
  }
  
  cat("Using soft-thresholding power:", power, "\n")
  
  # Option 1: Step-by-step approach (your current method)
  step_by_step_approach <- function() {
    cat("Building adjacency matrix...\n")
    adj_matrix <- adjacency(datExpr, power = power, type = "signed")
    
    cat("Computing TOM similarity...\n")
    TOM <- TOMsimilarity(adj_matrix, TOMType = "signed")  # Match network type
    dissTOM <- 1 - TOM
    
    cat("Hierarchical clustering...\n")
    geneTree <- hclust(as.dist(dissTOM), method = "average")
    
    cat("Dynamic module detection...\n")
    dynamicMods <- cutreeDynamic(dendro = geneTree, 
                                 distM = dissTOM,
                                 deepSplit = deepSplit, 
                                 pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize,
                                 cutHeight = detectCutHeight)
    
    dynamicColors <- labels2colors(dynamicMods)
    
    cat("Module merging...\n")
    MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
    MEs <- MEList$eigengenes
    MEDiss <- 1 - cor(MEs)
    METree <- hclust(as.dist(MEDiss), method = "average")
    
    merge <- mergeCloseModules(datExpr, dynamicColors, 
                               cutHeight = mergeCutHeight, verbose = 3)
    
    return(list(
      colors = merge$colors,
      MEs = merge$newMEs,
      geneTree = geneTree,
      dynamicColors = dynamicColors,
      TOM = TOM
    ))
  }
  
  # Option 2: blockwiseModules approach (recommended for large datasets)
  blockwise_approach <- function() {
    cat("Using blockwiseModules approach...\n")
    
    net <- blockwiseModules(datExpr,
                            power = power,
                            networkType = "signed",
                            TOMType = "signed",
                            minModuleSize = minModuleSize,
                            deepSplit = deepSplit,
                            pamRespectsDendro = FALSE,
                            mergeCutHeight = mergeCutHeight,
                            numericLabels = TRUE,
                            saveTOMs = FALSE,  # Set TRUE if you want to save TOM
                            maxBlockSize = maxBlockSize,
                            verbose = 3)
    
    return(list(
      colors = labels2colors(net$colors),
      MEs = net$MEs,
      dendrograms = net$dendrograms,
      blockGenes = net$blockGenes
    ))
  }
  
  # Choose approach based on dataset size
  n_genes <- ncol(datExpr)
  if (n_genes > maxBlockSize) {
    cat("Large dataset detected. Using blockwise approach.\n")
    result <- blockwise_approach()
  } else {
    cat("Using step-by-step approach.\n")
    result <- step_by_step_approach()
  }
  
  # Summary statistics
  cat("\n=== Results Summary ===\n")
  module_counts <- table(result$colors)
  cat("Number of modules:", length(module_counts), "\n")
  cat("Module sizes:\n")
  print(sort(module_counts, decreasing = TRUE))
  
  grey_genes <- sum(result$colors == "grey")
  cat("Genes not assigned to modules (grey):", grey_genes, 
      "(", round(100 * grey_genes / n_genes, 1), "%)\n")
  
  return(result)
}


if (iteration_number < 1 || iteration_number > nrow(param_grid)) {
  stop("Invalid iteration number. Must be between 1 and ", nrow(param_grid))
}

cat("\n=== ITERATION", iteration_number, "of", nrow(param_grid), "===\n")
cat("Parameters:", 
    "power =", param_grid$power[iteration_number],
    "minModuleSize =", param_grid$minModuleSize[iteration_number],
    "deepSplit =", param_grid$deepSplit[iteration_number],
    "mergeCutHeight =", param_grid$mergeCutHeight[iteration_number], "\n")

tryCatch({
  WGCNA.dat <- run_optimized_wgcna(
    t(expr_PFAS), 
    power = param_grid$power[iteration_number],
    minModuleSize = param_grid$minModuleSize[iteration_number],
    deepSplit = param_grid$deepSplit[iteration_number],
    mergeCutHeight = param_grid$mergeCutHeight[iteration_number]
  )
  
  color_ME_mapping <- create_numeric_module_mapping(WGCNA.dat)
  
  sids <- names(data.frame(expr_PFAS))
  covs_PFAS_subset <- covs_PFAS[covs_PFAS$SampleID %in% sids, ]
  
  cors <- correlate_modules_traits(
    t(expr_PFAS),
    moduleColors = WGCNA.dat$colors,
    datTraits = covs_PFAS_subset,
    cols = names(covs_PFAS)[c(53:60, 62:69)],
    MEs = WGCNA.dat$MEs,
    plot_heatmap = FALSE,
    plot_individual = FALSE
  )
  
  sig.cors <- extract_significant_associations(cors, cor_threshold = 0.2)
  
  output_file <- paste0("/rsrch5/home/epi/stbresnahan/scratch/PFAS/bootstrapped/WGCNA_output_", 
                        iteration_number, ".RData")
  save(
    file = output_file,
    list = c("WGCNA.dat", "color_ME_mapping", "cors", "sig.cors", "get_genes_by_module_id")
  )
  
  cat("Saved results to:", output_file, "\n")
  
}, error = function(e) {
  cat("ERROR in iteration", iteration_number, ":", conditionMessage(e), "\n")
  quit(status = 1)
})

cat("\n=== ITERATION", iteration_number, "COMPLETE ===\n")