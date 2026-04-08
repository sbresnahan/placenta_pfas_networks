.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
library(WGCNA)

load("input_for_WGCNA.RData")

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
    TOM <- TOMsimilarity(adj_matrix, TOMType = "signed")
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
                            saveTOMs = FALSE,
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
  
  traits_selected <- datTraits[, cols, drop = FALSE]
  
  for (col in colnames(traits_selected)) {
    if (is.factor(traits_selected[[col]])) {
      if (nlevels(traits_selected[[col]]) == 2) {
        traits_selected[[col]] <- as.numeric(traits_selected[[col]]) - 1
        cat("Converted binary factor", col, "to numeric (0/1)\n")
      } else if (nlevels(traits_selected[[col]]) > 2) {
        cat("Multi-level factor detected for", col, "- creating dummy variables\n")
        dummy_vars <- model.matrix(~ . - 1, data = traits_selected[col])
        colnames(dummy_vars) <- paste0(col, "_", colnames(dummy_vars))
        traits_selected <- traits_selected[, !colnames(traits_selected) %in% col, drop = FALSE]
        traits_selected <- cbind(traits_selected, dummy_vars)
      }
    } else if (is.character(traits_selected[[col]])) {
      traits_selected[[col]] <- as.factor(traits_selected[[col]])
      traits_selected[[col]] <- as.numeric(traits_selected[[col]]) - 1
      cat("Converted character", col, "to numeric\n")
    }
  }
  
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
  
  if (calculate_MEs) {
    cat("Calculating module eigengenes...\n")
    MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs <- orderMEs(MEs0)
    
    module_sizes <- table(moduleColors)
    large_modules <- names(module_sizes)[module_sizes >= min_module_size]
    
    if (min_module_size > 1) {
      keep_MEs <- paste0("ME", large_modules)
      if ("MEgrey" %in% names(MEs)) keep_MEs <- c(keep_MEs, "MEgrey")
      MEs <- MEs[, names(MEs) %in% keep_MEs, drop = FALSE]
      cat("Filtered to", ncol(MEs), "modules with >=", min_module_size, "genes\n")
    }
  } else {
    cat("Using provided module eigengenes\n")
    MEs <- orderMEs(MEs)
  }
  
  nModules <- ncol(MEs)
  modNames <- substring(names(MEs), 3)
  
  cat("Analyzing", nModules, "modules:", paste(modNames, collapse = ", "), "\n")
  
  cat("Computing module-trait correlations...\n")
  moduleTraitCor <- cor(MEs, traits_selected, use = "pairwise.complete.obs", method = cor_method)
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  moduleTraitPvalue_adj <- matrix(p.adjust(as.vector(moduleTraitPvalue), method = p_adjust_method),
                                  nrow = nrow(moduleTraitPvalue),
                                  ncol = ncol(moduleTraitPvalue))
  rownames(moduleTraitPvalue_adj) <- rownames(moduleTraitPvalue)
  colnames(moduleTraitPvalue_adj) <- colnames(moduleTraitPvalue)
  
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
  
  if (calculate_MEs) {
    cat("Computing gene-trait significance...\n")
    geneTraitSignificance <- list()
    GSPvalue <- list()
    
    for (trait in colnames(traits_selected)) {
      trait_clean <- make.names(trait)
      
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
  
  significant_correlations <- sum(moduleTraitPvalue_adj < 0.05)
  strong_correlations <- sum(abs(moduleTraitCor) > 0.5 & moduleTraitPvalue_adj < 0.05)
  
  cat("Summary:\n")
  cat("- Significant correlations (adj. p < 0.05):", significant_correlations, "out of", nModules * nTraits, "\n")
  cat("- Strong significant correlations (|r| > 0.5, adj. p < 0.05):", strong_correlations, "\n")
  
  if (plot_heatmap) {
    cat("Creating correlation heatmap ordered to separate positive/negative correlations...\n")
    
    mean_cors_modules <- rowMeans(moduleTraitCor)
    mean_cors_traits <- colMeans(moduleTraitCor)
    
    positive_modules <- which(mean_cors_modules >= 0)
    negative_modules <- which(mean_cors_modules < 0)
    
    pos_order <- positive_modules[order(mean_cors_modules[positive_modules], decreasing = TRUE)]
    neg_order <- negative_modules[order(mean_cors_modules[negative_modules], decreasing = FALSE)]
    
    module_order <- c(pos_order, neg_order)
    
    positive_traits <- which(mean_cors_traits >= 0)
    negative_traits <- which(mean_cors_traits < 0)
    
    pos_trait_order <- positive_traits[order(mean_cors_traits[positive_traits], decreasing = TRUE)]
    neg_trait_order <- negative_traits[order(mean_cors_traits[negative_traits], decreasing = FALSE)]
    
    trait_order <- c(pos_trait_order, neg_trait_order)
    
    moduleTraitCor_ordered <- moduleTraitCor[module_order, trait_order, drop = FALSE]
    moduleTraitPvalue_adj_ordered <- moduleTraitPvalue_adj[module_order, trait_order, drop = FALSE]
    
    textMatrix_ordered <- paste(signif(moduleTraitCor_ordered, 2), "\n(",
                                signif(moduleTraitPvalue_adj_ordered, 1), ")", sep = "")
    dim(textMatrix_ordered) <- dim(moduleTraitCor_ordered)
    
    heatmap_width <- max(8, min(20, 2 + nTraits * 1.2))
    heatmap_height <- max(6, min(16, 2 + nModules * 0.4))
    
    if (!is.null(output_dir) && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    if (!is.null(output_dir)) {
      pdf(file.path(output_dir, "module_trait_heatmap_ordered.pdf"), 
          width = heatmap_width, height = heatmap_height)
    } else {
      sizeGrWindow(heatmap_width, heatmap_height)
    }
    
    par(mar = c(8, 10, 4, 4))
    
    ordered_module_labels <- rownames(moduleTraitCor)[module_order]
    ordered_trait_labels <- colnames(moduleTraitCor)[trait_order]
    
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
      
      par(mfrow = c(1, 2))
      
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
  
  results <- list(
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    moduleTraitPvalue_adj = moduleTraitPvalue_adj,
    MEs = MEs,
    geneModuleMembership = geneModuleMembership,
    MMPvalue = MMPvalue,
    geneTraitSignificance = geneTraitSignificance,
    GSPvalue = GSPvalue,
    moduleColors = if(calculate_MEs) moduleColors else NULL,
    traits_analyzed = colnames(traits_selected),
    cor_method = cor_method,
    p_adjust_method = p_adjust_method,
    nSamples = nSamples,
    nGenes = nGenes,
    nModules = nModules,
    nTraits = nTraits,
    summary_stats = list(
      significant_correlations = significant_correlations,
      strong_correlations = strong_correlations,
      module_sizes = if(calculate_MEs) table(moduleColors) else "Not available"
    )
  )
  
  cat("Analysis complete!\n")
  return(results)
}

extract_significant_associations <- function(correlation_results, 
                                             p_threshold = 0.05, 
                                             cor_threshold = 0.3) {
  
  cor_matrix <- correlation_results$moduleTraitCor
  pval_matrix <- correlation_results$moduleTraitPvalue
  
  significant_idx <- which(abs(cor_matrix) >= cor_threshold & 
                             pval_matrix <= p_threshold, arr.ind = TRUE)
  
  if (nrow(significant_idx) == 0) {
    cat("No significant associations found with current thresholds\n")
    return(data.frame())
  }
  
  results_df <- data.frame(
    Module = rownames(cor_matrix)[significant_idx[, 1]],
    Trait = colnames(cor_matrix)[significant_idx[, 2]],
    Correlation = cor_matrix[significant_idx],
    P_value = correlation_results$moduleTraitPvalue[significant_idx],
    Adj_P_value = pval_matrix[significant_idx],
    stringsAsFactors = FALSE
  )
  
  results_df <- results_df[order(abs(results_df$Correlation), decreasing = TRUE), ]
  
  cat("Found", nrow(results_df), "significant module-trait associations\n")
  return(results_df)
}

create_numeric_module_mapping <- function(wgcna_result) {
  
  moduleColors <- wgcna_result$colors
  MEs <- wgcna_result$MEs
  
  color_info <- table(moduleColors)
  unique_colors <- names(color_info)
  
  ME_labels <- colnames(MEs)
  numeric_labels <- as.numeric(gsub("ME", "", ME_labels))
  
  if ("numericLabels" %in% names(wgcna_result)) {
    numeric_assignments <- wgcna_result$numericLabels
  } else {
    library(WGCNA)
    color_to_numeric <- c("grey" = 0)
    standard_colors <- standardColors()
    
    used_colors <- unique_colors[unique_colors != "grey"]
    
    color_sizes <- color_info[used_colors]
    sorted_colors <- names(sort(color_sizes, decreasing = TRUE))
    
    for (i in seq_along(sorted_colors)) {
      color_to_numeric[sorted_colors[i]] <- i
    }
    
    numeric_assignments <- color_to_numeric[moduleColors]
  }
  
  mapping_list <- list()
  
  for (me_label in ME_labels) {
    numeric_label <- as.numeric(gsub("ME", "", me_label))
    
    corresponding_color <- names(color_info)[sapply(names(color_info), function(color) {
      any(numeric_assignments == numeric_label & moduleColors == color)
    })]
    
    if (length(corresponding_color) == 0) {
      corresponding_color <- "unknown"
      module_size <- 0
    } else {
      corresponding_color <- corresponding_color[1]
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
  
  mapping_df <- do.call(rbind, mapping_list)
  rownames(mapping_df) <- NULL
  
  mapping_df <- mapping_df[order(mapping_df$numeric_label), ]
  
  mapping_df$genes_available <- !is.na(mapping_df$module_size) & mapping_df$module_size > 0
  
  cat("=== MODULE MAPPING CREATED ===\n")
  cat("Total modules:", nrow(mapping_df), "\n")
  cat("Total genes:", length(moduleColors), "\n\n")
  
  print(mapping_df)
  
  return(mapping_df)
}

get_genes_by_module_id <- function(wgcna_result, module_id, gene_names = NULL, mapping=NULL) {
  
  moduleColors <- wgcna_result$colors
  
  if (grepl("^ME", module_id)) {
    target_numeric <- as.numeric(gsub("ME", "", module_id))
    if(is.null(mapping)){
      mapping <- create_numeric_module_mapping(wgcna_result)
    }
    target_color <- mapping$module_color[mapping$numeric_label == target_numeric]
  } else if (is.numeric(module_id)) {
    if(is.null(mapping)){
      mapping <- create_numeric_module_mapping(wgcna_result)
    }
    target_color <- mapping$module_color[mapping$numeric_label == module_id]
  } else {
    target_color <- module_id
  }
  
  if (length(target_color) == 0 || !target_color %in% moduleColors) {
    stop("Module identifier '", module_id, "' not found")
  }
  
  gene_indices <- which(moduleColors == target_color)
  
  if (is.null(gene_names)) {
    return(gene_indices)
  } else {
    return(gene_names[gene_indices])
  }
}

powers <- 1:10

# Run soft-thresholding analysis
sft <- pickSoftThreshold(t(expr_PFAS), powerVector = powers, verbose = 5)

# Show mean connectivity for each power
mean_connectivity <- sft$fitIndices[, c("Power", "mean.k.")]
print("Mean connectivity by power:")
print(mean_connectivity)

# Optional: visualize scale-free fit and mean connectivity
par(mfrow = c(1,2))
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
     type = "b", xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     main = "Scale independence")
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     type = "b", xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity")

# Single run
WGCNA.dat <- run_optimized_wgcna(t(expr_PFAS), power=3, minModuleSize = 50)

sids <- names(data.frame(expr_PFAS))
covs_PFAS <- covs_PFAS[covs_PFAS$SampleID%in%sids,]

cors <- correlate_modules_traits(t(expr_PFAS),
                                 moduleColors=WGCNA.dat$colors,
                                 datTraits=covs_PFAS,
                                 cols=names(covs_PFAS)[c(53:60,62:69)],
                                 MEs=WGCNA.dat$MEs,
                                 output_dir="/rsrch5/scratch/epi/stbresnahan/PFAS")

sig.cors <- extract_significant_associations(cors,cor_threshold=0.2)

color_ME_mapping <- create_numeric_module_mapping(WGCNA.dat)

save(file="WGCNA_output.RData",list=c("WGCNA.dat","color_ME_mapping",
                                      "cors","sig.cors","get_genes_by_module_id"))

# Save functions for bootstrap script
save(
  file = "WGCNA_functions.RData",
  list = c("run_optimized_wgcna", 
           "correlate_modules_traits", 
           "extract_significant_associations", 
           "create_numeric_module_mapping", 
           "get_genes_by_module_id")
)

# Create parameter grid for bootstrapping
param_grid <- expand.grid(
  power = 3:5,
  minModuleSize = c(30, 50, 75, 100),
  deepSplit = 2:4,
  mergeCutHeight = c(0.15, 0.20, 0.25, 0.30)
)

param_grid <- cbind(iteration = 1:nrow(param_grid), param_grid)

save(param_grid, file = "WGCNA_parameter_grid.RData")
write.csv(param_grid, "WGCNA_parameter_grid.csv", row.names = FALSE)

cat("Created parameter grid with", nrow(param_grid), "iterations\n")