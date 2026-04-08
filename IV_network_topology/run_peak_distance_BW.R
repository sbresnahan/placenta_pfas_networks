#!/usr/bin/env Rscript

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
setwd("/rsrch5/scratch/epi/stbresnahan/PFAS")

suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
})

## ---- Parse command line argument ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript run_peak_distance.R <iteration>")
iter_val <- as.integer(args[1])
message("Running iteration: ", iter_val)

message("Loading files...")
adjacency <- readRDS("/rsrch5/home/epi/stbresnahan/scratch/PFAS_mediation/adjacency.rds")
load("/rsrch5/home/epi/stbresnahan/scratch/PFAS/data_for_distance_analysis.RDS")  # loads all_module_mediation_df, pfas_order, tpte_values

## ---- Subset modules for this iteration ----
dt_mediation <- as.data.table(all_module_mediation_df)
rm(all_module_mediation_df)
gc()
dt_mediation <- dt_mediation[iteration == iter_val]

modules <- unique(dt_mediation$module)

## ---- Peak distance function using data.table ----
compute_peak_distance_dt <- function(modules, pfas_order, adjacency, dt_mediation, tpte_values) {
  
  get_peak_or_mean <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(mean(x, na.rm = TRUE))
    d <- density(x, na.rm = TRUE)
    d$x[which.max(d$y)]
  }
  
  results <- vector("list", length = 0)
  res_idx <- 1
  
  for (module_name in modules) {
    message("Processing module: ", module_name)
    
    # Genes in module
    genes_in_module <- dt_mediation[module == module_name, unique(transcript_id)]
    
    # Subset adjacency and threshold weak edges
    me_adj <- adjacency[genes_in_module, genes_in_module, drop = FALSE]
    me_adj[me_adj < 0.1] <- 0
    
    # Keep only nodes with ≥1 edge
    connected_genes <- rowSums(me_adj > 0) > 0
    me_adj <- me_adj[connected_genes, connected_genes, drop = FALSE]
    
    if (nrow(me_adj) < 3) next
    g <- graph_from_adjacency_matrix(me_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    
    # Identify hub
    centrality_vals <- eigen_centrality(g, weights = E(g)$weight)$vector
    most_central_gene <- names(which.max(centrality_vals))
    
    for (this_exposure_window in c("child", "mat")) {
      for (pfas in pfas_order) {
        
        # Get significant mediators (data.table logic)
        sig_genes <- dt_mediation[
          module == module_name &
            exposure_window == this_exposure_window &
            PFAS == pfas &
            significant_strict == TRUE &
            transcript_id %in% V(g)$name,
          unique(transcript_id)
        ]
        
        if (length(sig_genes) == 0) next
        
        # Compute graph distances from hub to sig genes
        dist_matrix <- distances(g, v = most_central_gene, to = sig_genes, weights = 1 / E(g)$weight)
        distances_vec <- as.numeric(dist_matrix)
        
        # Compute peak (or mean if too few)
        peak_distance <- get_peak_or_mean(distances_vec)
        
        # Store result
        results[[res_idx]] <- data.table(
          module = module_name,
          PFAS = pfas,
          exposure_label = ifelse(this_exposure_window == "child", "Fetal", "Maternal"),
          central_hub = most_central_gene,
          peak_distance = peak_distance
        )
        res_idx <- res_idx + 1
      }
    }
  }
  
  result_dt <- rbindlist(results)
  result_dt <- merge(result_dt, tpte_values, by = "PFAS", all.x = TRUE)
  setorder(result_dt, exposure_label, PFAS)
  
  return(result_dt)
}

## ---- Run function ----
transcript_peak_distances <- compute_peak_distance_dt(
  modules = modules,
  pfas_order = pfas_order,
  adjacency = adjacency,
  dt_mediation = dt_mediation,
  tpte_values = tpte_values
)

## ---- Save output ----
out_dir <- "/rsrch5/home/epi/stbresnahan/scratch/PFAS/peak_distances"
out_file <- file.path(out_dir, paste0("transcript_peak_distances_", iter_val, ".rds"))
saveRDS(transcript_peak_distances, out_file)
message("Saved: ", out_file)
