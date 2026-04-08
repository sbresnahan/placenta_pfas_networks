#!/usr/bin/env Rscript

.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))
setwd("/rsrch5/scratch/epi/stbresnahan/PFAS")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(igraph)
  library(data.table)
})

## ---- Parse command line argument ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_distance_analysis.R <iteration>")
}

iter_val <- as.integer(args[1])

message("Running iteration: ", iter_val)

message("Loading files...")

adjacency <- readRDS("/rsrch5/home/epi/stbresnahan/scratch/PFAS_mediation/adjacency_GA.rds")
load("/rsrch5/home/epi/stbresnahan/scratch/PFAS/data_for_distance_analysis_GA.RDS")

message("Running analysis...")

## ---- Subset modules for this iteration ----
# Convert once to data.table and filter iteration
dt_mediation <- as.data.table(all_module_mediation_df)
rm(all_module_mediation_df)
gc()
dt_mediation <- dt_mediation[iteration == iter_val]

# Get unique modules
modules <- unique(dt_mediation$module)

## ---- Run distance computation ----
compute_transcript_distances_dt <- function(modules, 
                                            pfas_order,
                                            adjacency,
                                            dt_mediation,
                                            tpte_values) {
  
  # Set key for fast lookups
  setkey(dt_mediation, module, exposure_window, PFAS, transcript_id)
  
  # Pre-split by module for fast gene lookup
  dt_by_module <- split(dt_mediation, by = "module", keep.by = FALSE)
  
  results <- vector("list", length = 0)  # store results
  res_idx <- 1
  
  counter <- 0
  
  for (module_name in modules) {
    
    counter <- counter+1
    
    message("Processing module ",counter," of ",length(modules))
    
    # Genes in module
    genes_in_module <- unique(dt_by_module[[module_name]]$transcript_id)
    
    # Subset adjacency matrix
    me_adj <- adjacency[genes_in_module, genes_in_module, drop = FALSE]
    me_adj[me_adj < 0.1] <- 0  # threshold weak edges
    
    # Build graph
    g <- graph_from_adjacency_matrix(me_adj, mode = "undirected", weighted = TRUE, diag = FALSE)
    if (vcount(g) < 3) next
    
    # Compute hub: eigenvector centrality (or switch to strength() for speed)
    centrality_vals <- eigen_centrality(g, weights = E(g)$weight)$vector
    most_central_gene <- names(which.max(centrality_vals))
    
    for (this_exposure_window in c("mat", "child")) {
      for (pfas in pfas_order) {
        
        # Fast data.table lookup
        key_dt <- dt_mediation[
          module == module_name &
            exposure_window == this_exposure_window &
            PFAS == pfas &
            significant_strict == TRUE &
            transcript_id %in% V(g)$name,
          .(transcript_id)
        ]
        
        sig_genes <- unique(key_dt$transcript_id)
        n_sig <- length(sig_genes)
        if (n_sig == 0) next
        
        # Distances
        dist_matrix <- distances(g, v = most_central_gene, to = sig_genes, weights = 1 / E(g)$weight)
        distances_vec <- as.numeric(dist_matrix)
        
        # Store results
        for (i in seq_len(n_sig)) {
          results[[res_idx]] <- list(
            module = module_name,
            exposure_label = ifelse(this_exposure_window == "child", "Fetal", "Maternal"),
            PFAS = pfas,
            transcript_id = sig_genes[i],
            distance = distances_vec[i]
          )
          res_idx <- res_idx + 1
        }
      }
    }
  }
  
  # Combine results and merge TPTE values
  result_dt <- rbindlist(results)
  result_dt <- merge(result_dt, tpte_values, by = "PFAS", all.x = TRUE)
  setorder(result_dt, exposure_label, PFAS)
  
  return(result_dt)
}

## ---- Run the function ----
transcript_distances <- compute_transcript_distances_dt(
  modules = modules,
  pfas_order = pfas_order,
  adjacency = adjacency,
  dt_mediation = dt_mediation,
  tpte_values = tpte_values
)

## ---- Save output ----
out_dir <- "/rsrch5/home/epi/stbresnahan/scratch/PFAS/transcript_distances_GA"

out_file <- file.path(
  out_dir,
  paste0("transcript_distances_", iter_val, ".rds")
)

saveRDS(transcript_distances, out_file)

message("Saved: ", out_file)
