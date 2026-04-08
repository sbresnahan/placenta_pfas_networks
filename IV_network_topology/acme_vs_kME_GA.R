########## Mediation for all transcripts for all PFAS

########## LOAD ALL REQUIRED LIBRARIES ##########

# Core tidyverse and data manipulation
library(tidyverse)
library(dplyr)

# WGCNA and gene expression analysis
library(WGCNA)

# Gene annotation and GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(rrvgo)
library(GOSemSim)

# Network analysis and visualization
library(igraph)
library(network)
library(sna)
library(ggraph)
library(tidygraph)
library(GGally)
library(ggnet)

# Plotting and visualization
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(RColorBrewer)
library(viridis)

# Parallel processing and progress bars
library(parallel)
library(pbapply)

########## END LIBRARY LOADING ##########

load("input_for_WGCNA.RData")
load("WGCNA_output.RData")
expr_PFAS0 <- expr_PFAS
expr_PFAS <- readRDS("expr_PFAS_all.rds")

get_module_genes <- function(module_id,color_ME_mapping,WGCNA.dat) {
  target_color <- color_ME_mapping$module_color[color_ME_mapping$ME_label == module_id]
  gene_idx <- which(WGCNA.dat$colors == target_color)
  return(rownames(expr_PFAS)[gene_idx])
}

calculate_module_eigengene <- function(module_genes, kme_threshold = 0.95) {
  # Get expression for module genes
  module_expr <- expr_PFAS[module_genes, , drop = FALSE]
  
  # Calculate initial eigengene
  if (nrow(module_expr) >= 3) {
    pca <- prcomp(t(module_expr), center = TRUE, scale. = FALSE)
    initial_me <- pca$x[, 1]
  } else {
    initial_me <- colMeans(module_expr)
  }
  
  # Calculate kME (correlation with eigengene)
  kme_values <- cor(t(module_expr), initial_me, use = "pairwise.complete.obs")
  
  # Filter genes by kME threshold
  keep_genes <- abs(kme_values) >= kme_threshold
  filtered_genes <- module_genes[keep_genes]
  
  if (length(filtered_genes) < 3) {
    # Use all genes if too few pass threshold
    filtered_genes <- module_genes
  }
  
  # Recalculate eigengene with filtered genes
  module_expr_filtered <- expr_PFAS[filtered_genes, , drop = FALSE]
  
  if (length(filtered_genes) >= 3) {
    pca_final <- prcomp(t(module_expr_filtered), center = TRUE, scale. = FALSE)
    final_me <- pca_final$x[, 1]
  } else {
    final_me <- colMeans(module_expr_filtered)
  }
  
  return(list(
    eigengene = final_me,
    n_genes_total = length(module_genes),
    n_genes_used = length(filtered_genes)
  ))
}

# Define all PFAS compounds
pfas_compounds_base <- c("PFBA", "PFBS", "PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUnDA")

# Define modules
module_names <- names(WGCNA.dat$MEs)

# Get all transcripts for all modules
all_module_transcripts_list <- list()

for (module_name in module_names) {
  module_genes <- get_module_genes(module_name,color_ME_mapping,WGCNA.dat)
  
  if (length(module_genes) == 0) next
  
  module_expr <- expr_PFAS[module_genes, ]
  
  # Calculate module eigengene
  me_result <- calculate_module_eigengene(module_genes, kme_threshold = 0)
  me_values <- me_result$eigengene
  
  # Calculate kME for each gene
  kme_values <- cor(t(module_expr), me_values, use = "pairwise.complete.obs")
  
  # Create data frame
  module_transcripts <- data.frame(
    transcript_id = rownames(kme_values),
    kME = as.numeric(kme_values),
    module = module_name
  ) %>%
    filter(abs(kME) > 0.0)
  
  all_module_transcripts_list[[length(all_module_transcripts_list) + 1]] <- module_transcripts
}

# Combine all module transcripts
module_transcripts_base <- bind_rows(all_module_transcripts_list)
module_transcripts_base1 <- readRDS("module_transcripts_base1.rds")

module_transcripts_base <- module_transcripts_base[
  na.omit(match(module_transcripts_base1$transcript_id, module_transcripts_base$transcript_id)),
]

# Get unique modules in each dataset
old_modules <- unique(module_transcripts_base$module)
new_modules <- unique(module_transcripts_base1$module)

results <- list()

# Loop through old and new modules to count transcript overlaps
for (old_mod in old_modules) {
  transcripts_old <- module_transcripts_base$transcript_id[module_transcripts_base$module == old_mod]
  
  for (new_mod in new_modules) {
    transcripts_new <- module_transcripts_base1$transcript_id[module_transcripts_base1$module == new_mod]
    n_overlap <- sum(transcripts_old %in% transcripts_new)
    
    if (n_overlap > 0) {
      results[[length(results) + 1]] <- data.frame(
        old = old_mod,
        new = new_mod,
        overlap = n_overlap
      )
    }
  }
}

# Combine overlap results into a dataframe
module_mapping_df <- do.call(rbind, results)

# Create a named vector for remapping
mapping_vec <- setNames(module_mapping_df$new, module_mapping_df$old)

# Apply new module names to module_transcripts_base
module_transcripts_base$module <- mapping_vec[module_transcripts_base$module]

color_ME_mapping$ME_label <- mapping_vec[color_ME_mapping$ME_label]


# Duplicate for ALL PFAS and PFAS Sources
all_module_transcripts_df <- expand_grid(
  module_transcripts_base,
  PFAS = pfas_compounds_base,
  exposure_window = c("mat", "child")
) %>%
  mutate(PFAS_column = paste0(PFAS, "_", exposure_window))

cat("Total rows in all_module_transcripts_df:", nrow(all_module_transcripts_df), "\n")
cat("These are already the full combinations to test\n")

# Create batch assignments (1000 per batch)
batch_size <- 1000
n_batches <- ceiling(nrow(all_module_transcripts_df) / batch_size)

all_module_transcripts_df <- all_module_transcripts_df %>%
  mutate(batch_id = rep(1:n_batches, each = batch_size, length.out = nrow(all_module_transcripts_df)))

cat("Total batches needed:", n_batches, "\n")


# UPDATE 3/21/26 - need to filter out non-spontaneous labors
stats.labor <- read.csv("labor_onset.csv",header=T)
covs_PFAS <- left_join(covs_PFAS,stats.labor)

covs_PFAS <- covs_PFAS[covs_PFAS$labour_onset=="1_Spontaneous",]
expr_PFAS <- expr_PFAS[,covs_PFAS$SampleID]

# Save inputs for mediation analysis
# saveRDS(all_module_transcripts_df, "all_module_transcripts_df.rds")
# saveRDS(expr_PFAS, "expr_PFAS_GA.rds")
# saveRDS(covs_PFAS, "covs_PFAS_GA.rds")
# saveRDS(module_transcripts_base, "module_transcripts_base.rds")


########## BIRTH WEIGHT MEDIATION ANALYSIS ##########

########### STEP 0: Performed on HPC

########### STEP 1: Merge all batch results from HPC

cat("=== Merging batch results from HPC ===\n")

# Get list of all result files
result_files <- list.files("WGCNA_mediation_results_GA", pattern = "mediation_batch_.*\\.rds", full.names = TRUE)
cat("Found", length(result_files), "batch result files\n")

# Read and combine all batches
all_batches <- lapply(result_files, function(file) {
  cat("Reading:", basename(file), "\n")
  readRDS(file)
})

# Combine into single dataframe
all_module_mediation_df <- bind_rows(all_batches)
cat("Total results after merging:", nrow(all_module_mediation_df), "\n")

# Load module_transcripts_base for kME values
module_transcripts_base <- readRDS("module_transcripts_base.rds")

# Add significance
all_module_mediation_df <- all_module_mediation_df %>%
  mutate(
    exposure_label = ifelse(exposure_window == "mat", "Maternal", "Fetal"),
    significant_strict = ACME_p < 0.05 & 
      ((ACME_CI_lower > 0 & ACME_CI_upper > 0) | 
         (ACME_CI_lower < 0 & ACME_CI_upper < 0))
  )

# Merge with kME values (signed, not absolute)
all_module_mediation_df <- all_module_mediation_df %>%
  left_join(module_transcripts_base %>% dplyr::select(transcript_id, module, kME), 
            by = c("transcript_id", "module"))

cat("Completed module transcript mediation:", nrow(all_module_mediation_df), "results\n")
cat("Significant (strict):", sum(all_module_mediation_df$significant_strict, na.rm = TRUE), "\n")

# Save merged results
# saveRDS(all_module_mediation_df, "all_pfas_mediation_results.rds")
save.image("acme_vs_kME_GA.RData")


# Update 1/12/26 - added bootstrapping across parameter grid, load in results and merge

# Directory containing bootstrap results
bootstrap_dir <- "/Users/stbresnahan/Desktop/PFAS/bootstrapped"

# Get all RData files
rdata_files <- list.files(bootstrap_dir, pattern = "^WGCNA_output_\\d+\\.RData$", full.names = TRUE)

# Initialize list to store all results
all_module_transcripts_list <- list()

# Loop through each iteration
for (rdata_file in rdata_files) {
  
  # Extract iteration number from filename
  current_iteration <- as.numeric(gsub(".*WGCNA_output_(\\d+)\\.RData", "\\1", rdata_file))
  
  cat("Processing iteration", current_iteration, "...\n")
  
  tryCatch({
    # Load the WGCNA results for this iteration
    load(rdata_file)
    
    # Get module eigengenes and module names
    MEs <- WGCNA.dat$MEs
    # samples <- dimnames(MEs)[[1]]
    module_names <- colnames(MEs)
    
    # Loop through each module
    for (module_name in module_names) {
      
      # Get genes in this module
      module_genes <- get_module_genes(module_name,color_ME_mapping,WGCNA.dat)
      
      if (length(module_genes) == 0) {
        cat("  Module", module_name, "has no genes, skipping...\n")
        next
      }
      
      # Get expression data for genes in this module
      module_expr <- expr_PFAS[module_genes, , drop = FALSE]
      
      # Get module eigengene values
      me_result <- calculate_module_eigengene(module_genes, kme_threshold = 0)
      me_values <- me_result$eigengene
      
      # Calculate kME for each gene
      kme_values <- cor(t(module_expr), me_values, use = "pairwise.complete.obs")
      
      # Create data frame
      module_transcripts <- data.frame(
        transcript_id = rownames(kme_values),
        kME = as.numeric(kme_values),
        module = module_name,
        module_size = length(module_genes),
        iteration = current_iteration,
        stringsAsFactors = FALSE
      )
      
      # Add to list
      all_module_transcripts_list[[length(all_module_transcripts_list) + 1]] <- module_transcripts
    }
    
    cat("  Completed iteration", current_iteration, "- processed", length(module_names), "modules\n")
    
  }, error = function(e) {
    cat("  ERROR in iteration", current_iteration, ":", conditionMessage(e), "\n")
  })
}

# Combine all module transcripts
cat("\nCombining results from all iterations...\n")
module_transcripts_all <- bind_rows(all_module_transcripts_list)

# Summary statistics
cat("\nSummary:\n")
cat("Total rows:", nrow(module_transcripts_all), "\n")
cat("Unique transcripts:", length(unique(module_transcripts_all$transcript_id)), "\n")
cat("Iterations processed:", length(unique(module_transcripts_all$iteration)), "\n")
cat("Total module assignments:", nrow(module_transcripts_all), "\n")

# Fix module names
# module_transcripts_all$module <- mapping_vec[module_transcripts_all$module]
# module_transcripts_all <- module_transcripts_all[!is.na(module_transcripts_all$module),]

# module_transcripts_all <- module_transcripts_all[!module_transcripts_all$module=="ME0",]
module_transcripts_all <- module_transcripts_all[module_transcripts_all$module_size>49,]

all_module_mediation_df <- all_module_mediation_df[all_module_mediation_df$transcript_id%in%module_transcripts_all$transcript_id,]

# Save results
# save(module_transcripts_all, file = "/Users/stbresnahan/Desktop/PFAS/module_transcripts_bootstrap.RData")
# write.csv(module_transcripts_all, 
#           file = "/Users/stbresnahan/Desktop/PFAS/module_transcripts_bootstrap.csv", 
#           row.names = FALSE)



save.image("acme_vs_kME_GA.RData")

########### SAVE FOR RESTARTING ANALYSES HERE ###################
########## Compare ACME between DTEs and hubs

# ---- Setup ----
results_dir <- "/Users/stbresnahan/Desktop/PFAS/results"
pfas_compounds <- c("PFBA", "PFBS", "PFDA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFUnDA")

# ---- Process DTEs ----
dte_summary_list <- list()

for (pfas in pfas_compounds) {
  for (exposure in c("mat", "child")) {
    dte_file <- file.path(results_dir, paste0(pfas, "_", exposure, "_DTE.tsv"))
    if (!file.exists(dte_file)) next
    
    dte_data <- read_tsv(dte_file, show_col_types = FALSE)
    
    dte_sig <- dte_data %>%
      filter(FDR < 0.1, abs(log2FC) > 1) %>%
      mutate(
        PFAS = pfas,
        exposure_window = if_else(exposure == "mat", "Maternal", "Child")
      ) %>%
      select(PFAS, exposure_window, Transcript)
    
    dte_summary_list[[length(dte_summary_list) + 1]] <- dte_sig
  }
}

dte_summary <- bind_rows(dte_summary_list)

# ---- Merge DTEs with ACME (for the same PFAS) ----
dte_with_mediation <- dte_summary %>%
  left_join(
    all_module_mediation_df %>%
      mutate(exposure_window = case_when(
        exposure_window == "mat" ~ "Maternal",
        exposure_window == "child" ~ "Child",
        TRUE ~ exposure_window
      )) %>%
      select(PFAS, exposure_window, transcript_id, ACME, ACME_p, ACME_CI_lower, ACME_CI_upper),
    by = c("PFAS", "exposure_window", "Transcript" = "transcript_id")
  ) %>%
  mutate(
    mediator_significant = ACME_p < 0.05 & (ACME_CI_lower * ACME_CI_upper > 0),
    category = "DTE"
  ) %>%
  select(PFAS, exposure_window, transcript_id = Transcript, ACME, mediator_significant, category)

dte_with_mediation <- dte_with_mediation[!is.na(dte_with_mediation$ACME),]

# ---- Identify top 10 hub per module*iteration ----
top_hubs <- module_transcripts_all %>%
  group_by(module, iteration) %>%
  slice_max(order_by = abs(kME), n = 10, with_ties = FALSE) %>%
  ungroup()

# Keep all PFAS × exposure_window rows for top hubs
hub_with_mediation <- all_module_mediation_df %>%
  filter(transcript_id %in% top_hubs$transcript_id) %>%
  mutate(
    exposure_window = case_when(
      exposure_window == "mat" ~ "Maternal",
      exposure_window == "child" ~ "Child",
      TRUE ~ exposure_window
    ),
    mediator_significant = ACME_p < 0.05 & (ACME_CI_lower * ACME_CI_upper > 0),
    category = ifelse(mediator_significant, "Hub (Mediator)", "Hub (Non-mediator)")
  ) %>%
  select(PFAS, exposure_window, transcript_id, ACME, mediator_significant, category)

# ---- Combine DTEs and Hubs ----
mediation_DTE_vs_Hub_ACME <- bind_rows(dte_with_mediation, hub_with_mediation)

# ---- Absolute ACME ----
mediation_DTE_vs_Hub_ACME <- mediation_DTE_vs_Hub_ACME %>%
  mutate(
    ACME_abs = abs(ACME),
    category = factor(category, levels = c("DTE", "Hub (Non-mediator)", "Hub (Mediator)"))
  )

set.seed(123)
n_bootstrap <- 100
n_sample <- 100
bootstrap_results <- data.frame(
  iteration = integer(),
  diff_NonMed_vs_DTE = numeric(),
  diff_Med_vs_DTE = numeric(),
  diff_Med_vs_NonMed = numeric(),
  pval_NonMed_vs_DTE = numeric(),
  pval_Med_vs_DTE = numeric(),
  pval_Med_vs_NonMed = numeric()
)
for (i in 1:n_bootstrap) {
  # Sample with replacement from each category
  bootstrap_sample <- mediation_DTE_vs_Hub_ACME %>%
    group_by(category) %>%
    slice_sample(n = n_sample, replace = TRUE) %>%
    ungroup()
  
  # ANOVA and Tukey HSD
  anova_test <- aov(ACME_abs ~ category, data = bootstrap_sample)
  posthoc <- TukeyHSD(anova_test)
  
  # Store differences and p-values
  bootstrap_results <- rbind(bootstrap_results, data.frame(
    iteration = i,
    diff_NonMed_vs_DTE = posthoc$category["Hub (Non-mediator)-DTE", "diff"],
    diff_Med_vs_DTE = posthoc$category["Hub (Mediator)-DTE", "diff"],
    diff_Med_vs_NonMed = posthoc$category["Hub (Mediator)-Hub (Non-mediator)", "diff"],
    pval_NonMed_vs_DTE = posthoc$category["Hub (Non-mediator)-DTE", "p adj"],
    pval_Med_vs_DTE = posthoc$category["Hub (Mediator)-DTE", "p adj"],
    pval_Med_vs_NonMed = posthoc$category["Hub (Mediator)-Hub (Non-mediator)", "p adj"]
  ))
}
# Calculate means across bootstrap iterations
mean_diff_NonMed_vs_DTE <- mean(bootstrap_results$diff_NonMed_vs_DTE)
mean_diff_Med_vs_DTE <- mean(bootstrap_results$diff_Med_vs_DTE)
mean_diff_Med_vs_NonMed <- mean(bootstrap_results$diff_Med_vs_NonMed)
p_val_NonMed_vs_DTE <- mean(bootstrap_results$pval_NonMed_vs_DTE)
p_val_Med_vs_DTE <- mean(bootstrap_results$pval_Med_vs_DTE)
p_val_Med_vs_NonMed <- mean(bootstrap_results$pval_Med_vs_NonMed)

# Create significance stars
sig_stars_1 <- case_when(
  p_val_NonMed_vs_DTE < 0.001 ~ "***",
  p_val_NonMed_vs_DTE < 0.01 ~ "**",
  p_val_NonMed_vs_DTE < 0.05 ~ "*",
  TRUE ~ "ns"
)

sig_stars_2 <- case_when(
  p_val_Med_vs_DTE < 0.001 ~ "***",
  p_val_Med_vs_DTE < 0.01 ~ "**",
  p_val_Med_vs_DTE < 0.05 ~ "*",
  TRUE ~ "ns"
)

sig_stars_3 <- case_when(
  p_val_Med_vs_NonMed < 0.001 ~ "***",
  p_val_Med_vs_NonMed < 0.01 ~ "**",
  p_val_Med_vs_NonMed < 0.05 ~ "*",
  TRUE ~ "ns"
)

# ---- Count number per group for axis labels ----
counts <- mediation_DTE_vs_Hub_ACME %>%
  dplyr::count(category) %>%
  dplyr::mutate(label = paste0(category, "\n(n=", n, ")"))

# ---- Plot ----
y_max <- max(mediation_DTE_vs_Hub_ACME$ACME_abs, na.rm = TRUE)

g1 <- ggplot(mediation_DTE_vs_Hub_ACME, aes(x = category, y = ACME_abs)) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  geom_jitter(
    data = filter(mediation_DTE_vs_Hub_ACME, mediator_significant == FALSE),
    aes(x = category, y = ACME_abs),
    position = position_jitter(width = 0.2),
    color = "grey50", alpha = 0.05, size = 2
  ) +
  geom_jitter(
    data = filter(mediation_DTE_vs_Hub_ACME, mediator_significant == TRUE),
    aes(x = category, y = ACME_abs),
    position = position_jitter(width = 0.2),
    color = "#C4663E", alpha = 0.2, size = 2
  ) +
  scale_x_discrete(labels = counts$label) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "|ACME|",
       title="Module hubs, rather than differentially expressed\ntranscripts (DETs), mediate PFAS-gestational age effects",
       subtitle="WGCNA parameter space: 144 iterations\nHubs: top 10 transcripts per module per iteration by |kME|\nMediators: p < 0.05, CIs exclude 0") +
  theme(legend.position = "none",
        plot.title=element_text(size=13,face="bold"),
        plot.subtitle=element_text(size=12)) +
  # Significance bar 1: DTE vs Hub (Non-mediator)
  annotate("segment", x = 1, xend = 2, y = y_max * 1.08, yend = y_max * 1.08,
           color = "black", size = 0.5) +
  annotate("text", x = 1.5, y = y_max * 1.12,
           label = sig_stars_1, fontface = "bold", size = 6) +
  # Significance bar 2: DTE vs Hub (Mediator)
  annotate("segment", x = 1, xend = 3, y = y_max * 1.24, yend = y_max * 1.24,
           color = "black", size = 0.5) +
  annotate("text", x = 2, y = y_max * 1.28,
           label = sig_stars_2, fontface = "bold", size = 6) +
  # Significance bar 3: Hub (Non-mediator) vs Hub (Mediator)
  annotate("segment", x = 2, xend = 3, y = y_max * 1.14, yend = y_max * 1.14,
           color = "black", size = 0.5) +
  annotate("text", x = 2.5, y = y_max * 1.18,
           label = sig_stars_3, fontface = "bold", size = 6)
g1

ggsave("dte_vs_hub_ACME_GA.pdf", g1, width = 6.75, height = 6, units = "in", dpi = 300)

save.image("acme_vs_kME_GA.RData")



########## Average Direct Effect (ADE) forest plot
tpte_values <- data.frame(
  PFAS = c("PFBA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFBS", "PFHxS", "PFOS"),
  mean_tpte = c(1.97, 1.17, 0.99, 0.43, 0.46, 2.04, 0.64, 0.36)
)

# Calculate average direct effect by PFAS and PFAS Source
ade_effects <- all_module_mediation_df %>%
  filter(!is.na(ADE), !is.na(ADE_CI_lower), !is.na(ADE_CI_upper), 
         significant_strict == TRUE) %>%
  mutate(
    ADE_SE = (ADE_CI_upper - ADE_CI_lower) / (2 * 1.96),
    weight = 1 / (ADE_SE^2)  # Inverse variance weight
  ) %>%
  group_by(PFAS, exposure_label) %>%
  summarise(
    # Weighted mean
    ADE = sum(ADE * weight) / sum(weight),
    # Pooled SE
    ADE_SE = sqrt(1 / sum(weight)),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    ADE_CI_lower = ADE - 1.96 * ADE_SE,
    ADE_CI_upper = ADE + 1.96 * ADE_SE
  ) %>%
  left_join(tpte_values, by = "PFAS") %>%
  arrange(mean_tpte) %>%
  mutate(PFAS = factor(PFAS, levels = unique(PFAS[order(mean_tpte)])))

p_ade_forest <- ggplot(ade_effects, aes(x = PFAS, y = ADE, color = exposure_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(
    aes(ymin = ADE_CI_lower, ymax = ADE_CI_upper),
    position = position_dodge(width = 0.6),
    width = 0.8,
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  labs(
    x = "PFAS",
    y = "Average Direct Effect (ADE) ± 95% CI",
    title = "Direct effects of PFAS on gestational age\nvary by exposure source and compound"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + coord_flip()
p_ade_forest

saveRDS(ade_effects,"Figures/ade_effects_GA.RDS")

ggsave("p_ade_forest_GA.pdf", p_ade_forest, width = 6.75, height = 6, units = "in", dpi = 300)





################# PLOT TRENDS ####################

########### STEP 1: Prepare histogram data

# Prepare histogram data
all_module_mediation_df0 <- all_module_mediation_df
all_module_mediation_df$module <- NULL
all_module_mediation_df$kME <- NULL

all_module_mediation_df <- left_join(all_module_mediation_df,module_transcripts_all,
                                     by="transcript_id",relationship = "many-to-many")

histogram_data <- all_module_mediation_df %>%
  mutate(
    abs_kME = abs(kME),
    significance_group = ifelse(significant_strict, "Significant", "Non-significant")
  )

cat("\n=== Summarizing peak |kME| per iteration ===\n")

combined_results <- histogram_data %>%
  group_by(iteration, exposure_label, PFAS, significance_group) %>%
  summarise(
    peak_x = if (sum(!is.na(abs_kME)) >= 2) {
      d <- density(abs_kME)
      d$x[which.max(d$y)]
    } else {
      NA_real_
    },
    median_x = median(abs_kME, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


########### STEP 2: Load TPTE and module eigengene correlations

# Get ME*PFAS cors
pfas_vars <- c("PFBA", "PFBS", "PFHxS", "PFOS", "PFOA", "PFNA", "PFDA", "PFUnDA")

pfas_map <- list(
  mat   = paste0(pfas_vars, "_mat"),
  child = paste0(pfas_vars, "_child")
)

me_pfas_cors_list <- list()

for (rdata_file in rdata_files) {
  
  iteration <- as.numeric(gsub(".*WGCNA_output_(\\d+)\\.RData", "\\1", rdata_file))
  cat("Processing iteration", iteration, "...\n")
  
  tryCatch({
    
    load(rdata_file)
    
    MEs <- WGCNA.dat$MEs
    MEs <- MEs[row.names(MEs)%in%covs_PFAS$SampleID,]
    
    covs_sub <- covs_PFAS
    
    for (window in names(pfas_map)) {
      
      for (pfas in pfas_vars) {
        
        pfas_var <- paste0(pfas, "_", window)
        
        if (!pfas_var %in% colnames(covs_sub)) next
        
        keep <- !is.na(covs_sub[[pfas_var]])
        
        ## require BOTH maternal and child PFAS
        if (window == "mat") {
          keep <- keep & !is.na(covs_sub[[paste0(pfas, "_child")]])
        } else {
          keep <- keep & !is.na(covs_sub[[paste0(pfas, "_mat")]])
        }
        
        if (sum(keep) < 10) next  # safety
        
        ME_sub   <- MEs[keep, , drop = FALSE]
        PFAS_sub <- covs_sub[[pfas_var]][keep]
        
        for (ME in colnames(ME_sub)) {
          
          ct <- suppressWarnings(
            cor.test(ME_sub[[ME]], PFAS_sub)
          )
          
          me_pfas_cors_list[[length(me_pfas_cors_list) + 1]] <-
            data.frame(
              PFAS = pfas,
              exposure_window = window,
              ME = ME,
              r = unname(ct$estimate),
              p_value = ct$p.value,
              ci_lower = ct$conf.int[1],
              ci_upper = ct$conf.int[2],
              n = sum(keep),
              iteration = iteration,
              stringsAsFactors = FALSE
            )
        }
      }
    }
    
    cat("  Completed iteration", iteration, "\n")
    
  }, error = function(e) {
    cat("  ERROR in iteration", iteration, ":", conditionMessage(e), "\n")
  })
}

me_pfas_cors_df <- do.call(rbind, me_pfas_cors_list)

# Select columns
me_pfas_cors_df <- me_pfas_cors_df %>%
  dplyr::select(PFAS, exposure_window, ME, r, p_value, n)

# Rename columns manually
names(me_pfas_cors_df)[3] <- "module"
names(me_pfas_cors_df)[4] <- "ME_PFAS_cor"
row.names(me_pfas_cors_df) <- NULL

# Calculate mean absolute correlation across modules for each PFAS-exposure combo
me_pfas_summary <- me_pfas_cors_df %>%  # Use me_pfas_cors_df, not me_pfas_cors
  mutate(exposure_label = ifelse(exposure_window == "mat", "Maternal", "Fetal")) %>%
  group_by(PFAS, exposure_label) %>%
  summarise(
    mean_ME_PFAS_cor = mean(abs(ME_PFAS_cor), na.rm = TRUE),
    median_ME_PFAS_cor = median(abs(ME_PFAS_cor), na.rm = TRUE),
    max_ME_PFAS_cor = max(abs(ME_PFAS_cor), na.rm = TRUE),
    n_modules = n(),
    n_significant = sum(p_value < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

cat("Summary statistics:\n")
print(me_pfas_summary)

########### STEP 3: Merge peak data with TPTE and ME correlations

# Extract peak data for significant mediators only
peak_data <- combined_results %>%
  filter(significance_group == "Significant") %>%
  dplyr::select(PFAS, exposure_label, peak_x, median_x, n, iteration)

# Merge with TPTE
peak_data_with_tpte <- peak_data %>%
  left_join(tpte_values, by = "PFAS")

# Merge with ME correlations
peak_data_full <- peak_data_with_tpte %>%
  left_join(me_pfas_summary, by = c("PFAS", "exposure_label"))

########### VISUALIZATION 1: Density plots for all PFAS

# Order PFAS by increasing TPTE for faceting
pfas_ordered_by_tpte <- tpte_values[order(tpte_values$mean_tpte),"PFAS"]

# Calculate peak of each density distribution
density_peaks <- histogram_data %>%
  group_by(PFAS, exposure_label, significance_group) %>%
  summarise(
    peak_x = density(abs_kME, na.rm = TRUE)$x[which.max(density(abs_kME, na.rm = TRUE)$y)],
    .groups = "drop"
  ) %>%
  mutate(
    PFAS_ordered = factor(PFAS, levels = pfas_ordered_by_tpte)
  )

histogram_data <- histogram_data %>%
  mutate(PFAS_ordered = factor(PFAS, levels = pfas_ordered_by_tpte))

plot_histogram_data <- all_module_mediation_df %>%
  group_by(transcript_id, PFAS, exposure_label, significant_strict) %>%
  summarise(
    abs_kME = mean(abs(kME), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    significance_group = ifelse(significant_strict, "Significant", "Non-significant"),
    PFAS_ordered = factor(PFAS, levels = pfas_ordered_by_tpte)  # order by TPTE
  ) %>%
  select(abs_kME, significance_group, PFAS_ordered, exposure_label)

p_kde_enhanced <- plot_histogram_data %>%
  ggplot(aes(x = abs_kME, fill = significance_group, color = significance_group)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(data = density_peaks %>% filter(significance_group == "Significant"), 
             aes(xintercept = peak_x, color = significance_group), 
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE, alpha = 1) +
  scale_fill_manual(values = c("Significant" = "#C4663E", "Non-significant" = "grey50")) +
  scale_color_manual(values = c("Significant" = "#C4663E", "Non-significant" = "grey50")) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), labels = c("0", "0.5", "1")) +
  facet_grid(exposure_label ~ PFAS_ordered, scales = "free_y") +
  labs(
    x = "Eigengene correlation (|kME|)",
    y = "Density",
    title = "Co-expression strength of mediators scales with\ntransplacental transfer efficiency (TPTE)",
    subtitle = "PFAS ordered by TPTE (cord:maternal serum ratio);\ndashed lines indicate co-expression strength (peak |kME|)",
    fill = "Mediation Status",
    color = "Mediation Status"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    strip.text = element_text(size = 7, face = "bold"),
    panel.spacing = unit(0.3, "lines"),
    axis.text = element_text(size = 7),
    panel.grid = element_line(color = "grey90", linewidth = 0.2)
  )

p_kde_enhanced

ggsave("p_kde_all_pfas_GA.pdf", p_kde_enhanced, width = 5.25, height = 4.5, units = "in", dpi = 300)


########### VISUALIZATION 2: TPTE vs Peak |kME|
cat("\n=== Creating TPTE vs peak |kME| plot ===\n")
library(dplyr)
library(ggplot2)
library(ggrepel)

# Summarize peak |kME| across iterations
peak_data_summary <- peak_data_full %>%
  group_by(PFAS, exposure_label, mean_tpte) %>%
  summarise(
    peak_x_mean = mean(peak_x, na.rm = TRUE),
    peak_x_sd   = sd(peak_x, na.rm = TRUE),
    n_iter     = n(),
    peak_x_se  = peak_x_sd / sqrt(n_iter),
    .groups = "drop"
  )

# Recalculate correlation stats using the averaged peak_x
tpte_stats_for_plot <- peak_data_summary %>%
  group_by(exposure_label) %>%
  summarise(
    r = cor(mean_tpte, peak_x_mean, use = "complete.obs"),
    r_squared = summary(lm(peak_x_mean ~ mean_tpte))$r.squared,
    p_value = summary(lm(peak_x_mean ~ mean_tpte))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value)
  )

# Determine label positions
label_y_start <- 1
label_y_spacing <- max(peak_data_summary$peak_x_mean, na.rm = TRUE) * 0.15

# Plot with error bars
p_tpte_vs_peak <- peak_data_summary %>%
  ggplot(aes(x = mean_tpte, y = peak_x_mean, color = exposure_label, fill = exposure_label)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_errorbar(aes(ymin = peak_x_mean - peak_x_se, ymax = peak_x_mean + peak_x_se),
                width = 0.08, alpha = 0.6, linewidth = 0.8) +
  geom_smooth(aes(group = exposure_label), method = "lm", se = TRUE, alpha = 0.1, linewidth = 1.2) +
  geom_text_repel(aes(label = PFAS), 
                  size = 3.5, 
                  show.legend = FALSE,
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0) +
  # Add correlation labels
  geom_label(data = tpte_stats_for_plot %>% filter(exposure_label == "Fetal"),
             aes(x = min(peak_data_summary$mean_tpte, na.rm = TRUE) + 0.1, 
                 y = label_y_start,
                 label = label),
             hjust = 0, vjust = 1, size = 3.5, 
             label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE,
             color = "#4682B4", fill = "white") +
  geom_label(data = tpte_stats_for_plot %>% filter(exposure_label == "Maternal"),
             aes(x = min(peak_data_summary$mean_tpte, na.rm = TRUE) + 0.1, 
                 y = label_y_start - label_y_spacing,
                 label = label),
             hjust = 0, vjust = 1, size = 3.5,
             label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE,
             color = "#F08080", fill = "white") +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  scale_fill_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  labs(
    x = "TPTE",
    y = "Co-expression strength (peak |kME|)",
    title = "Co-expression strength of fetal PFAS-gestational age mediator\ntranscripts scales with transplacental transfer efficiency",
    subtitle = "Points show mean peak |kME| across iterations; error bars = standard error"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_tpte_vs_peak)

ggsave("p_tpte_vs_peak_kme_GA.pdf", p_tpte_vs_peak, width = 6.75, height = 6, units = "in", dpi = 300)


GA.COEX <- peak_data_summary %>%
  group_by(exposure_label) %>%
  do({
    mod <- lm(peak_x_mean ~ mean_tpte, data = .)
    newdata <- data.frame(mean_tpte = seq(min(.$mean_tpte), max(.$mean_tpte), length.out = 100))
    preds <- predict(mod, newdata, interval = "confidence")
    cbind(newdata, as.data.frame(preds))
  }) %>%
  ungroup()

GA.COEX.SLOPES <- peak_data_summary %>%
  group_by(exposure_label) %>%
  do({
    mod <- lm(peak_x_mean ~ mean_tpte, data = .)
    
    coef_summary <- summary(mod)$coefficients
    slope <- coef_summary["mean_tpte", 1]
    se <- coef_summary["mean_tpte", 2]
    pval <- coef_summary["mean_tpte", 4]
    
    ci <- confint(mod)["mean_tpte", ]
    
    data.frame(
      slope = slope,
      slope_se = se,
      slope_lwr = ci[1],
      slope_upr = ci[2],
      slope_pval = pval,
      significant = pval < 0.05
    )
  }) %>%
  ungroup()
GA.COEX.SLOPES$outcome <- "Gestational age"
GA.COEX.SLOPES$metric <- "Co-expression strength"

saveRDS(GA.COEX,"ga_coex.RDS")
saveRDS(GA.COEX.SLOPES,"ga_coex_slopes.RDS")

########### VISUALIZATION: TPTE vs Average Absolute ACME of Significant Mediators
# Prepare data - significant mediators only, keep individual transcript values
all_module_mediation_df0$mean_tpte <- NULL

acme_tpte_sig <- all_module_mediation_df0 %>%
  filter(significant_strict == TRUE) %>%
  left_join(tpte_values, by = "PFAS") %>%
  filter(!is.na(mean_tpte), !is.na(ACME)) %>%
  mutate(
    abs_ACME = abs(ACME),
    acme_direction = if_else(ACME > 0, "Positive ACME", "Negative ACME")
  )

# Order PFAS by TPTE
pfas_order <- tpte_values %>% arrange(mean_tpte) %>% pull(PFAS)
acme_tpte_sig$PFAS <- factor(acme_tpte_sig$PFAS, levels = pfas_order)

# Overall ANOVA testing interaction between PFAS and PFAS Source
anova_model <- aov(abs_ACME ~ PFAS * exposure_label, data = acme_tpte_sig)
anova_summary <- summary(anova_model)
print(anova_summary)

# Perform t-tests for each PFAS and direction
pval_data <- acme_tpte_sig %>%
  group_by(PFAS, acme_direction) %>%
  summarise(
    p_value = t.test(abs_ACME[exposure_label == "Maternal"], 
                     abs_ACME[exposure_label == "Fetal"])$p.value,
    y_max = max(abs_ACME, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    y_bracket = y_max * 1.1,
    sig_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    x_num = as.numeric(PFAS),
    xmin = x_num - 0.2,
    xmax = x_num + 0.2
  ) %>%
  filter(sig_label != "ns")

# Count points per PFAS and exposure
count_data <- acme_tpte_sig %>%
  group_by(PFAS, acme_direction) %>%
  summarise(
    n_maternal = sum(exposure_label == "Maternal"),
    n_fetal = sum(exposure_label == "Fetal"),
    .groups = "drop"
  ) %>%
  mutate(
    x_num = as.numeric(PFAS),
    y_pos = if_else(acme_direction == "Negative ACME", 0.16, 0.23)
  ) %>%
  pivot_longer(cols = starts_with("n_"), names_to = "exposure_type", values_to = "n") %>%
  mutate(
    exposure_label = if_else(exposure_type == "n_maternal", "Maternal", "Fetal"),
    label = paste0("n=", n),
    y_offset = case_when(
      exposure_label == "Maternal" & acme_direction == "Positive ACME" ~ y_pos + 0.02,
      exposure_label == "Maternal" & acme_direction == "Negative ACME" ~ y_pos + 0.015,
      TRUE ~ y_pos
    )
  )

p_acme_boxplot <- ggplot(acme_tpte_sig, aes(x = PFAS, y = abs_ACME, color = exposure_label)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA, fill = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    alpha = 0.1, size = 1
  ) +
  geom_segment(data = pval_data,
               aes(x = xmin, xend = xmax, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_text(data = pval_data,
            aes(x = x_num, y = y_bracket * 1.03, label = sig_label),
            inherit.aes = FALSE, size = 5) +
  geom_text(data = count_data,
            aes(x = x_num, y = y_offset, label = label, color = exposure_label),
            vjust = 0.5, hjust = 0.5, size = 5, fontface = "bold") +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  facet_wrap(~acme_direction, nrow = 2, scales = "free_y") +
  labs(
    x = "PFAS (ordered by TPTE)",
    y = "|ACME| of Significant Mediators",
    title = "PFAS-gestational age mediation effects vary by exposure source and compound"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "right",
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold")
  )
print(p_acme_boxplot)

ggsave("p_tpte_vs_mean_acme_GA.pdf", p_acme_boxplot, width = 7.25, height = 7, units = "in", dpi = 300)


# Ridgeline plot version:
library(ggridges)

# Count points per PFAS and exposure
count_data <- acme_tpte_sig %>%
  group_by(PFAS) %>%
  summarise(
    n_maternal = sum(exposure_label == "Maternal"),
    n_fetal = sum(exposure_label == "Fetal"),
    .groups = "drop"
  )

p_acme_ridges <- ggplot(acme_tpte_sig, aes(x = ACME, y = PFAS, fill = exposure_label)) +
  geom_density_ridges(
    alpha = 0.6,
    scale = 0.9,
    quantile_lines = TRUE,
    quantiles = 2
  ) +
  # Fetal labels on left (x = -0.2)
  geom_text(data = count_data,
            aes(x = -0.24, y = PFAS, label = paste0("n = ", n_fetal)),
            hjust = 1, vjust = -0.5, size = 3, fontface = "bold",
            color = "#4682B4", inherit.aes = FALSE) +
  # Maternal labels on right (x = 0.2)
  geom_text(data = count_data,
            aes(x = 0.24, y = PFAS, label = paste0("n = ", n_maternal)),
            hjust = 0, vjust = -0.5, size = 3, fontface = "bold",
            color = "#F08080", inherit.aes = FALSE) +
  scale_fill_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  labs(
    x = "ACME of Significant Mediators",
    y = "PFAS",
    title = "PFAS-gestational age effects mediated by placental transcript\nexpression vary by exposure source and compound"
  ) + theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top")
print(p_acme_ridges)

ggsave("p_tpte_vs_mean_acme_ridges_GA.pdf", p_acme_ridges, width = 6.75, height = 6,, units = "in", dpi = 300)


# Scatterplot of significant mediators by PFAS and PFAS Source
# PFAS ordered by TPTE, colored by PFAS Source

# Count significant mediators by PFAS and exposure_window
mediator_counts <- all_module_mediation_df %>%
  filter(significant_strict == TRUE) %>%
  distinct(transcript_id, PFAS, exposure_window) %>%
  group_by(PFAS, exposure_window) %>%
  summarise(n_significant = n(), .groups = "drop")

# Merge with tpte values and order PFAS by tpte (ascending)
mediator_counts <- mediator_counts %>%
  left_join(tpte_values, by = "PFAS") %>%
  arrange(mean_tpte) %>%
  mutate(PFAS = factor(PFAS, levels = unique(PFAS)))

# Create cleaner labels for PFAS Source
mediator_counts <- mediator_counts %>%
  mutate(exposure_label = case_when(
    exposure_window == "mat" ~ "Maternal",
    exposure_window == "child" ~ "Fetal",
    TRUE ~ exposure_window
  ))

# Calculate correlations for each PFAS Source
cor_stats <- mediator_counts %>%
  group_by(exposure_label) %>%
  summarise(
    r = cor(mean_tpte, n_significant),
    r2 = cor(mean_tpte, n_significant)^2,
    p = cor.test(mean_tpte, n_significant)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(
    label = sprintf("r = %.2f, R² = %.2f, p = %.3f", r, r2, p),
    color = c("#4682B4", "#F08080")[match(exposure_label, c("Fetal", "Maternal"))]
  )

# Create positions for labels (upper left quadrant)
label_y_start <- max(mediator_counts$n_significant)
label_y_spacing <- max(mediator_counts$n_significant) * 0.12

# Dodge points slightly for visibility
mediator_counts <- mediator_counts %>%
  mutate(mean_tpte_dodged = ifelse(exposure_window == "mat", 
                                   mean_tpte - 0.02, 
                                   mean_tpte + 0.02))

p <- ggplot(mediator_counts, aes(x = mean_tpte_dodged, y = n_significant, color = exposure_label, fill = exposure_label)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_smooth(aes(x = mean_tpte), method = "lm", se = TRUE, alpha = 0.1, linewidth = 1) +
  geom_text_repel(aes(label = PFAS), 
                  size = 3.5, 
                  show.legend = FALSE,
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0) +
  scale_color_manual(values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
                     name = "PFAS Source") +
  scale_fill_manual(values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
                    name = "PFAS Source") +
  # Add correlation labels with borders
  geom_label(data = cor_stats %>% filter(exposure_label == "Fetal"),
             aes(x = min(mediator_counts$mean_tpte) + 0.1, 
                 y = label_y_start,
                 label = label, color = exposure_label),
             hjust = 0, vjust = 1, size = 3.5, 
             label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE,
             color = "#4682B4", fill = "white") +
  geom_label(data = cor_stats %>% filter(exposure_label == "Maternal"),
             aes(x = min(mediator_counts$mean_tpte) + 0.1, 
                 y = label_y_start - label_y_spacing,
                 label = label, color = exposure_label),
             hjust = 0, vjust = 1, size = 3.5,
             label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE,
             color = "#F08080", fill = "white") +
  # scale_y_continuous(limits=c(-500,1800),breaks=c(0,500,1000,1500)) +
  labs(x = "TPTE",
       y = "Count of Significant Mediators",
       title = "Transplacental transfer efficiency predicts transcriptional\nmediation of fetal PFAS-gestational age effects") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

# Display the plot
print(p)

ggsave("mediator_scatterplot_GA.pdf", p, width = 6.75, height = 6, units="in", dpi = 300)


GA.MC <- mediator_counts %>%
  group_by(exposure_label) %>%
  do({
    mod <- lm(n_significant ~ mean_tpte, data = .)
    newdata <- data.frame(mean_tpte = seq(min(.$mean_tpte), max(.$mean_tpte), length.out = 100))
    preds <- predict(mod, newdata, interval = "confidence")
    cbind(newdata, as.data.frame(preds))
  }) %>%
  ungroup()

GA.MC.SLOPES <- mediator_counts %>%
  group_by(exposure_label) %>%
  do({
    mod <- lm(n_significant ~ mean_tpte, data = .)
    
    coef_summary <- summary(mod)$coefficients
    slope <- coef_summary["mean_tpte", 1]
    se <- coef_summary["mean_tpte", 2]
    pval <- coef_summary["mean_tpte", 4]
    
    ci <- confint(mod)["mean_tpte", ]
    
    data.frame(
      slope = slope,
      slope_se = se,
      slope_lwr = ci[1],
      slope_upr = ci[2],
      slope_pval = pval,
      significant = pval < 0.05
    )
  }) %>%
  ungroup()
GA.MC.SLOPES$outcome <- "Gestational age"
GA.MC.SLOPES$metric <- "Mediator count"

saveRDS(GA.MC,"ga_mc.RDS")
saveRDS(GA.MC.SLOPES,"ga_mc_slopes.RDS")



########### VISUALIZATION: Boxplots of Enrichment ORs of Significant Mediators by Module


# Calculate proportion of significant mediators and Fisher's exact test for each module-PFAS-exposure combination
module_enrichment <- all_module_mediation_df %>%
  group_by(module, PFAS, exposure_window, exposure_label) %>%
  summarise(
    n_transcripts = n(),
    n_significant = sum(significant_strict, na.rm = TRUE),
    prop_significant = n_significant / n_transcripts,
    .groups = "drop"
  )

# Calculate overall proportion significant for Fisher's test background
overall_stats <- all_module_mediation_df %>%
  group_by(PFAS, exposure_window) %>%
  summarise(
    total_transcripts = n(),
    total_significant = sum(significant_strict, na.rm = TRUE),
    .groups = "drop"
  )

# Join and calculate Fisher's exact test for each module
module_enrichment_with_fisher <- module_enrichment %>%
  left_join(overall_stats, by = c("PFAS", "exposure_window")) %>%
  mutate(
    # Fisher's exact test for enrichment
    # Contingency table:
    # [sig in module, not sig in module]
    # [sig outside module, not sig outside module]
    fisher_result = purrr::pmap(
      list(n_significant, n_transcripts, total_significant, total_transcripts),
      function(sig_module, total_module, sig_total, total_all) {
        sig_outside <- sig_total - sig_module
        not_sig_module <- total_module - sig_module
        not_sig_outside <- (total_all - total_module) - sig_outside
        
        contingency <- matrix(
          c(sig_module, not_sig_module,
            sig_outside, not_sig_outside),
          nrow = 2,
          byrow = TRUE
        )
        test <- fisher.test(contingency, alternative = "greater")
        list(p_value = test$p.value, odds_ratio = test$estimate)
      }
    )
  ) %>%
  unnest_wider(fisher_result) %>%
  mutate(
    # Adjust p-values for multiple testing within each PFAS Source
    fisher_fdr = p.adjust(p_value, method = "BH")
  ) %>%
  mutate(
    # Create significance labels
    sig_label = case_when(
      fisher_fdr < 0.001 ~ "***",
      fisher_fdr < 0.01 ~ "**",
      fisher_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Order PFAS by TPTE
module_enrichment_with_fisher <- module_enrichment_with_fisher %>%
  left_join(tpte_values, by = "PFAS") %>%
  arrange(mean_tpte) %>%
  mutate(PFAS = factor(PFAS, levels = unique(PFAS[order(mean_tpte)])))

glm_df <- module_enrichment_with_fisher %>%
  filter(
    is.finite(odds_ratio),
    odds_ratio > 0,
    is.finite(mean_tpte),
    exposure_label %in% c("Maternal", "Fetal")
  ) %>%
  mutate(
    log_odds_ratio = log(odds_ratio),
    exposure_label = factor(exposure_label, levels = c("Fetal", "Maternal"))
  )

glm_tp <- glm(
  log_odds_ratio ~ mean_tpte * exposure_label,
  data = glm_df,
  family = gaussian()
)

summary(glm_tp)

p_module_boxplot <- module_enrichment_with_fisher %>%
  filter(
    is.finite(odds_ratio),
    odds_ratio > 0
  ) %>%
  ggplot(aes(x = PFAS, y = odds_ratio, color = exposure_label)) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    outlier.shape = NA,
    linewidth = 0.4
  ) +
  geom_jitter(
    aes(color = exposure_label),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.75
    ),
    alpha = 0.2,
    size = 1,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "gray50"
  ) +
  scale_y_log10() +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  labs(
    x = "PFAS (ordered by TPTE)",
    y = "Enrichment Odds Ratio (log10)",
    title = "Co-expression network enrichment for PFAS-gestational age mediators\nvaries by exposure source and compound",
    subtitle = "Points indicate fisher's exact test ORs per iteration x module"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Display
p_module_boxplot

ggsave("p_module_boxplot_GA.pdf", p_module_boxplot, width = 7.25, height = 6.75, units="in", dpi = 300)



########### VISUALIZATION: Average Proportion Significant vs Peak |kME|
cat("\n=== Creating scatterplot of average proportion significant vs peak |kME| ===\n")

library(ggplot2)
library(dplyr)
library(ggrepel)

# Calculate proportion of significant mediators by module-PFAS-exposure
prop_by_module <- all_module_mediation_df %>%
  group_by(module, PFAS, exposure_window, exposure_label,iteration) %>%
  summarise(
    n_transcripts = n(),
    n_significant = sum(significant_strict, na.rm = TRUE),
    prop_significant = n_significant / n_transcripts,
    .groups = "drop"
  )

# Calculate average proportion across modules for each PFAS-exposure combination
avg_prop_data <- prop_by_module %>%
  group_by(PFAS, exposure_window, exposure_label, iteration) %>%
  summarise(
    mean_prop_significant = mean(prop_significant, na.rm = TRUE),
    median_prop_significant = median(prop_significant, na.rm = TRUE),
    sd_prop_significant = sd(prop_significant, na.rm = TRUE),
    se_prop_significant = sd(prop_significant, na.rm = TRUE) / sqrt(n()),
    n_modules = n(),
    .groups = "drop"
  )

# Get peak |kME| data from peak_data_full
# peak_x in peak_data_full is the peak |kME| for each PFAS-exposure combination
peak_kme_data <- peak_data_full %>%
  dplyr::select(PFAS, exposure_label, peak_kME = peak_x, n_significant = n, iteration)

# Merge the data for plotting
plot_data <- avg_prop_data %>%
  left_join(peak_kme_data, by = c("PFAS", "exposure_label", "iteration")) %>%
  filter(!is.na(peak_kME))

plot_data_summary <- plot_data %>%
  group_by(PFAS, exposure_label) %>%
  summarise(
    peak_kME_mean = mean(peak_kME, na.rm = TRUE),
    peak_kME_sd   = sd(peak_kME, na.rm = TRUE),
    n_iter        = n(),
    peak_kME_se   = peak_kME_sd / sqrt(n_iter),
    mean_prop_significant_mean = mean(mean_prop_significant, na.rm = TRUE),
    mean_prop_significant_sd   = sd(mean_prop_significant, na.rm = TRUE),
    mean_prop_significant_se   = mean_prop_significant_sd / sqrt(n_iter),
    .groups = "drop"
  )

cor_stats_summary <- plot_data_summary %>%
  group_by(exposure_label) %>%
  summarise(
    r = cor(peak_kME_mean, mean_prop_significant_mean, use = "complete.obs"),
    r_squared = summary(lm(mean_prop_significant_mean ~ peak_kME_mean))$r.squared,
    p_value = summary(lm(mean_prop_significant_mean ~ peak_kME_mean))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value)
  )

# Determine label positions
label_x_start <- 0.1
label_y_start <- 0.01
label_y_spacing <- label_y_start * 0.15

# Plot
library(scales)

p_avg_prop_vs_peak_kme <- plot_data_summary %>%
  ggplot(aes(x = peak_kME_mean, y = mean_prop_significant_mean, color = exposure_label, fill = exposure_label)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_prop_significant_mean - mean_prop_significant_se,
                    ymax = mean_prop_significant_mean + mean_prop_significant_se),
                width = 0.02, alpha = 0.6) +
  geom_errorbarh(aes(xmin = peak_kME_mean - peak_kME_se,
                     xmax = peak_kME_mean + peak_kME_se),
                 alpha = 0.6, height = 0.001) +
  geom_smooth(aes(group = exposure_label), method = "lm", se = TRUE, alpha = 0.1, linewidth = 1.2) +
  geom_text_repel(aes(label = PFAS), 
                  size = 3.5, 
                  show.legend = FALSE,
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0) +
  # Add correlation labels
  geom_label(data = cor_stats_summary %>% filter(exposure_label == "Fetal"),
             aes(x = label_x_start, 
                 y = label_y_start,
                 label = label),
             hjust = 0, vjust = 1, size = 3.5, 
             label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE,
             color = "#4682B4", fill = "white") +
  geom_label(data = cor_stats_summary %>% filter(exposure_label == "Maternal"),
             aes(x = label_x_start, 
                 y = label_y_start - label_y_spacing,
                 label = label),
             hjust = 0, vjust = 1, size = 3.5,
             label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE,
             color = "#F08080", fill = "white") +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  scale_fill_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Mean Peak |kME| of Significant Mediators ± SE",
    y = "Mean Proportion Significant ± SE",
    title = "Average module enrichment for fetal PFAS-gestational age\nmediator transcripts scales with peak module membership"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_avg_prop_vs_peak_kme)

ggsave("p_avg_prop_sig_vs_peak_kme_GA.pdf", p_avg_prop_vs_peak_kme, 
       width = 6.75, height = 6.25, units = "in", dpi = 300)


save.image("acme_vs_kME_GA.RData")


