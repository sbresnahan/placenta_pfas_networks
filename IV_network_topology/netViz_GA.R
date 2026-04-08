# Visualize specific module networks


########## LOAD ALL REQUIRED LIBRARIES ##########

# Core tidyverse and data manipulation
library(tidyverse)
library(dplyr)

# WGCNA and gene expression analysis
library(WGCNA)
library(DESeq2)
library(RUVSeq)

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


load("~/Desktop/PFAS/acme_vs_kME_GA.RData")

expr_PFAS <- readRDS("expr_PFAS_all.rds")

# UPDATE 03/26/21 - filter out non-spontaneous labors
expr_PFAS <- expr_PFAS[,covs_PFAS$SampleID]

# adjacency <- adjacency(datExpr = t(expr_PFAS), power = 3, type = "signed") # Run on HPC
# saveRDS(adjacency,"adjacency_GA.rds")

adjacency <- readRDS("adjacency_GA.rds")

# Remove everything else
keep_objects <- c("all_module_mediation_df", "tpte_values", "pfas_order", 
                  "peak_data_full", "prop_by_module",
                  "expr_PFAS","adjacency")
rm(list = setdiff(ls(), keep_objects))

save.image("netViz_GA.RData")

################ Formal analysis of proximity to central node #################
library(tidyverse)
library(igraph)
library(data.table)

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


##### Run on HPC
# save(list=c("all_module_mediation_df","tpte_values",
#             "pfas_order","compute_transcript_distances_dt"),
#      file="data_for_distance_analysis.RDS")

##### Merge results

files <- list.files(
  "transcript_distances_GA",
  pattern = "^transcript_distances_[0-9]+\\.rds$",
  full.names = TRUE
)

transcript_distances <- files %>%
  lapply(function(f) {
    readRDS(f) %>%
      mutate(
        iteration = as.integer(
          str_extract(basename(f), "[0-9]+")
        )
      )
  }) %>%
  bind_rows()

transcript_distances <- transcript_distances %>%
  inner_join(
    all_module_mediation_df %>%
      select(
        PFAS,
        module,
        exposure_label,
        transcript_id,
        iteration,
        significant_strict
      ),
    by = c(
      "PFAS",
      "module",
      "exposure_label",
      "transcript_id",
      "iteration"
    )
  ) %>%
  filter(significant_strict)

# Order PFAS by increasing TPTE
pfas_ordered_by_tpte <- tpte_values %>% arrange(mean_tpte) %>% pull(PFAS)

# Create histogram data
histogram_data_distance <- transcript_distances %>%
  # Average abs_kME across iterations for each transcript-PFAS-exposure-significance
  group_by(transcript_id, PFAS, exposure_label) %>%
  summarise(
     distance = mean(distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Set PFAS factor levels for consistent plotting
  mutate(PFAS_ordered = factor(PFAS, levels = pfas_ordered_by_tpte))

# Create density plot
p_kde_distance <- histogram_data_distance %>%
  ggplot(aes(x = distance, fill = exposure_label, color = exposure_label)) +
  geom_density(linewidth = 0.8, alpha = 0.3) +
  scale_fill_manual(values = c("Maternal" = "#F08080", "Fetal" = "#4682B4")) +
  scale_color_manual(values = c("Maternal" = "#F08080", "Fetal" = "#4682B4")) +
  facet_wrap(~ PFAS_ordered, scales = "free_y", ncol = 4) +
  labs(
    x = "Distance from Central Hub",
    y = "Density",
    title = "Distribution of Network Distance for Significant Mediators by PFAS",
    subtitle = "PFAS ordered by increasing transplacental transfer efficiency (TPTE)",
    fill = "PFAS Source",
    color = "PFAS Source"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    axis.text = element_text(size = 8)
  )

print(p_kde_distance)
ggsave("p_kde_distance_mediators_GA.pdf", p_kde_distance, width = 12, height = 8, units = "in", dpi = 300)



########### VISUALIZATION: TPTE vs Peak Distance (closest to hub)
cat("\n=== Creating TPTE vs minimum distance plot ===\n")

# Function to get density peak
get_density_peak <- function(x) {
  if(length(x) < 2) return(NA_real_)
  d <- density(x, na.rm = TRUE)
  d$x[which.max(d$y)]
}

peak_distance_by_iter <- transcript_distances %>%
  group_by(iteration, PFAS, exposure_label, mean_tpte) %>%
  summarise(
    peak_distance = get_density_peak(distance),
    n_transcripts = n(),
    .groups = "drop"
  ) %>%
  filter(!is.na(peak_distance))

peak_distance_summary <- peak_distance_by_iter %>%
  group_by(PFAS, exposure_label, mean_tpte) %>%
  summarise(
    mean_peak_distance = mean(peak_distance),
    sd_peak_distance   = sd(peak_distance),
    se_peak_distance   = sd_peak_distance / sqrt(n()),
    n_iterations       = n(),
    .groups = "drop"
  )

tpte_stats_distance <- peak_distance_summary %>%
  group_by(exposure_label) %>%
  summarise(
    r = cor(mean_tpte, mean_peak_distance, use = "complete.obs"),
    r_squared = summary(lm(mean_peak_distance ~ mean_tpte))$r.squared,
    p_value = summary(lm(mean_peak_distance ~ mean_tpte))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value)
  )

label_y_start <- max(peak_distance_summary$mean_peak_distance, na.rm = TRUE)
label_y_spacing <- diff(range(peak_distance_summary$mean_peak_distance, na.rm = TRUE)) * 0.15

p_tpte_vs_peak_distance <- peak_distance_summary %>%
  ggplot(aes(
    x = mean_tpte,
    y = mean_peak_distance,
    color = exposure_label,
    fill  = exposure_label
  )) +
  geom_errorbar(
    aes(
      ymin = mean_peak_distance - se_peak_distance,
      ymax = mean_peak_distance + se_peak_distance
    ),
    width = 0.08,
    alpha = 0.6
  ) +
  geom_point(size = 4, alpha = 0.85) +
  geom_smooth(
    aes(group = exposure_label),
    method = "lm",
    se = TRUE,
    alpha = 0.1,
    linewidth = 1.2
  ) +
  geom_text_repel(
    aes(label = PFAS),
    size = 3.5,
    show.legend = FALSE,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3
  ) +
  geom_label(
    data = tpte_stats_distance %>% filter(exposure_label == "Fetal"),
    aes(x = 1.4, y = label_y_start, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3.5,
    color = "#4682B4",
    fill = "white"
  ) +
  geom_label(
    data = tpte_stats_distance %>% filter(exposure_label == "Maternal"),
    aes(x = 1.4, y = label_y_start - label_y_spacing, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3.5,
    color = "#F08080",
    fill = "white"
  ) +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  scale_fill_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  labs(
    x = "TPTE (Transplacental Transfer Efficiency)",
    y = "Peak distance from central hub",
    title = "Network centrality of PFAS–gestational age mediator transcripts\ndoes not scale with transplacental transfer efficiency",
    subtitle = "Points show mean ± SE across iterations"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_tpte_vs_peak_distance)

ggsave("p_tpte_vs_peak_distance_GA.pdf", p_tpte_vs_peak_distance, width = 6.75, height = 6, units = "in", dpi = 300)


GA.NC <- peak_distance_summary %>%
  group_by(exposure_label) %>%
  do({
    mod <- lm(mean_peak_distance ~ mean_tpte, data = .)
    newdata <- data.frame(mean_tpte = seq(min(.$mean_tpte), max(.$mean_tpte), length.out = 100))
    preds <- predict(mod, newdata, interval = "confidence")
    cbind(newdata, as.data.frame(preds))
  }) %>%
  ungroup()

GA.NC.SLOPES <- peak_distance_summary %>%
  group_by(exposure_label) %>%
  do({
    mod <- lm(mean_peak_distance ~ mean_tpte, data = .)
    
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
GA.NC.SLOPES$outcome <- "Gestational age"
GA.NC.SLOPES$metric <- "Centrality"

saveRDS(GA.NC,"ga_nc.RDS")
saveRDS(GA.NC.SLOPES,"ga_nc_slopes.RDS")



########### VISUALIZATION: Peak |kME| vs Peak Distance
cat("\n=== Creating scatterplot of peak |kME| vs peak distance ===\n")

# Merge peak |kME| with peak distance (iteration-level)
peak_kme_data <- peak_distance_by_iter %>%
  left_join(
    peak_data_full %>% select(PFAS, exposure_label, peak_x, iteration),
    by = c("PFAS", "exposure_label", "iteration")
  ) %>%
  rename(peak_kME = peak_x)

# Summarize across iterations for mean ± SE
plot_data_kme_summary <- peak_kme_data %>%
  group_by(PFAS, exposure_label) %>%
  summarise(
    mean_peak_kME = mean(peak_kME),
    se_peak_kME   = sd(peak_kME)/sqrt(n()),
    mean_peak_distance = mean(peak_distance),
    se_peak_distance   = sd(peak_distance)/sqrt(n()),
    n_iterations       = n(),
    .groups = "drop"
  )

cor_stats_kme_distance <- plot_data_kme_summary %>%
  group_by(exposure_label) %>%
  summarise(
    r = cor(mean_peak_kME, mean_peak_distance),
    r_squared = summary(lm(mean_peak_distance ~ mean_peak_kME))$r.squared,
    p_value = summary(lm(mean_peak_distance ~ mean_peak_kME))$coefficients[2, 4],
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value)
  )

label_x_start <- max(plot_data_kme_summary$mean_peak_kME, na.rm = TRUE) * 0.95
label_y_start <- max(plot_data_kme_summary$mean_peak_distance, na.rm = TRUE)
label_y_spacing <- diff(range(plot_data_kme_summary$mean_peak_distance, na.rm = TRUE)) * 0.15

p_peak_kme_vs_peak_distance <- ggplot(
  plot_data_kme_summary,
  aes(x = mean_peak_kME, y = mean_peak_distance, color = exposure_label, fill = exposure_label)
) +
  geom_errorbar(
    aes(
      ymin = mean_peak_distance - se_peak_distance,
      ymax = mean_peak_distance + se_peak_distance
    ),
    width = 0.03, alpha = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = mean_peak_kME - se_peak_kME,
      xmax = mean_peak_kME + se_peak_kME
    ),
    height = 0.15, alpha = 0.6
  ) +
  geom_point(size = 4, alpha = 0.85) +
  geom_smooth(aes(group = exposure_label), method = "lm", se = TRUE, alpha = 0.1, linewidth = 1.2) +
  geom_text_repel(
    aes(label = PFAS),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  geom_label(
    data = cor_stats_kme_distance %>% filter(exposure_label == "Fetal"),
    aes(x = label_x_start, y = label_y_start, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    size = 3.5,
    color = "#4682B4",
    fill = "white"
  ) +
  geom_label(
    data = cor_stats_kme_distance %>% filter(exposure_label == "Maternal"),
    aes(x = label_x_start, y = label_y_start - label_y_spacing, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    size = 3.5,
    color = "#F08080",
    fill = "white"
  ) +
  scale_color_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  scale_fill_manual(
    values = c("Maternal" = "#F08080", "Fetal" = "#4682B4"),
    name = "PFAS Source"
  ) +
  labs(
    x = "Peak |kME|",
    y = "Peak distance from central hub",
    title = "Network centrality of fetal PFAS-gestational age\nmediators scales with module membership",
    subtitle = "Points show mean ± SE across iterations"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_peak_kme_vs_peak_distance)

ggsave("p_peak_kme_vs_peak_distance_GA.pdf", p_peak_kme_vs_peak_distance, 
       width = 6.75, height = 6, units = "in", dpi = 300)



######### Diff D vs tpte

library(igraph)
library(dplyr)
library(tidyr)
library(ggraph)
library(tibble)

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

###### Run on HPC

##### Merge results

files <- list.files(
  "peak_distances_GA",
  pattern = "^transcript_peak_distances_[0-9]+\\.rds$",
  full.names = TRUE
)

# Load all iteration-level files
transcript_peak_distances <- files %>%
  lapply(function(f) {
    readRDS(f) %>%
      mutate(
        iteration = as.integer(str_extract(basename(f), "[0-9]+"))
      )
  }) %>%
  bind_rows()

# Compute |Fetal - Maternal| distances per module × PFAS × iteration
iteration_abs_diff <- transcript_peak_distances %>%
  select(module, PFAS, exposure_label, peak_distance, mean_tpte, iteration) %>%
  pivot_wider(
    names_from = exposure_label,
    values_from = peak_distance
  ) %>%
  mutate(abs_diff = abs(Fetal - Maternal)) %>%
  filter(is.finite(abs_diff))

# Filter valid modules (modules with ≥2 significant mediators in both exposures)
valid_modules <- prop_by_module %>%
  filter(n_significant >= 2) %>%
  group_by(PFAS, module) %>%
  summarise(n_exposures = n_distinct(exposure_label), .groups = "drop") %>%
  filter(n_exposures == 2) %>%
  select(PFAS, module)

iteration_abs_diff <- iteration_abs_diff %>%
  inner_join(valid_modules, by = c("PFAS", "module")) %>%
  mutate(PFAS = factor(PFAS, levels = tpte_values$PFAS[order(tpte_values$mean_tpte)]))

# Summarise across iterations: mean ± SE per PFAS
pfas_peaks <- iteration_abs_diff %>%
  group_by(PFAS, mean_tpte) %>%
  summarise(
    peak_abs_diff = mean(abs_diff, na.rm = TRUE),
    se = sd(abs_diff, na.rm = TRUE) / sqrt(n()),
    n_iterations = n(),
    .groups = "drop"
  )

# Correlation / regression statistics
reg_stats <- lm(peak_abs_diff ~ mean_tpte, data = pfas_peaks)
cor_label <- sprintf(
  "r = %.3f, R² = %.3f, p = %.3f",
  cor(pfas_peaks$mean_tpte, pfas_peaks$peak_abs_diff),
  summary(reg_stats)$r.squared,
  summary(reg_stats)$coefficients[2, 4]
)

# Determine label position
label_x <- min(pfas_peaks$mean_tpte) + 0.1
label_y <- 2.25

# Plot: mean ± SE per PFAS with regression line
gdtpte <- ggplot(pfas_peaks, aes(x = mean_tpte, y = peak_abs_diff)) +
  geom_point(size = 3, color = "black") +
  geom_errorbar(aes(ymin = peak_abs_diff - se, ymax = peak_abs_diff + se), width = 0.05) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_label_repel(aes(label = PFAS),
                   box.padding = 0.5,
                   point.padding = 0.3,
                   segment.color = "grey50",
                   segment.size = 0.3,
                   min.segment.length = 0,
                   size = 3.5,
                   fill = "white",
                   color = "black") +
  annotate("label",
           x = label_x, y = label_y,
           label = cor_label,
           hjust = 0, vjust = 1,
           size = 3.5,
           fill = "white", color = "black",
           label.padding = unit(0.3, "lines")) +
  theme_minimal(base_size = 12) +
  theme(plot.title=element_text(face="bold")) +
  labs(
    x = "Mean TPTE",
    y = expression(paste("Mean |", D[f], " - ", D[m], "| ± SE")),
    title = "Distance between PFAS-gestational age\nmediator hubs does not scale with TPTE"
  )

gdtpte

ggsave("gdtpte_GA.pdf", gdtpte, 
       width = 6.75, height = 6, units = "in", dpi = 300)


GA.COMP.SLOPES <- {
  mod <- lm(peak_abs_diff ~ mean_tpte, data = pfas_peaks)
  
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
}
GA.COMP.SLOPES$outcome <- "Gestational age"
GA.COMP.SLOPES$metric <- "Compartmentalization"

saveRDS(GA.COMP,"ga_comp.RDS")
saveRDS(GA.COMP.SLOPES,"ga_comp_slopes.RDS")





########## Peak distance vs ACME
library(dplyr)

# --- Helper function for safe SE calculation ---
safe_se <- function(x) {
  if(length(x) <= 1) return(0)
  sd(x, na.rm = TRUE) / sqrt(length(x))
}

# --- Step 1: Summarise ACME per iteration per PFAS × exposure × module ---
acme_summary_iter <- all_module_mediation_df %>%
  filter(significant_strict) %>%
  group_by(PFAS, exposure_label, iteration, module) %>%
  summarise(
    mean_abs_acme = mean(abs(ACME), na.rm = TRUE),
    n_sig = n(),
    .groups = "drop"
  )

# Then compute mean ± SE across iterations per PFAS × exposure
acme_summary <- acme_summary_iter %>%
  group_by(PFAS, exposure_label, iteration) %>%
  summarise(
    mean_abs_acme_iter = mean(mean_abs_acme, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(PFAS, exposure_label) %>%
  summarise(
    mean_abs_acme = mean(mean_abs_acme_iter, na.rm = TRUE),
    se_abs_acme   = safe_se(mean_abs_acme_iter),
    n_iterations  = n(),
    .groups = "drop"
  )

# --- Step 2: Summarise peak distance per iteration per PFAS × exposure × module ---
peak_summary_iter <- transcript_peak_distances %>%
  group_by(PFAS, exposure_label, iteration, module) %>%
  summarise(
    mean_peak_distance = mean(peak_distance, na.rm = TRUE),
    n_transcripts = n(),
    .groups = "drop"
  )

# Then compute mean ± SE across iterations per PFAS × exposure
peak_summary <- peak_summary_iter %>%
  group_by(PFAS, exposure_label, iteration) %>%
  summarise(
    mean_peak_distance_iter = mean(mean_peak_distance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(PFAS, exposure_label) %>%
  summarise(
    mean_peak_distance = mean(mean_peak_distance_iter, na.rm = TRUE),
    se_peak_distance   = safe_se(mean_peak_distance_iter),
    n_iterations       = n(),
    .groups = "drop"
  )

# --- Step 3: Combine ACME and peak distance summaries ---
pfas_summary <- peak_summary %>%
  left_join(acme_summary, by = c("PFAS", "exposure_label"))

# --- Step 4: Compute regression / correlation statistics ---
pfas_stats <- pfas_summary %>%
  group_by(exposure_label) %>%
  summarise(
    r         = cor(mean_abs_acme, mean_peak_distance, use = "complete.obs"),
    r_squared = summary(lm(mean_peak_distance ~ mean_abs_acme))$r.squared,
    p_value   = summary(lm(mean_peak_distance ~ mean_abs_acme))$coefficients[2, 4],
    .groups   = "drop"
  ) %>%
  mutate(label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value))

# --- Step 5: Dynamic label positioning ---
label_y_start <- 3.5
label_y_spacing <- (max(pfas_summary$mean_peak_distance, na.rm = TRUE) -
                      min(pfas_summary$mean_peak_distance, na.rm = TRUE)) * 0.18

# --- Step 6: Plot ---
p_acme_vs_peak_distance <- ggplot(pfas_summary, aes(
  x = mean_abs_acme, y = mean_peak_distance,
  color = exposure_label, fill = exposure_label
)) +
  
  # Error bars reflect SE across iterations
  geom_errorbar(aes(ymin = mean_peak_distance - se_peak_distance,
                    ymax = mean_peak_distance + se_peak_distance),
                width = 0.005, alpha = 0.3) +
  geom_errorbarh(aes(xmin = mean_abs_acme - se_abs_acme,
                     xmax = mean_abs_acme + se_abs_acme),
                 height = 0.2, alpha = 0.3) +
  
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(aes(group = exposure_label), method = "lm", se = TRUE,
              alpha = 0.15, linewidth = 1.2) +
  
  geom_text_repel(aes(label = PFAS),
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0,
                  show.legend = FALSE) +
  
  geom_label(data = pfas_stats %>% filter(exposure_label == "Fetal"),
             aes(x = 0.02, y = label_y_start, label = label),
             hjust = 0, vjust = 1,
             size = 3.8, label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE, color = "#4682B4", fill = "white") +
  
  geom_label(data = pfas_stats %>% filter(exposure_label == "Maternal"),
             aes(x = 0.02, y = label_y_start - label_y_spacing, label = label),
             hjust = 0, vjust = 1,
             size = 3.8, label.padding = unit(0.3, "lines"),
             inherit.aes = FALSE, color = "#F08080", fill = "white") +
  
  scale_color_manual(values = c("Maternal" = "#F08080", "Fetal" = "#4682B4")) +
  scale_fill_manual(values = c("Maternal" = "#F08080", "Fetal" = "#4682B4")) +
  
  labs(
    x = "Average |ACME| ± SE",
    y = "Average peak distance between\nmediator and network hubs ± SE",
    title = "PFAS-gestational age mediation strength\ndoes not scale with network centrality"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "top")

p_acme_vs_peak_distance

ggsave("p_acme_vs_peak_distance_avg_GA.pdf", p_acme_vs_peak_distance,
       width = 6.75, height = 6, units = "in", dpi = 300)

