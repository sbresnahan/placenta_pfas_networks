library(ggplot2)
library(dplyr)
library(ggtext)
library(ggforce)
library(ggh4x)


covs_PFAS <- readRDS("covs_PFAS.rds")

covs_PFAS <- covs_PFAS[covs_PFAS$SampleID != "J1042", ]

pfas <- c("PFBA","PFOA","PFNA","PFDA","PFUnDA","PFBS","PFHxS","PFOS")

get_cor_stats <- function(var){
       
         child_var <- paste0(var, "_child")
         mat_var   <- paste0(var, "_mat")
         
           # Keep complete cases
           df <- covs_PFAS[complete.cases(covs_PFAS[, c(child_var, mat_var)]), ]
           
             test <- cor.test(df[[child_var]], df[[mat_var]], method = "pearson")
             
               data.frame(
                     PFAS = var,
                     N    = nrow(df),
                     r    = unname(test$estimate),
                     R2   = unname(test$estimate)^2,
                     p    = test$p.value
                 )
}

results <- do.call(rbind, lapply(pfas, get_cor_stats))

order_vec <- c("PFOS", "PFDA", "PFUnDA", "PFHxS", 
               "PFNA", "PFOA", "PFBA", "PFBS")



plot_data <- do.call(rbind, lapply(pfas, function(var) {
  child_var <- paste0(var, "_child")
  mat_var   <- paste0(var, "_mat")
  df <- covs_PFAS[complete.cases(covs_PFAS[, c(child_var, mat_var)]), ]
  data.frame(PFAS = var, child = df[[child_var]], mat = df[[mat_var]])
}))


p_label <- ifelse(results$p < 0.001, "p < 0.001", paste0("p=", round(results$p, 2)))

results$label <- paste0("**", results$PFAS, "**<br>r=", round(results$r, 2),
                        ", R²=", round(results$R2, 2),
                        ", ", p_label)

label_map <- setNames(results$label, results$PFAS)

plot_data$PFAS <- factor(plot_data$PFAS, levels = order_vec)
plot_data$facet_label <- factor(label_map[as.character(plot_data$PFAS)],
                                levels = label_map[order_vec])

g1 <- ggplot(plot_data, aes(x = mat, y = child)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  facet_wrap(~ facet_label, scales = "free", nrow = 2) +
  labs(x = "Maternal PFAS z-score", y = "Fetal PFAS z-score") +
  theme_bw() +
  theme(strip.text = element_markdown())

g1

ggsave("PFAS_cors.pdf",g1,width=10.75,height=5.75,units="in",dpi=300)



# Final summary figure as forest plot of slopes & CIs

bw_coex <- readRDS("bw_coex_slopes.RDS")
names(bw_coex)[1] <- "exposure_window"
bw_coex$outcome <- "Birth weight"
bw_coex$metric <- "Co-expression strength"

ga_coex <- readRDS("ga_coex_slopes.RDS")
names(ga_coex)[1] <- "exposure_window"
ga_coex$outcome <- "Gestational age"
ga_coex$metric <- "Co-expression strength"

bw_mc <- readRDS("bw_mc_slopes.RDS")
names(bw_mc)[1] <- "exposure_window"
bw_mc$outcome <- "Birth weight"
bw_mc$metric <- "Mediator count"

ga_mc <- readRDS("ga_mc_slopes.RDS")
names(ga_mc)[1] <- "exposure_window"
ga_mc$outcome <- "Gestational age"
ga_mc$metric <- "Mediator count"

bw_nc <- readRDS("bw_nc_slopes.RDS")
names(bw_nc)[1] <- "exposure_window"
bw_nc$outcome <- "Birth weight"
bw_nc$metric <- "Centrality"

ga_nc <- readRDS("ga_nc_slopes.RDS")
names(ga_nc)[1] <- "exposure_window"
ga_nc$outcome <- "Gestational age"
ga_nc$metric <- "Centrality"

bw_comp <- readRDS("bw_comp_slopes.RDS")
bw_comp$exposure_window <- NA
bw_comp$outcome <- "Birth weight"
bw_comp$metric <- "Compartmentalization"
bw_comp <- bw_comp[,names(bw_nc)]

ga_comp <- readRDS("ga_comp_slopes.RDS")
ga_comp$exposure_window <- NA
ga_comp$outcome <- "Gestational age"
ga_comp$metric <- "Compartmentalization"
ga_comp <- ga_comp[,names(ga_nc)]

slopes_fitted <- rbind(bw_coex,ga_coex,
                       bw_mc,ga_mc,
                       bw_nc,ga_nc,
                       bw_comp,ga_comp)
row.names(slopes_fitted) <- NULL

slopes_fitted <- slopes_fitted %>%
  mutate(group = ifelse(is.na(exposure_window),
                        metric,
                        paste0(metric, " (", exposure_window, ")")),
         is_compartment = is.na(exposure_window))

forest_data <- data.frame(slopes_fitted)
forest_data[forest_data$metric=="Centrality","slope"] <- -forest_data[forest_data$metric=="Centrality","slope"]
forest_data[forest_data$metric=="Centrality","slope_lwr"] <- -forest_data[forest_data$metric=="Centrality","slope_lwr"]
forest_data[forest_data$metric=="Centrality","slope_upr"] <- -forest_data[forest_data$metric=="Centrality","slope_upr"]

# Reorder the metric factor so facets appear in desired order
forest_data$metric <- factor(
  forest_data$metric,
  levels = c("Mediator count", "Co-expression strength", "Centrality", "Compartmentalization")
)

p <- ggplot(filter(forest_data, !is_compartment), 
            aes(x = slope, y = exposure_window, color = exposure_window)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = slope_lwr, xmax = slope_upr), height = 0.2) +
  geom_point(aes(shape = significant), size = 3) +
  geom_errorbarh(data = filter(forest_data, is_compartment),
                 aes(x = slope, y = 1.5, xmin = slope_lwr, xmax = slope_upr),
                 height = 0.2, color = "black") +
  geom_point(data = filter(forest_data, is_compartment),
             aes(x = slope, y = 1.5, shape = significant),
             color = "black", size = 3) +
  facet_grid2(outcome ~ metric, scales = "free_x") +
  scale_color_manual(
    values = c("Fetal" = "#4682B4", "Maternal" = "#F08080"),
    name = NULL
  ) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), guide = "none") +  # hide significance legend
  labs(x = "Slope (change per unit TPTE)", y = NULL) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p

ggsave("network_metrics_forest.pdf", p, width = 8.5, height = 2.75, units="in", dpi = 300)







# Compare mediators

# Libraries
library(dplyr)
library(ggplot2)

covs_PFAS <- readRDS("covs_PFAS.rds")

# Get list of all result files (BW)
result_files <- list.files("WGCNA_mediation_results", pattern = "mediation_batch_.*\\.rds", full.names = TRUE)
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
all_module_mediation_df_bw <- all_module_mediation_df %>%
  left_join(module_transcripts_base %>% dplyr::select(transcript_id, module, kME), 
            by = c("transcript_id", "module"))



# Get list of all result files (GA)
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
all_module_mediation_df_GA <- all_module_mediation_df



# Merge
all_module_mediation_df_bw$outcome <- "birth_weight"
all_module_mediation_df_GA$outcome <- "gestational_age"

mediation_results <- rbind(all_module_mediation_df_bw,all_module_mediation_df_GA)


library(UpSetR)

# Prepare data for UpSet plot
upset_data <- mediation_results %>%
  filter(significant_strict == TRUE) %>%
  mutate(group = paste(exposure_label, outcome, sep = "_")) %>%
  select(transcript_id, group) %>%
  distinct() %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = group,
    values_from = present,
    values_fill = 0
  ) %>%
  as.data.frame()

# Remove transcript_id column for the plot (keep as rownames if needed)
rownames(upset_data) <- upset_data$transcript_id
upset_data <- upset_data[, -1]

# Create UpSet plot

pdf("upset_plot_significant_transcripts.pdf", width = 7.75, height = 5)

upset(
  upset_data,
  sets = colnames(upset_data),
  order.by = "freq",
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = 1.5,
  point.size = 3,
  line.size = 1,
  main.bar.color = "black",
  sets.bar.color = "black"
)

dev.off()

# Find transcripts present in all sets
upset_data %>%
  rownames_to_column("transcript_id") %>%
  filter(rowSums(select(., -transcript_id)) == ncol(upset_data)) %>%
  pull(transcript_id) # [1] "ENST00000544179.1" "ENST00000683927.1" "ENST00000509648.5"

classif <- read.table("/Users/stbresnahan/Desktop/Manuscripts/Placenta/Submission/lr-placenta-transcriptome-main_freeze/Assembly/ESPRESSO_corrected_SC_filtered_classification.txt",header=T)

# Correlation of gestational_age and birth_weight

covs_PFAS.sub <- covs_PFAS[covs_PFAS$labour_onset=="1_Spontaneous",]

# Calculate regression statistics
model <- lm(c_weight_birth_kg ~ GA, data = covs_PFAS.sub)
r <- cor(covs_PFAS.sub$GA, covs_PFAS.sub$c_weight_birth_kg, use = "complete.obs")
r_squared <- summary(model)$r.squared
p_value <- summary(model)$coefficients[2, 4]

# Create label
stat_label <- sprintf("r = %.3f, R² = %.3f, p = %.3e", r, r_squared, p_value)

# Create plot
p <- ggplot(covs_PFAS.sub, aes(x = GA, y = c_weight_birth_kg)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  annotate("label",
           x = -2.8,
           y = max(covs_PFAS$c_weight_birth_kg, na.rm = TRUE),
           label = stat_label,
           hjust = 0, vjust = 1,
           size = 4,
           fill = "white",
           color = "black",
           label.padding = unit(0.5, "lines")) +
  labs(x = "Gestational Age (weeks) z-score",
       y = "Birth Weight (kg) z-score",
       title = "Birth weight is positively correlated with gestational age") +
  theme_classic(base_size = 12)
p

ggsave("bw_GA_corr.pdf",p,width=5.5,height=5.5,units="in")



# Proportion unique genes by TPTE

tx2g <- data.frame(transcript_id=classif$isoform,gene_id=classif$associated_gene)
mediation_results <- left_join(mediation_results,tx2g)

library(dplyr)
library(ggplot2)
library(ggrepel)

tpte_lookup <- tibble(
  PFAS = c("PFBA", "PFOA", "PFNA", "PFDA", "PFUnDA", "PFBS", "PFHxS", "PFOS"),
  tpte = c(1.97, 1.17, 0.99, 0.43, 0.46, 2.04, 0.64, 0.36)
)

results_list <- list()

for (out in unique(mediation_results$outcome)) {
  for (pfas in unique(mediation_results$PFAS)) {
    
    # Get significant mediator transcripts with their gene mappings
    mat_df <- mediation_results %>%
      filter(outcome == out, PFAS == pfas, exposure_window == "mat", significant_strict == TRUE) %>%
      select(transcript_id, gene_id)
    
    child_df <- mediation_results %>%
      filter(outcome == out, PFAS == pfas, exposure_window == "child", significant_strict == TRUE) %>%
      select(transcript_id, gene_id)
    
    if (nrow(mat_df) == 0 || nrow(child_df) == 0) next
    
    # Unique genes from each set
    mat_genes <- unique(mat_df$gene_id)
    child_genes <- unique(child_df$gene_id)
    shared_genes <- intersect(mat_genes, child_genes)
    union_genes <- union(mat_genes, child_genes)
    
    # For Fisher's test: use all genes that have transcripts tested as background
    all_genes <- mediation_results %>%
      filter(outcome == out, PFAS == pfas) %>%
      pull(gene_id) %>%
      unique()
    
    # 2x2 contingency table:
    # Gene is mat mediator (yes/no) x Gene is child mediator (yes/no)
    a <- length(shared_genes)                                    # mat=yes, child=yes
    b <- length(setdiff(mat_genes, child_genes))                 # mat=yes, child=no
    c <- length(setdiff(child_genes, mat_genes))                 # mat=no, child=yes
    d <- length(setdiff(all_genes, union_genes))                 # mat=no, child=no
    
    fisher_res <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
    
    results_list[[length(results_list) + 1]] <- tibble(
      outcome = out,
      PFAS = pfas,
      n_transcripts_mat = nrow(mat_df),
      n_transcripts_child = nrow(child_df),
      n_genes_mat = length(mat_genes),
      n_genes_child = length(child_genes),
      n_shared = length(shared_genes),
      n_union = length(union_genes),
      n_background = length(all_genes),
      jaccard = length(shared_genes) / length(union_genes),
      odds_ratio = fisher_res$estimate,
      log_or = log(fisher_res$estimate),
      fisher_p = fisher_res$p.value
    )
  }
}

overlap_results <- bind_rows(results_list) %>%
  left_join(tpte_lookup, by = "PFAS")

# Plot 1: Jaccard index
for (out in unique(overlap_results$outcome)) {
  
  plot_data <- overlap_results %>% filter(outcome == out)
  
  tpte_stats <- plot_data %>%
    summarise(
      r = cor(tpte, jaccard, use = "complete.obs"),
      r_squared = summary(lm(jaccard ~ tpte))$r.squared,
      p_value = summary(lm(jaccard ~ tpte))$coefficients[2, 4]
    ) %>%
    mutate(label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value))
  
  label_y_start <- max(plot_data$jaccard, na.rm = TRUE)
  
  p <- plot_data %>%
    ggplot(aes(x = tpte, y = jaccard)) +
    geom_point(size = 4, alpha = 0.85, color = "#666666") +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 1.2, color = "#333333", fill = "#333333") +
    geom_text_repel(
      aes(label = PFAS),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    geom_label(
      data = tpte_stats,
      aes(x = min(plot_data$tpte, na.rm = TRUE), y = label_y_start, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 3.5,
      color = "#333333",
      fill = "white"
    ) +
    labs(
      x = "TPTE (Transplacental Transfer Efficiency)",
      y = "Gene overlap (Jaccard index)",
      title = paste0("Shared mediator genes between maternal and fetal exposure\n(", out, ")"),
      subtitle = "Jaccard = shared genes / union of genes from significant mediator transcripts"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
}

# Plot 2: Fisher's exact test log odds ratios
for (out in unique(overlap_results$outcome)) {
  
  plot_data <- overlap_results %>% 
    filter(outcome == out, is.finite(log_or))
  
  if (nrow(plot_data) < 3) next
  
  tpte_stats <- plot_data %>%
    summarise(
      r = cor(tpte, log_or, use = "complete.obs"),
      r_squared = summary(lm(log_or ~ tpte))$r.squared,
      p_value = summary(lm(log_or ~ tpte))$coefficients[2, 4]
    ) %>%
    mutate(label = sprintf("r = %.3f, R² = %.3f, p = %.3f", r, r_squared, p_value))
  
  label_y_start <- max(plot_data$log_or, na.rm = TRUE)
  
  p <- plot_data %>%
    ggplot(aes(x = tpte, y = log_or)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 4, alpha = 0.85, color = "#666666") +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.1, linewidth = 1.2, color = "#333333", fill = "#333333") +
    geom_text_repel(
      aes(label = paste0(PFAS, ifelse(fisher_p < 0.05, "*", ""))),
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50",
      segment.size = 0.3
    ) +
    geom_label(
      data = tpte_stats,
      aes(x = min(plot_data$tpte, na.rm = TRUE), y = label_y_start, label = label),
      inherit.aes = FALSE,
      hjust = 0, vjust = 1,
      size = 3.5,
      color = "#333333",
      fill = "white"
    ) +
    labs(
      x = "TPTE (Transplacental Transfer Efficiency)",
      y = "Log odds ratio (Fisher's exact test)",
      title = paste0("Gene-level overlap between maternal and fetal mediators\n(", out, ")"),
      subtitle = "Positive = more gene sharing than expected by chance; * = p < 0.05"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
}

# Summary table
overlap_results %>%
  select(outcome, PFAS, tpte, n_transcripts_mat, n_transcripts_child, 
         n_genes_mat, n_genes_child, n_shared, n_background, jaccard, log_or, fisher_p) %>%
  arrange(outcome, desc(tpte))
