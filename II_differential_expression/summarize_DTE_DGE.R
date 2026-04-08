library(tidyverse)
library(ggplot2)
library(ggview)

load("tx2g_ESPRESSO_assembly_SC_filtered.RData")

tx2g <- tx2g.assembly

# Define the directory and PFAS compounds
results_dir <- "/Users/stbresnahan/Desktop/PFAS/results"
pfas_compounds <- c("PFBA", "PFBS", "PFDA", "PFHxS", "PFOA", "PFOS", "PFNA", "PFUnDA")


# MaternalS

# Initialize empty list to store results
summary_list <- list()

# Process each PFAS compound
for (pfas in pfas_compounds) {
  
  # Read DTE file
  dte_file <- file.path(results_dir, paste0(pfas, "_mat_DTE.tsv"))
  dte_data <- read_tsv(dte_file, show_col_types = FALSE)
  
  # Filter significant DTE results
  dte_sig <- dte_data %>%
    filter(FDR < 0.1, abs(log2FC) > 1)
  
  # Convert transcripts to genes using tx2g
  dte_genes <- dte_sig %>%
    left_join(tx2g, by = c("Transcript" = "transcript_id")) %>%
    pull(gene_id) %>%
    unique() %>%
    length()
  
  # Read DGE file
  dge_file <- file.path(results_dir, paste0(pfas, "_mat_DGE.tsv"))
  dge_data <- read_tsv(dge_file, show_col_types = FALSE)
  
  # Filter significant DGE results
  dge_sig <- dge_data %>%
    filter(FDR < 0.1, abs(log2FC) > 1) %>%
    pull(`Ensembl`) %>%
    unique() %>%
    length()
  
  # Add to summary
  summary_list[[length(summary_list) + 1]] <- data.frame(
    PFAS = pfas,
    level = "DTE",
    significant = dte_genes
  )
  
  summary_list[[length(summary_list) + 1]] <- data.frame(
    PFAS = pfas,
    level = "DGE",
    significant = dge_sig
  )
}

# Combine into final table
pfas_summary.mat <- bind_rows(summary_list)

# Create grouped bar plot
g1 <- ggplot(pfas_summary.mat, aes(x = PFAS, y = log2(significant+1), fill = level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Differential expression in placenta by maternal PFAS exposure",
    y = "Significant Effects (log2)",
    fill = "Analysis Level",
    subtitle = "Significant = FDR < 10%, |log2FC| > 1"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank()
  ) +
  scale_fill_manual(
    values = c("DTE" = "#2E86AB", "DGE" = "#A23B72"),
    labels = c("Genes with gene-level effects","Genes with isoform-level effects")
  )
g1




# CHILD EXPOSURES

# Initialize empty list to store results
summary_list <- list()

# Process each PFAS compound
for (pfas in pfas_compounds) {
  
  # Read DTE file
  dte_file <- file.path(results_dir, paste0(pfas, "_child_DTE.tsv"))
  dte_data <- read_tsv(dte_file, show_col_types = FALSE)
  
  # Filter significant DTE results
  dte_sig <- dte_data %>%
    filter(FDR < 0.1, abs(log2FC) > 1)
  
  # Convert transcripts to genes using tx2g
  dte_genes <- dte_sig %>%
    left_join(tx2g, by = c("Transcript" = "transcript_id")) %>%
    pull(gene_id) %>%
    unique() %>%
    length()
  
  # Read DGE file
  dge_file <- file.path(results_dir, paste0(pfas, "_child_DGE.tsv"))
  dge_data <- read_tsv(dge_file, show_col_types = FALSE)
  
  # Filter significant DGE results
  dge_sig <- dge_data %>%
    filter(FDR < 0.1, abs(log2FC) > 1) %>%
    pull(`Ensembl`) %>%
    unique() %>%
    length()
  
  # Add to summary
  summary_list[[length(summary_list) + 1]] <- data.frame(
    PFAS = pfas,
    level = "DTE",
    significant = dte_genes
  )
  
  summary_list[[length(summary_list) + 1]] <- data.frame(
    PFAS = pfas,
    level = "DGE",
    significant = dge_sig
  )
}

# Combine into final table
pfas_summary.child <- bind_rows(summary_list)

# Create grouped bar plot
g2 <- ggplot(pfas_summary.child, aes(x = PFAS, y = log2(significant+1), fill = level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Differential expression in placenta by fetal PFAS exposure",
    y = "Significant Effects (log2)",
    fill = "Analysis Level",
    subtitle = "Significant = FDR < 10%, |log2FC| > 1"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank()
  ) +
  scale_fill_manual(
    values = c("DTE" = "#2E86AB", "DGE" = "#A23B72"),
    labels = c("Genes with gene-level effects","Genes with isoform-level effects")
  )
g2




# Combined
pfas_summary.mat$group <- "Maternal"
pfas_summary.child$group <- "Fetal"
pfas_summary <- rbind(pfas_summary.mat,pfas_summary.child)


# Create grouped bar plot
library(ggbreak)
library(dplyr)
library(ggplot2)

# Transform the data
# 0-50 takes up 2/3 of space (scale factor = 1)
# 250-300 takes up 1/3 of space (scale factor = 50/3 per unit ≈ 16.67)
# So map 250-300 to 50-53 range (compressed)

pfas_summary_cut <- pfas_summary %>%
  mutate(
    significant_display = case_when(
      significant <= 50 ~ significant,
      significant >= 250 ~ 50 + (significant - 250) * (50/3) / 50,  # Map 250-300 to 50-83.33
      TRUE ~ NA_real_
    )
  )

# Calculate the actual displayed range
# 0-50 = 2/3, so total range should be 75 (50 + 25)
# Map 250-300 to 50-75

pfas_summary_cut <- pfas_summary %>%
  mutate(
    significant_display = case_when(
      significant <= 50 ~ significant,
      significant >= 250 ~ 50 + (significant - 250) / 2,  # Map 250-300 to 50-75
      TRUE ~ NA_real_
    )
  )

# Load TPTE data
tpte <- read.csv("pfas_cors.csv")
tpte_values <- tpte %>%
  dplyr::select(PFAS, mean_tpte) %>%
  distinct()

# Order PFAS by mean_tpte
pfas_order <- tpte_values %>%
  arrange(mean_tpte) %>%
  pull(PFAS)

# Add TPTE values to pfas_summary_cut and set factor levels
pfas_summary_cut <- pfas_summary_cut %>%
  left_join(tpte_values, by = "PFAS") %>%
  mutate(PFAS = factor(PFAS, levels = pfas_order))

g3 <- ggplot(pfas_summary_cut, aes(x = PFAS, y = significant_display, fill = level)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  # Add horizontal white line at break
  annotate("segment", x = -Inf, xend = Inf, y = 50, yend = 50,
           linetype = "solid", color = "white", linewidth = 4) +
  # Add single // symbol on the y-axis only
  annotate("text", x = 0.6, y = 50, label = "//", 
           size = 5, angle = 60, fontface = "bold", color = "black") +
  scale_y_continuous(
    breaks = c(0, 25, 50, 62.5, 75),
    labels = c("0", "25", "50", "275", "300"),
    limits = c(0, 77),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "PFAS-associated differential expression in placenta\nis more extensive for maternal than for fetal exposures",
    subtitle = "FDR < 10%, |log2FC| > 1",
    y = "Significant Effects",
    fill = "Analysis Level") +
  scale_fill_manual(
    values = c("DGE" = "#A23B72", "DTE" = "#2E86AB"),
    labels = c("Genes with\ngene-level effects", "Genes with\nisoform-level effects")
  ) +
  guides(fill = guide_legend(byrow = TRUE)) +
  facet_wrap(~group) +
  theme_minimal() +
  theme(
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10, hjust=0),
    plot.subtitle = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_blank(),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # facet borders,
    axis.title.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90")
  ) +
  canvas(3.7, 3.2, units = "in")
g3

save_ggplot(g3, "PFAS_effects_comparison.pdf", dpi=300)



##### Trans-placental transfer efficiency vs DTE & DGE
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(ggrepel)

# Initialize results for both DTE and DGE
cor_results_dte <- list()
cor_results_dge <- list()

# Loop through PFAS compounds
for (pfas in pfas_compounds) {
  
  # ---- Process DTE ----
  dte_mat_file <- file.path(results_dir, paste0(pfas, "_mat_DTE.tsv"))
  dte_child_file <- file.path(results_dir, paste0(pfas, "_child_DTE.tsv"))
  
  if (file.exists(dte_mat_file) & file.exists(dte_child_file)) {
    dte_mat <- read_tsv(dte_mat_file, show_col_types = FALSE) %>%
      select(Transcript, log2FC_mat = log2FC)
    
    dte_child <- read_tsv(dte_child_file, show_col_types = FALSE) %>%
      select(Transcript, log2FC_child = log2FC)
    
    joined_dte <- inner_join(dte_mat, dte_child, by = "Transcript")
    
    if (nrow(joined_dte) >= 5) {
      cor_test_dte <- cor.test(joined_dte$log2FC_mat, joined_dte$log2FC_child, method = "pearson")
      
      cor_results_dte[[pfas]] <- tibble(
        PFAS = pfas,
        estimate = cor_test_dte$estimate[[1]],
        p_value = cor_test_dte$p.value,
        conf_low = cor_test_dte$conf.int[1],
        conf_high = cor_test_dte$conf.int[2],
        n_features = nrow(joined_dte),
        analysis = "DTE"
      )
    }
  }
  
  # ---- Process DGE ----
  dge_mat_file <- file.path(results_dir, paste0(pfas, "_mat_DGE.tsv"))
  dge_child_file <- file.path(results_dir, paste0(pfas, "_child_DGE.tsv"))
  
  if (file.exists(dge_mat_file) & file.exists(dge_child_file)) {
    dge_mat <- read_tsv(dge_mat_file, show_col_types = FALSE) %>%
      select(Ensembl, log2FC_mat = log2FC)
    
    dge_child <- read_tsv(dge_child_file, show_col_types = FALSE) %>%
      select(Ensembl, log2FC_child = log2FC)
    
    joined_dge <- inner_join(dge_mat, dge_child, by = "Ensembl")
    
    if (nrow(joined_dge) >= 5) {
      cor_test_dge <- cor.test(joined_dge$log2FC_mat, joined_dge$log2FC_child, method = "pearson")
      
      cor_results_dge[[pfas]] <- tibble(
        PFAS = pfas,
        estimate = cor_test_dge$estimate[[1]],
        p_value = cor_test_dge$p.value,
        conf_low = cor_test_dge$conf.int[1],
        conf_high = cor_test_dge$conf.int[2],
        n_features = nrow(joined_dge),
        analysis = "DGE"
      )
    }
  }
}

# Combine results
correlation_summary <- bind_rows(
  bind_rows(cor_results_dte),
  bind_rows(cor_results_dge)
)

# Merge TPTE values
tpte <- read.csv("pfas_cors.csv")
tpte_values <- tpte %>%
  dplyr::select(PFAS, mean_tpte) %>%
  distinct()

correlation_summary <- left_join(correlation_summary, tpte_values, by = "PFAS")

# Prepare data for plotting
tpte_long <- correlation_summary %>%
  dplyr::select(PFAS, mean_tpte, estimate, conf_low, conf_high, analysis) %>%
  rename(r = estimate)

# Fit separate regressions for DTE and DGE
reg_stats <- tpte_long %>%
  group_by(analysis) %>%
  do({
    model <- lm(r ~ mean_tpte, data = .)
    glance_model <- glance(model)
    cor_test <- cor.test(.$r, .$mean_tpte, method = "pearson")
    data.frame(
      analysis = unique(.$analysis),
      pearson_r = cor_test$estimate,
      r_squared = glance_model$r.squared,
      reg_p_value = glance_model$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    x_pos = max(tpte_long$mean_tpte) * 0.95,
    y_pos = ifelse(analysis == "DTE", 0.7, 0.55)
  )

# Recode analysis levels to DEGs and DETs
tpte_long <- tpte_long %>%
  mutate(analysis = recode(analysis, "DGE" = "DEGs", "DTE" = "DETs"))

reg_stats <- reg_stats %>%
  mutate(analysis = recode(analysis, "DGE" = "DEGs", "DTE" = "DETs"))

# Define colors
analysis_colors <- c("DEGs" = "#A23B72", "DETs" = "#2E86AB")

# ---- Final Plot ----
g_tpte_combined <- ggplot(tpte_long, aes(x = mean_tpte, y = r, color = analysis, fill = analysis)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 1) +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high),
                width = 0.1, linewidth = 0.6,
                position = position_dodge(width = 0.05)) +
  geom_point(size = 2.5, alpha = 1,
             position = position_dodge(width = 0.05)) +
  geom_text_repel(aes(label = PFAS),
                  size = 3.5,
                  max.overlaps = Inf,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  segment.size = 0.3,
                  min.segment.length = 0,
                  show.legend = FALSE) +
  geom_label(data = reg_stats,
             aes(x = x_pos + 0.1, y = y_pos + 0.05,
                 label = sprintf("%s: r = %.3f, R² = %.3f, p = %.3f",
                                 analysis, pearson_r, r_squared, reg_p_value),
                 color = analysis),
             fill = "white",
             size = 3.5,
             fontface = "bold",
             label.size = 0.5,
             label.padding = unit(0.3, "lines"),
             hjust = 1,
             show.legend = FALSE) +
  scale_color_manual(values = analysis_colors, name = "Analysis") +
  scale_fill_manual(values = analysis_colors, name = "Analysis") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  labs(
    x = "TPTE",
    y = "Correlation (r)",
    title = "Maternal-fetal correlation of transcript effect estimates\ndecreases with PFAS transplacental transfer efficiency"
  )

# Display and save
g_tpte_combined

ggsave("TPTE_DTE_DGE_combined.pdf", g_tpte_combined, width = 6.75, height = 6, units = "in", dpi = 300)
