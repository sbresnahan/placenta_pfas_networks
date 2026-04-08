# Overlapping DGE/DTE results between Elkin placental explants with GUSTO term placentae

setwd("/Users/stbresnahan/Desktop/PFAS")

library(ggplot2)
library(dplyr)
library(ggview)

# Thresholds
lfc.threshold=1
p.threshold=0.05


# Helper to extract contingency counts

counts_df <- data.frame()

extract_counts <- function(explant_sig, tissue_sig, overlap, universe, label) {
  universe    <- unique(universe)
  explant_sig <- intersect(explant_sig, universe)
  tissue_sig  <- intersect(tissue_sig, universe)
  overlap     <- intersect(overlap, universe)
  
  in_explant <- universe %in% explant_sig
  in_tissue  <- universe %in% tissue_sig
  in_overlap <- universe %in% overlap
  
  data.frame(
    concordant   = sum(in_overlap),
    explant_only = sum(in_explant & !in_overlap),
    tissue_only  = sum(in_tissue & !in_overlap),
    neither      = sum(!in_explant & !in_tissue),
    label        = label
  )
}

# Load data

load("tx2g_ESPRESSO_assembly_SC_filtered.RData")

# intersect_feature <- function(x,y,q.threshold=0.1,lfc.threshold=.05){
#   return(intersect(x[abs(x$log2FC)>lfc.threshold&x$P<q.threshold,1],
#                    y[abs(y$log2FC)>lfc.threshold&y$P<q.threshold,1]))
# }

intersect_feature <- function(x, y, q.threshold = 0.1, lfc.threshold = 0) {
  # Get features with significant effects in x (positive direction)
  x_pos <- x[x$log2FC > lfc.threshold & x$P < q.threshold, 1]
  # Get features with significant effects in x (negative direction)
  x_neg <- x[x$log2FC < -lfc.threshold & x$P < q.threshold, 1]
  
  # Get features with significant effects in y (positive direction)
  y_pos <- y[y$log2FC > lfc.threshold & y$P < q.threshold, 1]
  # Get features with significant effects in y (negative direction)
  y_neg <- y[y$log2FC < -lfc.threshold & y$P < q.threshold, 1]
  
  # Return features that are significant in the same direction in both datasets
  same_direction <- c(intersect(x_pos, y_pos), intersect(x_neg, y_neg))
  
  return(same_direction)
}

run_fisher <- function(target, background, universe) {
  target <- intersect(target, universe)
  background <- intersect(background, universe)
  
  # construct contingency table
  in_target <- universe %in% target
  in_bg     <- universe %in% background
  
  mat <- matrix(c(
    sum(in_target & in_bg),       # overlap
    sum(in_target & !in_bg),      # in target only
    sum(!in_target & in_bg),      # in bg only
    sum(!in_target & !in_bg)      # neither
  ), nrow = 2, byrow = TRUE)
  
  fisher.test(mat)
}

# lr-assembly

## DGE
res.elkin <- read.table("Elkin2025_Placenta/results/drug_concentration_numeric_DGE_1.tsv",sep="\t",header=T)
res.gusto.mat <- read.table("results/PFBS_mat_DGE.tsv",sep="\t",header=T)
res.gusto.child <- read.table("results/PFBS_child_DGE.tsv",sep="\t",header=T)

elkin.sig <- res.elkin[abs(res.elkin$log2FC)>lfc.threshold&res.elkin$P<p.threshold,"Ensembl"]
gusto.mat.sig <- res.gusto.mat[abs(res.gusto.mat$log2FC)>lfc.threshold&res.gusto.mat$P<p.threshold,"Ensembl"]
gusto.child.sig <- res.gusto.child[abs(res.gusto.child$log2FC)>lfc.threshold&res.gusto.child$P<p.threshold,"Ensembl"]

ol.em.loose <- intersect_feature(res.elkin,res.gusto.mat)

ol.ec.loose <- intersect_feature(res.elkin,res.gusto.child)

ol.mc.loose <- intersect_feature(res.gusto.mat,res.gusto.child)

results <- list(
  em = run_fisher(ol.em.loose, elkin.sig, tx2g.assembly$gene_id),
  ec = run_fisher(ol.ec.loose, elkin.sig, tx2g.assembly$gene_id),
  mc = run_fisher(ol.mc.loose, elkin.sig, tx2g.assembly$gene_id)
)

res_df.g <- data.frame(
  contrast = names(results),
  p.value  = sapply(results, function(x) x$p.value),
  odds.ratio = sapply(results, function(x) x$estimate),
  conf.low  = sapply(results, function(x) x$conf.int[1]),
  conf.high = sapply(results, function(x) x$conf.int[2]),
  stringsAsFactors = FALSE
)

res_df.g

counts_df <- bind_rows(counts_df,
                       extract_counts(elkin.sig, gusto.mat.sig,   ol.em.loose, tx2g.assembly$gene_id, "lr-assembly|DGE|Explants vs Tissue (Maternal)"),
                       extract_counts(elkin.sig, gusto.child.sig, ol.ec.loose, tx2g.assembly$gene_id, "lr-assembly|DGE|Explants vs Tissue (Fetal)"))

## DTE
res.elkin <- read.table("Elkin2025_Placenta/results/drug_concentration_numeric_DTE_1.tsv",sep="\t",header=T)
res.gusto.mat <- read.table("results/PFBS_mat_DTE.tsv",sep="\t",header=T)
res.gusto.child <- read.table("results/PFBS_child_DTE.tsv",sep="\t",header=T)

elkin.sig <- res.elkin[abs(res.elkin$log2FC)>lfc.threshold&res.elkin$P<p.threshold,"Ensembl"]
gusto.mat.sig <- res.gusto.mat[abs(res.gusto.mat$log2FC)>lfc.threshold&res.gusto.mat$P<p.threshold,"Ensembl"]
gusto.child.sig <- res.gusto.child[abs(res.gusto.child$log2FC)>lfc.threshold&res.gusto.child$P<p.threshold,"Ensembl"]

ol.em.loose <- intersect_feature(res.elkin,res.gusto.mat)

ol.ec.loose <- intersect_feature(res.elkin,res.gusto.child)

ol.mc.loose <- intersect_feature(res.gusto.mat,res.gusto.child)

results <- list(
  em = run_fisher(ol.em.loose, elkin.sig, tx2g.assembly$gene_id),
  ec = run_fisher(ol.ec.loose, elkin.sig, tx2g.assembly$gene_id),
  mc = run_fisher(ol.mc.loose, elkin.sig, tx2g.assembly$gene_id)
)

res_df.tx <- data.frame(
  contrast = names(results),
  p.value  = sapply(results, function(x) x$p.value),
  odds.ratio = sapply(results, function(x) x$estimate),
  conf.low  = sapply(results, function(x) x$conf.int[1]),
  conf.high = sapply(results, function(x) x$conf.int[2]),
  stringsAsFactors = FALSE
)

res_df.tx

counts_df <- bind_rows(counts_df,
                       extract_counts(elkin.sig, gusto.mat.sig,   ol.em.loose, tx2g.assembly$gene_id, "lr-assembly|DTE|Explants vs Tissue (Maternal)"),
                       extract_counts(elkin.sig, gusto.child.sig, ol.ec.loose, tx2g.assembly$gene_id, "lr-assembly|DTE|Explants vs Tissue (Fetal)"))

# GENCODE

## DGE
res.elkin <- read.table("Elkin2025_Placenta/results_gencode/drug_concentration_numeric_DGE_1.tsv",sep="\t",header=T)
res.gusto.mat <- read.table("results_gencode/PFBS_mat_DGE.tsv",sep="\t",header=T)
res.gusto.child <- read.table("results_gencode/PFBS_child_DGE.tsv",sep="\t",header=T)

elkin.sig <- res.elkin[abs(res.elkin$log2FC)>lfc.threshold&res.elkin$P<p.threshold,"Ensembl"]
gusto.mat.sig <- res.gusto.mat[abs(res.gusto.mat$log2FC)>lfc.threshold&res.gusto.mat$P<p.threshold,"Ensembl"]
gusto.child.sig <- res.gusto.child[abs(res.gusto.child$log2FC)>lfc.threshold&res.gusto.child$P<p.threshold,"Ensembl"]

ol.em.loose <- intersect_feature(res.elkin,res.gusto.mat)

ol.ec.loose <- intersect_feature(res.elkin,res.gusto.child)

ol.mc.loose <- intersect_feature(res.gusto.mat,res.gusto.child)

results <- list(
  em = run_fisher(ol.em.loose, elkin.sig, tx2g.gencode$gene_id),
  ec = run_fisher(ol.ec.loose, elkin.sig, tx2g.gencode$gene_id),
  mc = run_fisher(ol.mc.loose, elkin.sig, tx2g.gencode$gene_id)
)

res_df.g.gencode <- data.frame(
  contrast = names(results),
  p.value  = sapply(results, function(x) x$p.value),
  odds.ratio = sapply(results, function(x) x$estimate),
  conf.low  = sapply(results, function(x) x$conf.int[1]),
  conf.high = sapply(results, function(x) x$conf.int[2]),
  stringsAsFactors = FALSE
)

res_df.g.gencode

counts_df <- bind_rows(counts_df,
                       extract_counts(elkin.sig, gusto.mat.sig,   ol.em.loose, tx2g.assembly$gene_id, "gencode|DGE|Explants vs Tissue (Maternal)"),
                       extract_counts(elkin.sig, gusto.child.sig, ol.ec.loose, tx2g.assembly$gene_id, "gencode|DGE|Explants vs Tissue (Fetal)"))

## DTE
res.elkin <- read.table("Elkin2025_Placenta/results_gencode/drug_concentration_numeric_DTE_1.tsv",sep="\t",header=T)
res.gusto.mat <- read.table("results_gencode/PFBS_mat_DTE.tsv",sep="\t",header=T)
res.gusto.child <- read.table("results_gencode/PFBS_child_DTE.tsv",sep="\t",header=T)

elkin.sig <- res.elkin[abs(res.elkin$log2FC)>lfc.threshold&res.elkin$P<p.threshold,"Ensembl"]
gusto.mat.sig <- res.gusto.mat[abs(res.gusto.mat$log2FC)>lfc.threshold&res.gusto.mat$P<p.threshold,"Ensembl"]
gusto.child.sig <- res.gusto.child[abs(res.gusto.child$log2FC)>lfc.threshold&res.gusto.child$P<p.threshold,"Ensembl"]

ol.em.loose <- intersect_feature(res.elkin,res.gusto.mat)

ol.ec.loose <- intersect_feature(res.elkin,res.gusto.child)

ol.mc.loose <- intersect_feature(res.gusto.mat,res.gusto.child)

results <- list(
  em = run_fisher(ol.em.loose, elkin.sig, tx2g.gencode$gene_id),
  ec = run_fisher(ol.ec.loose, elkin.sig, tx2g.gencode$gene_id),
  mc = run_fisher(ol.mc.loose, elkin.sig, tx2g.gencode$gene_id)
)

res_df.tx.gencode <- data.frame(
  contrast = names(results),
  p.value  = sapply(results, function(x) x$p.value),
  odds.ratio = sapply(results, function(x) x$estimate),
  conf.low  = sapply(results, function(x) x$conf.int[1]),
  conf.high = sapply(results, function(x) x$conf.int[2]),
  stringsAsFactors = FALSE
)

res_df.tx.gencode

counts_df <- bind_rows(counts_df,
                       extract_counts(elkin.sig, gusto.mat.sig,   ol.em.loose, tx2g.assembly$gene_id, "gencode|DTE|Explants vs Tissue (Maternal)"),
                       extract_counts(elkin.sig, gusto.child.sig, ol.ec.loose, tx2g.assembly$gene_id, "gencode|DTE|Explants vs Tissue (Fetal)"))

# Combine all results
add_meta_cols <- function(df, name) {
  df %>%
    mutate(
      annotation = ifelse(grepl("gencode", name), "gencode", "lr-assembly"),
      level      = ifelse(grepl("\\.g($|\\.)", name), "DGE", "DTE")
    )
}

res_df.g         <- add_meta_cols(res_df.g, "res_df.g")
res_df.tx        <- add_meta_cols(res_df.tx, "res_df.tx")
res_df.g.gencode <- add_meta_cols(res_df.g.gencode, "res_df.g.gencode")
res_df.tx.gencode<- add_meta_cols(res_df.tx.gencode, "res_df.tx.gencode")

res_all <- bind_rows(
  res_df.g,
  res_df.tx,
  res_df.g.gencode,
  res_df.tx.gencode,
)

row.names(res_all) <- NULL

res_all <- res_all %>%
  mutate(
    contrast   = factor(contrast, levels = c("em", "ec", "mc")),
    annotation = factor(annotation, levels = c("lr-assembly", "gencode")),
    level      = factor(level, levels = c("DGE", "DTE"))
  )

epsilon <- 1e-3  # small pseudocount

res_all <- res_all %>%
  mutate(
    odds.ratio = ifelse(odds.ratio <= 0, epsilon, odds.ratio),
    conf.low   = ifelse(conf.low <= 0, epsilon, conf.low),
    conf.high  = ifelse(conf.high <= 0, epsilon, conf.high)
  )

res_all <- res_all %>%
  mutate(
    contrast = recode(
      contrast,
      "em" = "Explants vs\nTissue (Maternal)",
      "ec" = "Explants vs\nTissue (Fetal)",
      "mc" = "Tissue vs Tissue\n(Maternal vs Fetal)"
    ),
    contrast = factor(
      contrast,
      levels = c(
        "Explants vs\nTissue (Maternal)",
        "Explants vs\nTissue (Fetal)",
        "Tissue vs Tissue\n(Maternal vs Fetal)"
      )
    )
  )

facet_labels <- c(
  "DTE" = "Genes with isoform-level effects",
  "DGE" = "Genes with gene-level effects"
)

res_all.sub <- res_all[-c(3,6,9,12),]
row.names(res_all.sub) <- NULL
res_all.sub[2,3] <- 0.01
res_all.sub[2,4] <- 0.01

p <- ggplot(res_all.sub, aes(x = odds.ratio, y = contrast, color = annotation)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    position = position_dodge(width = 0.5), height = 0.3
  ) +
  scale_color_manual(values=c("lr-assembly" = "#C4663E", "gencode" = "grey50"),name = "Annotation") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(name = "Log odds ratio (overlap enrichment)") +
  facet_wrap(~level,labeller = labeller(level = facet_labels),scales="free_x") +
  labs(subtitle = "FDR < 10%, concordant sign; enrichment vs. |log2FC| > 1, p < 0.05") +
  theme_minimal(base_size = 9) +
  theme(
    plot.title.position = "plot",
    legend.position = "top",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    plot.subtitle=element_text(size=7),
    legend.key.size = unit(0.3, "cm"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt"),
    plot.title = element_text(size = 10, hjust=0),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 7),
    strip.text = element_text(size = 6)
  ) +
  labs(title="Long-read assembly improves concordance between\nexperimental and observational PFBS effects") +
  canvas(3.7, 3.2, units = "in")

p

save_ggplot(p, "PFBS_explants_vs_tissue_P.pdf", dpi=300)



# Upset plot
counts_long <- counts_df %>%
  tidyr::separate(label, into = c("annotation", "level", "contrast"), sep = "\\|") %>%
  tidyr::pivot_longer(cols = c(concordant, explant_only, tissue_only, neither),
                      names_to = "set", values_to = "count") %>%
  mutate(set = recode(set,
                      "concordant"   = "Concordant",
                      "explant_only" = "Explants Only",
                      "tissue_only"  = "Tissue Only",
                      "neither"      = "Neither"
  ),
  set = factor(set, levels = c("Concordant", "Explants Only", "Tissue Only", "Neither")))

p2 <- ggplot(counts_long, aes(x = set, y = count + 1, fill = annotation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = count, vjust = ifelse(annotation == "lr-assembly", -0.3, -1.6)),
            position = position_dodge(width = 0.7), size = 3) +
  scale_fill_manual(values = c("lr-assembly" = "#C4663E", "gencode" = "grey50")) +
  scale_y_continuous(trans = "log2", labels = scales::label_log(base = 2), limits = c(1, 150000)) +
  facet_grid(level ~ contrast) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x     = element_text(angle = 35, hjust = 1),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "top",
    axis.title.x    = element_blank()
  ) +
  labs(y = "Count (log2 scale)", fill = "Annotation",
       title = "Contingency counts underlying Fisher enrichment tests",
       subtitle = "Numbers above bars show raw counts")
p2

ggsave("PFBS_explants_vs_tissue_counts.pdf",p2,width=6.75,height=4.75,units="in",dpi=300,bg="white")





# Sensitivity analysis

library(SummarizedExperiment)
library(DESeq2)
library(apeglm)
library(BiocParallel)

load('elkin_gusto/PFBS_gusto_elkin_DGE_DTE.RData')
load('elkin_gusto/se_elkin.RData')
load('elkin_gusto/se_gusto.RData')
tx2g = readRDS('elkin_gusto/tx2g_assembly.RDS')

assembly.DGE.elkin_sig = subset(assembly.DGE.elkin,
                                P < .05 &
                                  abs(log2FC) > 1)
assembly.DTE.elkin_sig = subset(assembly.DTE.elkin,
                                P < .05 &
                                  abs(log2FC) > 1)
gencode.DGE.elkin_sig = subset(gencode.DGE.elkin,
                               P < .05 &
                                 abs(log2FC) > 1)
gencode.DTE.elkin_sig = subset(gencode.DTE.elkin,
                               P < .05 &
                                 abs(log2FC) > 1)


# --- Gencode DGE ---
gencode_all_merge_DGE = merge(
  merge(gencode.DGE.elkin_sig,
        gencode.DGE.gusto.fetal[,c('Ensembl','log2FC','SE')],
        by = 'Ensembl'),
  gencode.DGE.gusto.maternal[,c('Ensembl','log2FC','SE')],
  by = 'Ensembl')
colnames(gencode_all_merge_DGE)[c(5,6)]   = c('log2FC_Elkin','SE_Elkin')
colnames(gencode_all_merge_DGE)[c(9,10)]  = c('log2FC_GUSTO_Fetal','SE_GUSTO_Fetal')
colnames(gencode_all_merge_DGE)[c(11,12)] = c('log2FC_GUSTO_Maternal','SE_GUSTO_Maternal')

# --- Gencode DTE ---
gencode_all_merge_DTE = merge(
  merge(gencode.DTE.elkin_sig,
        gencode.DTE.gusto.fetal[,c('Transcript','log2FC','SE')],
        by = 'Transcript'),
  gencode.DTE.gusto.maternal[,c('Transcript','log2FC','SE')],
  by = 'Transcript')
colnames(gencode_all_merge_DTE)[c(5,6)]   = c('log2FC_Elkin','SE_Elkin')
colnames(gencode_all_merge_DTE)[c(9,10)]  = c('log2FC_GUSTO_Fetal','SE_GUSTO_Fetal')
colnames(gencode_all_merge_DTE)[c(11,12)] = c('log2FC_GUSTO_Maternal','SE_GUSTO_Maternal')

# --- Assembly DGE ---
assembly_all_merge_DGE = merge(
  merge(assembly.DGE.elkin_sig,
        assembly.DGE.gusto.fetal[,c('Ensembl','log2FC','SE')],
        by = 'Ensembl'),
  assembly.DGE.gusto.maternal[,c('Ensembl','log2FC','SE')],
  by = 'Ensembl')
colnames(assembly_all_merge_DGE)[c(5,6)]   = c('log2FC_Elkin','SE_Elkin')
colnames(assembly_all_merge_DGE)[c(9,10)]  = c('log2FC_GUSTO_Fetal','SE_GUSTO_Fetal')
colnames(assembly_all_merge_DGE)[c(11,12)] = c('log2FC_GUSTO_Maternal','SE_GUSTO_Maternal')

# --- Assembly DTE ---
assembly_all_merge_DTE = merge(
  merge(assembly.DTE.elkin_sig,
        assembly.DTE.gusto.fetal[,c('Transcript','log2FC','SE')],
        by = 'Transcript'),
  assembly.DTE.gusto.maternal[,c('Transcript','log2FC','SE')],
  by = 'Transcript')
colnames(assembly_all_merge_DTE)[c(5,6)]   = c('log2FC_Elkin','SE_Elkin')
colnames(assembly_all_merge_DTE)[c(9,10)]  = c('log2FC_GUSTO_Fetal','SE_GUSTO_Fetal')
colnames(assembly_all_merge_DTE)[c(11,12)] = c('log2FC_GUSTO_Maternal','SE_GUSTO_Maternal')




library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

make_scatter <- function(df, title_label) {
  
  # Pivot to long so fetal & maternal are rows
  df_long <- bind_rows(
    df %>% transmute(
      log2FC_Elkin, SE_Elkin,
      log2FC_GUSTO = log2FC_GUSTO_Fetal,
      SE_GUSTO     = SE_GUSTO_Fetal,
      Exposure     = "Fetal"
    ),
    df %>% transmute(
      log2FC_Elkin, SE_Elkin,
      log2FC_GUSTO = log2FC_GUSTO_Maternal,
      SE_GUSTO     = SE_GUSTO_Maternal,
      Exposure     = "Maternal"
    )
  )
  
  # R² and p per exposure
  stats <- df_long %>%
    group_by(Exposure) %>%
    summarise(
      r2 = cor(log2FC_Elkin, log2FC_GUSTO, use = "complete.obs")^2,
      p  = cor.test(log2FC_Elkin, log2FC_GUSTO)$p.value,
      n  = n(),
      .groups = "drop"
    ) %>%
    mutate(label = sprintf("%s: R² = %.3f, p = %.2e, n = %d", Exposure, r2, p, n))
  
  annotation <- paste(stats$label, collapse = "\n")
  
  ggplot(df_long, aes(x = log2FC_Elkin, y = log2FC_GUSTO, color = Exposure)) +
    geom_errorbar(aes(ymin = log2FC_GUSTO - 1.96 * SE_GUSTO,
                      ymax = log2FC_GUSTO + 1.96 * SE_GUSTO),
                  width = 0, alpha = 0.2) +
    geom_errorbarh(aes(xmin = log2FC_Elkin - 1.96 * SE_Elkin,
                       xmax = log2FC_Elkin + 1.96 * SE_Elkin),
                   height = 0, alpha = 0.2) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
             label = annotation, size = 3, fontface = "italic") +
    labs(title = title_label,
         x = "log2FC (Elkin)", y = "log2FC (GUSTO)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

p1 <- make_scatter(assembly_all_merge_DGE, "Assembly — DGE")
p2 <- make_scatter(assembly_all_merge_DTE, "Assembly — DTE")
p3 <- make_scatter(gencode_all_merge_DGE, "Gencode — DGE")
p4 <- make_scatter(gencode_all_merge_DTE, "Gencode — DTE")

(p1 + p2) / (p3 + p4)








library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Threshold grids
p_thresholds   <- c(0.01, 0.025, 0.05, 0.1)
lfc_thresholds <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2)
compute_r2_sweep <- function(elkin_full, gusto_fetal, gusto_maternal,
                             ref_label, level_label) {
  
  results <- expand.grid(
    p_thr   = p_thresholds,
    lfc_thr = lfc_thresholds,
    stringsAsFactors = FALSE
  )
  
  out <- lapply(seq_len(nrow(results)), function(i) {
    p_thr   <- results$p_thr[i]
    lfc_thr <- results$lfc_thr[i]
    
    elkin_sig <- subset(elkin_full, P < p_thr & abs(log2FC) > lfc_thr)
    
    if (nrow(elkin_sig) < 5) {
      return(data.frame(
        p_thr = p_thr, lfc_thr = lfc_thr,
        Reference = ref_label, Level = level_label,
        Exposure = c("Fetal", "Maternal"),
        R2 = NA_real_, pval = NA_real_, Sign_Conc = NA_real_, n = 0L
      ))
    }
    
    mm1 <- merge(elkin_sig,
                 gusto_fetal[, c('Ensembl', 'log2FC', 'SE')],
                 by = 'Ensembl', suffixes = c('_Elkin', '_Fetal'))
    mm <- merge(mm1,
                gusto_maternal[, c('Ensembl', 'log2FC', 'SE')],
                by = 'Ensembl')
    colnames(mm)[colnames(mm) == 'log2FC'] <- 'log2FC_Maternal'
    colnames(mm)[colnames(mm) == 'SE']     <- 'SE_Maternal'
    
    if (nrow(mm) < 5) {
      return(data.frame(
        p_thr = p_thr, lfc_thr = lfc_thr,
        Reference = ref_label, Level = level_label,
        Exposure = c("Fetal", "Maternal"),
        R2 = NA_real_, pval = NA_real_, Sign_Conc = NA_real_, n = nrow(mm)
      ))
    }
    
    get_stats <- function(x, y, exposure) {
      valid <- complete.cases(x, y)
      if (sum(valid) < 5) return(data.frame(Exposure = exposure, R2 = NA_real_,
                                            pval = NA_real_, Sign_Conc = NA_real_,
                                            n = sum(valid)))
      ct <- cor.test(x[valid], y[valid])
      data.frame(Exposure  = exposure,
                 R2        = ct$estimate^2,
                 pval      = ct$p.value,
                 Sign_Conc = mean(sign(x[valid]) == sign(y[valid])),
                 n         = sum(valid))
    }
    
    bind_rows(
      get_stats(mm$log2FC_Elkin, mm$log2FC_Fetal, "Fetal"),
      get_stats(mm$log2FC_Elkin, mm$log2FC_Maternal, "Maternal")
    ) %>% mutate(p_thr = p_thr, lfc_thr = lfc_thr,
                 Reference = ref_label, Level = level_label)
  })
  
  bind_rows(out)
}

r2_df <- bind_rows(
  compute_r2_sweep(assembly.DGE.elkin, assembly.DGE.gusto.fetal,
                   assembly.DGE.gusto.maternal, "Assembly", "DGE"),
  compute_r2_sweep(assembly.DTE.elkin, assembly.DTE.gusto.fetal,
                   assembly.DTE.gusto.maternal, "Assembly", "DTE"),
  compute_r2_sweep(gencode.DGE.elkin, gencode.DGE.gusto.fetal,
                   gencode.DGE.gusto.maternal, "Gencode", "DGE"),
  compute_r2_sweep(gencode.DTE.elkin, gencode.DTE.gusto.fetal,
                   gencode.DTE.gusto.maternal, "Gencode", "DTE")
)

r2_df$LFC_label <- factor(paste0("|log2FC| > ", r2_df$lfc_thr),
                          levels = paste0("|log2FC| > ", lfc_thresholds))
r2_df$P_label   <- factor(paste0("P < ", r2_df$p_thr),
                          levels = paste0("P < ", p_thresholds))

# --- R² vs LFC threshold ---

r2_df <- r2_df %>%
  mutate(Reference = case_when(
    Reference == "Assembly" ~ "lr-assembly",
    Reference == "Gencode" ~ "gencode",
    TRUE ~ Reference
  ))

p2 <- ggplot(r2_df, aes(x = lfc_thr, y = R2,
                  color = Reference, linetype = Exposure)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(Level ~ P_label, scales = "free_y") +
  labs(title = "R² between explant and tissue effect sizes",
       x = "Explant |log2FC| threshold", y = "R²") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom", legend.box = "vertical") +
  scale_color_manual(values=c("lr-assembly" = "#C4663E", "gencode" = "grey50"),name = "Annotation")

p2

ggsave("PFBS_explants_vs_tissue_sensitivity.pdf",p2,width=6.75,height=4.75,units="in",dpi=300,bg="white")




