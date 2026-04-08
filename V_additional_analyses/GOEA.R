# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(stringr)

# tx2g mapping table
tx2g <- readRDS("tx2g_assembly.RDS")

# Read tx mediation results
mediation_results <- readRDS("all_module_mediation_df_tx.RDS")

# Filter for significant strict mediators
significant_mediators <- mediation_results %>%
  filter(significant_strict == TRUE, 
         exposure_window == "mat")

# Get unique transcript IDs
unique_transcripts <- unique(significant_mediators$transcript_id)

# Join with tx2g to get gene IDs
transcript_to_gene <- significant_mediators %>%
  select(transcript_id) %>%
  distinct() %>%
  left_join(tx2g, by = "transcript_id")

# Remove version numbers from gene_id to get Ensembl IDs
gene_list <- transcript_to_gene %>%
  mutate(gene_id_clean = str_remove(gene_id, "\\..*$")) %>%
  pull(gene_id_clean) %>%
  unique()

# Convert Ensembl to Entrez IDs (required for KEGG)
gene_mapping <- bitr(gene_list, 
                     fromType = "ENSEMBL",
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = org.Hs.eg.db)

entrez_ids <- gene_mapping$ENTREZID

# ============================================================================
# GO ENRICHMENT ANALYSIS
# ============================================================================

# Biological Process
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# Molecular Function
go_mf <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# Cellular Component
go_cc <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# ============================================================================
# KEGG ENRICHMENT ANALYSIS
# ============================================================================

kegg <- enrichKEGG(gene = entrez_ids,
                   organism = 'hsa',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2)

# Convert KEGG gene IDs to symbols for readability
kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# ============================================================================
# SIMPLIFY AND SUMMARIZE RESULTS (manuscript-ready)
# ============================================================================

# Function to simplify GO results by removing redundant terms
simplify_go <- function(go_result) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    return(NULL)
  }
  simplified <- clusterProfiler::simplify(go_result)
  return(simplified)
}


# Simplify GO results (removes redundant terms based on semantic similarity)
go_bp_simple <- simplify_go(go_bp)
go_mf_simple <- simplify_go(go_mf)
go_cc_simple <- simplify_go(go_cc)

# ============================================================================
# MANUSCRIPT-READY SUMMARY TABLES
# ============================================================================

# Function to create manuscript table
create_manuscript_table <- function(enrichment_result, top_n = 10, type = "GO") {
  if (is.null(enrichment_result) || nrow(enrichment_result@result) == 0) {
    return(data.frame(Message = "No significant enrichment"))
  }
  
  df <- as.data.frame(enrichment_result)
  
  # Select top N by adjusted p-value
  df <- df %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)
  
  # Add readable fold-enrichment
  df <- df %>%
    mutate(
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
      BgRatio_num = sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2])),
      FoldEnrichment = round(GeneRatio_num / BgRatio_num, 2)
    ) %>%
    select(-GeneRatio_num, -BgRatio_num)
  
  return(df)
}

# Generate tables
table_go_bp <- create_manuscript_table(go_bp_simple, top_n = 10, type = "GO-BP")
table_go_mf <- create_manuscript_table(go_mf_simple, top_n = 10, type = "GO-MF")
table_go_cc <- create_manuscript_table(go_cc_simple, top_n = 10, type = "GO-CC")
table_kegg <- create_manuscript_table(kegg, top_n = 10, type = "KEGG")

# ============================================================================
# GENERATE MANUSCRIPT TEXT SUMMARY
# ============================================================================

generate_text_summary <- function(go_bp, go_mf, go_cc, kegg, n_terms = 5) {
  
  summary_text <- "Functional enrichment analysis revealed:\n\n"
  
  # GO Biological Process
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    top_bp <- head(go_bp@result$Description, n_terms)
    summary_text <- paste0(summary_text, 
                           "**Gene Ontology - Biological Process:** ",
                           "Enriched terms included ",
                           paste(paste0("'", top_bp, "'"), collapse = ", "),
                           " (FDR < 0.2).\n\n")
  }
  
  # GO Molecular Function
  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    top_mf <- head(go_mf@result$Description, n_terms)
    summary_text <- paste0(summary_text,
                           "**Gene Ontology - Molecular Function:** ",
                           "Enriched functions included ",
                           paste(paste0("'", top_mf, "'"), collapse = ", "),
                           " (FDR < 0.2).\n\n")
  }
  
  # GO Cellular Component
  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    top_cc <- head(go_cc@result$Description, n_terms)
    summary_text <- paste0(summary_text,
                           "**Gene Ontology - Cellular Component:** ",
                           "Enriched components included ",
                           paste(paste0("'", top_cc, "'"), collapse = ", "),
                           " (FDR < 0.2).\n\n")
  }
  
  # KEGG Pathways
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    top_kegg <- head(kegg@result$Description, n_terms)
    summary_text <- paste0(summary_text,
                           "**KEGG Pathways:** ",
                           "Enriched pathways included ",
                           paste(paste0("'", top_kegg, "'"), collapse = ", "),
                           " (FDR < 0.2).")
  }
  
  return(summary_text)
}

# Generate text
manuscript_summary <- generate_text_summary(go_bp_simple, go_mf_simple, go_cc_simple, kegg, n_terms = 5)
cat(manuscript_summary)



# ============================================================================
# EXPORT RESULTS
# ============================================================================

# Save tables
write.csv(table_go_bp, "GO_BiologicalProcess_top10.csv", row.names = FALSE)
write.csv(table_go_mf, "GO_MolecularFunction_top10.csv", row.names = FALSE)
write.csv(table_go_cc, "GO_CellularComponent_top10.csv", row.names = FALSE)
write.csv(table_kegg, "KEGG_Pathways_top10.csv", row.names = FALSE)

# Save full results
write.csv(as.data.frame(go_bp_simple), "GO_BP_full.csv", row.names = FALSE)
write.csv(as.data.frame(go_mf_simple), "GO_MF_full.csv", row.names = FALSE)
write.csv(as.data.frame(go_cc_simple), "GO_CC_full.csv", row.names = FALSE)
write.csv(as.data.frame(kegg), "KEGG_full.csv", row.names = FALSE)

# ============================================================================
# VISUALIZATIONS (optional - for supplement)
# ============================================================================

# Dotplot for GO BP
if (!is.null(go_bp_simple) && nrow(go_bp_simple@result) > 0) {
  pdf("GO_BP_dotplot.pdf", width = 10, height = 8)
  print(dotplot(go_bp_simple, showCategory = 20) + 
          ggtitle("GO Biological Process Enrichment"))
  dev.off()
}

# Dotplot for KEGG
if (!is.null(kegg) && nrow(kegg@result) > 0) {
  pdf("KEGG_dotplot.pdf", width = 10, height = 8)
  print(dotplot(kegg, showCategory = 20) + 
          ggtitle("KEGG Pathway Enrichment"))
  dev.off()
}

# Combined plot
if (!is.null(go_bp_simple) && nrow(go_bp_simple@result) > 0) {
  pdf("enrichment_barplot.pdf", width = 12, height = 6)
  print(barplot(go_bp_simple, showCategory = 15))
  dev.off()
}

# ============================================================================
# PRINT SUMMARY
# ============================================================================

cat("\n=== ENRICHMENT ANALYSIS SUMMARY ===\n")
cat("GO BP terms found:", ifelse(is.null(go_bp_simple), 0, nrow(go_bp_simple@result)), "\n")
cat("GO MF terms found:", ifelse(is.null(go_mf_simple), 0, nrow(go_mf_simple@result)), "\n")
cat("GO CC terms found:", ifelse(is.null(go_cc_simple), 0, nrow(go_cc_simple@result)), "\n")
cat("KEGG pathways found:", ifelse(is.null(kegg), 0, nrow(kegg@result)), "\n")
cat("\nManuscript-ready text summary:\n")
cat(manuscript_summary)