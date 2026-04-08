---
output:
  html_document: default
  pdf_document: default
---

# Natural variation in transplacental transfer efficiency exposes distinct transcriptional network architectures of PFAS effects on birth weight and gestational age

## V. Additional analyses

**Overview**: This section describes supplementary analyses including Gene Ontology and KEGG pathway enrichment of significant mediator transcripts, maternal–fetal PFAS correlation characterization, summary forest plots of network topology metric slopes against TPTE, cross-outcome mediator overlap quantification, and the birth weight–gestational age correlation in the spontaneous labor subset.

### Gene Ontology and KEGG pathway enrichment

📜 [GOEA.R](V_additional_analyses/GOEA.R)

This R script performs functional enrichment analysis on the set of genes corresponding to significant mediator transcripts (ACME p < 0.05, confidence intervals excluding zero) identified in the maternal exposure window. Transcript IDs are mapped to Ensembl gene IDs via the lr-assembly transcript-to-gene table, then converted to Entrez IDs for enrichment testing. GO enrichment is run separately for Biological Process, Molecular Function, and Cellular Component ontologies using `clusterProfiler::enrichGO` with BH-adjusted p-value < 0.05 and q-value < 0.2 thresholds. KEGG pathway enrichment is performed via `enrichKEGG`. GO results are simplified by semantic similarity to remove redundant terms. The script exports top-10 manuscript-ready summary tables, full result tables, dotplots, and a text summary for each enrichment category.

### Supplementary visualizations and cross-outcome comparisons

📜 [extra.R](V_additional_analyses/extra.R)

This R script generates several supplementary analyses and figures. First, it characterizes maternal–fetal PFAS correlations by computing Pearson correlations and scatterplots for all 8 PFAS compounds between maternal and cord blood z-scores. Second, it assembles a summary forest plot of TPTE regression slopes across four network topology metrics (mediator count, co-expression strength, centrality, and compartmentalization) for both birth weight and gestational age outcomes, loading pre-computed slope estimates from the section IV analyses. Third, it merges birth weight and gestational age mediation results, quantifies cross-outcome mediator overlap using UpSet plots, and identifies transcripts that are significant mediators across all exposure window × outcome combinations. Fourth, it computes the gene-level Jaccard index and Fisher's exact test log odds ratios for maternal–fetal mediator gene sharing as a function of TPTE, separately for each outcome. Finally, it visualizes the birth weight–gestational age correlation in the spontaneous labor subset.

---
