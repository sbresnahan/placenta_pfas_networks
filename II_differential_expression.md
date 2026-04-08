---
output:
  html_document: default
  pdf_document: default
---

# Natural variation in transplacental transfer efficiency exposes distinct transcriptional network architectures of PFAS effects on birth weight and gestational age

## II. Differential expression analyses

**Overview**: Differential transcript and gene expression analyses with PFAS exposures in GUSTO placental tissue samples and patient-derived placental explants. All analyses use a RUVr + DESeq2 pipeline: residuals from an edgeR GLM fit are passed to RUVSeq's RUVr to estimate factors of unwanted variation, which are then included as covariates in a DESeq2 negative binomial model.

**Note**: The placenta reference transcriptome (lr-assembly) classification file used for transcript-to-gene mapping can be downloaded from our GitHub repository, [lr-placenta-transcriptomics/Assembly/ESPRESSO_corrected_SC_filtered_classification.txt.gz](https://github.com/sbresnahan/lr-placenta-transcriptomics/blob/main/Assembly/ESPRESSO_corrected_SC_filtered_classification.txt.gz)

### GUSTO analyses

**Purpose**: Test for differential transcript and gene expression associated with p=8 PFAS compounds measured in maternal blood and cord blood at delivery (16 exposures total) in term placental samples from the GUSTO birth cohort (N=124). Each PFAS is modeled as a continuous exposure in separate DESeq2 models adjusting for fetal sex and gestational age, with RUVr-estimated latent factors included to correct for unwanted technical variation.

#### GUSTO differential transcript expression (lr-assembly)

📜 [PFAS_DETx.R](II_differential_expression/PFAS_DETx.R)

This R script loads Salmon quantification estimates (corrected for overdispersion from read-to-transcript mapping ambiguity following Baldoni et al.) and PFAS exposure data, merges covariates into the SummarizedExperiment objects, and runs the RUVr + DESeq2 pipeline independently for each of 16 PFAS exposures (8 maternal, 8 cord blood) at the transcript level. Results, volcano plots, and CSV exports are generated per exposure.

#### GUSTO differential gene expression (lr-assembly)

📜 [PFAS_DEG.R](II_differential_expression/PFAS_DEG.R)

This R script follows the same data loading and covariate merging workflow as the transcript-level analysis, but operates on gene-level SummarizedExperiment objects (`gse_gusto.tx`). The RUVr + DESeq2 pipeline is run independently for each of 16 PFAS exposures adjusting for fetal sex and gestational age. No isoform aggregation is performed since the analysis is already at the gene level.

#### GUSTO differential transcript expression (GENCODE v45)

📜 [PFAS_DETx_gencode.R](II_differential_expression/PFAS_DETx_gencode.R)

This R script repeats the transcript-level differential expression analysis using GENCODE v45 quantification estimates instead of the lr-assembly, enabling direct comparison of annotation-dependent results. The same RUVr + DESeq2 pipeline, and covariates are applied.

#### GUSTO differential gene expression (GENCODE v45)

📜 [PFAS_DEG_gencode.R](II_differential_expression/PFAS_DEG_gencode.R)

This R script repeats the gene-level differential expression analysis using GENCODE v45 quantification estimates. The same RUVr + DESeq2 pipeline and covariates are applied as in the lr-assembly gene-level analysis.

### Explant analyses

**Purpose**: Test for differential transcript and gene expression associated with PFBS dose in patient-derived placental explants (N=18 tissue samples from 3 third-trimester and 3 second-trimester placentas) treated at 0 µM, 5 µM, and 20 µM for 24 hours. PFBS concentration is modeled as a scaled continuous variable in DESeq2 models adjusting for trimester (all-tissue analysis) or replicate (term-only analysis), with RUVr correction for unwanted variation. Analyses are performed for both all-tissue and term-only subsets at both transcript and gene levels.

#### Explant differential transcript and gene expression (lr-assembly)

📜 [Explant_DTE_DGE.R](II_differential_expression/Explant_DTE_DGE.R)

This R script imports Salmon quantification estimates against the lr-assembly for explant samples, constructs SummarizedExperiment objects at transcript and gene levels, and runs the RUVr + DESeq2 pipeline for PFBS dose-response. Transcript-level DTE is run on all tissue samples (adjusting for trimester, k=9 RUVr factors) and on the term-only subset (adjusting for replicate, k=2). Gene-level DGE follows the same design (all tissue: k=10; term-only: k=3).

#### Explant differential transcript and gene expression (GENCODE v45)

📜 [Explant_DTE_DGE_gencode.R](II_differential_expression/Explant_DTE_DGE_gencode.R)

This R script repeats the explant DTE and DGE analyses using GENCODE v45 quantification estimates instead of the lr-assembly, with the same RUVr + DESeq2 pipeline, dose-response design, and analytical subsets.

### Summary and cross-cohort validation

#### Summarize GUSTO differential expression results

📜 [summarize_DTE_DGE.R](II_differential_expression/summarize_DTE_DGE.R)

This R script summarizes the GUSTO differential expression results across all 8 PFAS compounds for both maternal and fetal exposure windows. For each PFAS, it counts the number of significant effects (FDR < 10%, |log2FC| > 1) at both the transcript level (genes with isoform-level effects, mapped via the lr-assembly transcript-to-gene table) and gene level (genes with gene-level effects), generating grouped bar plots by exposure window. A combined figure faceted by maternal vs. fetal exposure is produced with PFAS ordered by transplacental transfer efficiency (TPTE). The script also computes Pearson correlations between maternal and fetal log2FC effect estimates for each PFAS at both the DTE and DGE levels, and regresses these correlations against TPTE to test whether maternal–fetal concordance of transcriptional effects decreases with increasing transplacental transfer.

#### GUSTO–Explant overlap analysis

📜 [GUSTO_Explant_overlap.R](II_differential_expression/GUSTO_Explant_overlap.R)

This R script quantifies the concordance of PFBS differential expression results between the placental explant (experimental) and GUSTO tissue (observational) cohorts, separately for lr-assembly and GENCODE v45 annotations at both transcript and gene levels. Direction-aware overlap is computed by intersecting features that are significant and concordant in sign between datasets. Fisher's exact tests assess enrichment of concordant effects relative to the annotation-specific background, comparing explant vs. tissue (maternal), explant vs. tissue (fetal), and tissue maternal vs. fetal contrasts. The script also includes a sensitivity analysis that sweeps across p-value (0.01–0.1) and |log2FC| (0–2) thresholds, computing R² between explant and tissue effect sizes at each combination to evaluate robustness of cross-cohort concordance across both annotations and exposure windows.
