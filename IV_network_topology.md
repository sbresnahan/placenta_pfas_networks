---
output:
  html_document: default
  pdf_document: default
---

# Natural variation in transplacental transfer efficiency exposes distinct transcriptional network architectures of PFAS effects on birth weight and gestational age

## IV. Network topology analyses

**Overview**: This section describes analyses of co-expression network topology to characterize the structural position of PFAS-outcome mediator transcripts within WGCNA modules. Analyses include computation of hub-to-mediator graph distances across the 144-parameterization bootstrap grid, extraction of peak distances from hub-distance density distributions, comparison of mediation effect sizes between module hubs and differentially expressed transcripts/genes, and visualization of module subnetworks. All network analyses were performed independently for birth weight and gestational age outcomes, and at both the transcript and gene level.

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Primary working directories
DIR_WGCNA=/path/to/WGCNA_analysis              # WGCNA bootstrap outputs and adjacency matrices
DIR_DISTANCE=/path/to/distance_analysis         # Hub-to-mediator distance computation

# Input data files
# adjacency.rds / adjacency_GA.rds:                Signed weighted adjacency matrices (power=3)
# data_for_distance_analysis.RDS / _GA.RDS:        all_module_mediation_df, pfas_order, tpte_values
# WGCNA_output_<iteration>.RData:                  Per-iteration WGCNA results from bootstrap grid
```

### Hub-to-mediator graph distances: Birth weight

**Purpose**: For each WGCNA parameterization, construct an igraph network from the signed adjacency matrix within each module, identify the eigenvector centrality hub, and compute weighted shortest-path distances from the hub to every significant mediator transcript. Distances are computed per module × PFAS × exposure window combination and saved per iteration for downstream aggregation.

📜 [run_distance_analysis_BW.R](IV_network_topology/run_distance_analysis_BW.R)

This R script accepts an iteration number, loads the birth weight adjacency matrix and mediation results for that parameterization, and computes hub-to-mediator distances within each module. For each module, edges below 0.1 are thresholded to zero, a weighted undirected graph is built, the eigenvector centrality hub is identified, and shortest-path distances (weighted by inverse edge weight) from the hub to all significant mediators are computed for each PFAS × exposure window. Per-transcript distance results are saved to `transcript_distances/transcript_distances_<iteration>.rds`.

```bash
################################################################################
# Submit Hub-to-Mediator Distance Array Job — Birth Weight
################################################################################

# Number of iterations from parameter grid
N_ITERATIONS=144  # 3 powers × 4 minModuleSizes × 3 deepSplits × 4 mergeCutHeights

# Submit array job
bsub -J "distance_analysis[1-${N_ITERATIONS}]" \
     -W 3:00 \
     -n 4 \
     -M 120 \
     -R "rusage[mem=120]" \
     -o "${DIR_WGCNA}/logs/distance_%J_%I.out" \
     -e "${DIR_WGCNA}/logs/distance_%J_%I.err" \
     -q short \
     "cd ${DIR_WGCNA} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript run_distance_analysis_BW.R \${LSB_JOBINDEX}"

# -M 120 / rusage[mem=120]: 120 GB memory (full adjacency matrix held in memory)
```

### Hub-to-mediator graph distances: Gestational age

**Purpose**: Repeat hub-to-mediator distance computation for the gestational age outcome, using the GA-specific adjacency matrix and mediation results.

📜 [run_distance_analysis_GA.R](IV_network_topology/run_distance_analysis_GA.R)

This R script follows the same structure as the birth weight version, but loads the gestational age adjacency matrix (`adjacency_GA.rds`) and distance analysis inputs (`data_for_distance_analysis_GA.RDS`). Results are saved to `transcript_distances_GA/transcript_distances_<iteration>.rds`.

```bash
################################################################################
# Submit Hub-to-Mediator Distance Array Job — Gestational Age
################################################################################

N_ITERATIONS=144

bsub -J "distance_analysis[1-${N_ITERATIONS}]" \
     -W 3:00 \
     -n 4 \
     -M 120 \
     -R "rusage[mem=120]" \
     -o "${DIR_WGCNA}/logs/distance_%J_%I.out" \
     -e "${DIR_WGCNA}/logs/distance_%J_%I.err" \
     -q short \
     "cd ${DIR_WGCNA} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript run_distance_analysis_GA.R \${LSB_JOBINDEX}"
```

---

### Peak hub-distance extraction: Birth weight

**Purpose**: Summarize the per-transcript distance distributions from the previous step by extracting the mode (peak of kernel density estimate) of the hub-to-mediator distance distribution for each module × PFAS × exposure window, providing a single measure of how centrally or peripherally significant mediators are positioned relative to the module hub.

📜 [run_peak_distance_BW.R](IV_network_topology/run_peak_distance_BW.R)

This R script accepts an iteration number and, for each module, reconstructs the weighted igraph network, identifies the eigenvector centrality hub, computes shortest-path distances from the hub to all significant mediators, and extracts the peak of the distance density (or mean if fewer than 2 mediators). Results are merged with TPTE values and saved to `peak_distances/transcript_peak_distances_<iteration>.rds`.

```bash
################################################################################
# Submit Peak Distance Extraction Array Job — Birth Weight
################################################################################

N_ITERATIONS=144

bsub -J "distance_analysis[1-${N_ITERATIONS}]" \
     -W 3:00 \
     -n 4 \
     -M 120 \
     -R "rusage[mem=120]" \
     -o "${DIR_WGCNA}/logs/distance_%J_%I.out" \
     -e "${DIR_WGCNA}/logs/distance_%J_%I.err" \
     -q short \
     "cd ${DIR_WGCNA} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript run_peak_distance_BW.R \${LSB_JOBINDEX}"
```

### Peak hub-distance extraction: Gestational age

**Purpose**: Repeat peak hub-distance extraction for gestational age, using the GA-specific adjacency matrix and mediation results.

📜 [run_peak_distance_GA.R](IV_network_topology/run_peak_distance_GA.R)

This R script follows the same structure as the birth weight version, loading `adjacency_GA.rds` and `data_for_distance_analysis_GA.RDS`. Results are saved to `peak_distances_GA/transcript_peak_distances_<iteration>.rds`.

```bash
################################################################################
# Submit Peak Distance Extraction Array Job — Gestational Age
################################################################################

N_ITERATIONS=144

bsub -J "distance_analysis[1-${N_ITERATIONS}]" \
     -W 3:00 \
     -n 4 \
     -M 120 \
     -R "rusage[mem=120]" \
     -o "${DIR_WGCNA}/logs/distance_%J_%I.out" \
     -e "${DIR_WGCNA}/logs/distance_%J_%I.err" \
     -q short \
     "cd ${DIR_WGCNA} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript run_peak_distance_GA.R \${LSB_JOBINDEX}"
```

---

### Mediation effect sizes, co-expression topology, and TPTE: Birth weight (transcript level)

**Purpose**: Integrate mediation results with module membership (kME) and network topology to characterize the structural position of PFAS-birth weight mediator transcripts. This includes merging bootstrap WGCNA results with mediation outputs, comparing absolute ACME between differentially expressed transcripts and module hubs, computing ADE forest plots by PFAS and exposure source, modeling the relationship between TPTE and co-expression strength (peak |kME|) and mediator counts, and testing module-level enrichment for significant mediators via Fisher's exact test.

📜 [acme_vs_kME_BW.R](IV_network_topology/acme_vs_kME_BW.R)

This R script loads WGCNA outputs, bootstrap iteration results, and mediation batch results for birth weight. It computes kME values across all 144 parameterizations, merges with ACME results, and performs the following analyses: (1) bootstrap ANOVA comparing |ACME| between DETs, non-mediator hubs, and mediator hubs; (2) inverse-variance weighted ADE forest plots by PFAS and exposure source; (3) density plots of mediator |kME| distributions faceted by PFAS (ordered by TPTE); (4) linear regression of TPTE vs. peak |kME| and mediator count; (5) module enrichment odds ratios (Fisher's exact test) for significant mediators by PFAS; and (6) scatterplots of average module enrichment proportion vs. peak |kME|.

### Mediation effect sizes, co-expression topology, and TPTE: Birth weight (gene level)

📜 [acme_vs_kME_BW_gene.R](IV_network_topology/acme_vs_kME_BW_gene.R)

This R script repeats the birth weight mediation–topology integration at the gene level, loading gene-level WGCNA outputs (`input_for_WGCNA_gene.RData`), gene-level bootstrap results (`WGCNA_output_<iteration>_gene.RData`), and gene-level mediation batches. The same analyses are performed: |ACME| comparison between DEGs and module hubs, ADE forest plots, kME density distributions, TPTE regression models, and module enrichment testing.

### Mediation effect sizes, co-expression topology, and TPTE: Gestational age

📜 [acme_vs_kME_GA.R](IV_network_topology/acme_vs_kME_GA.R)

This R script repeats the mediation–topology integration for gestational age at the transcript level. It loads GA-specific expression and covariate data (filtered to spontaneous labor onset), merges GA mediation batch results, and performs the same suite of analyses: |ACME| hub vs. DET comparison, ADE forest plots, kME density distributions, TPTE vs. peak |kME| and mediator count regression, and module enrichment testing.

---

### Network visualization: Birth weight (transcript level)

📜 [netViz_BW.R](IV_network_topology/netViz_BW.R)

This R script generates module subnetwork visualizations for the birth weight transcript-level analysis. For each module, it constructs an igraph network from the signed adjacency matrix, identifies the eigenvector centrality hub, computes hub-to-node distances, and renders force-directed (Fruchterman-Reingold) layouts using ggraph with nodes colored by mediation significance and sized by |kME|, and edges weighted by adjacency strength. Accompanying distance histograms show the distribution of hub-to-mediator distances.

### Network visualization: Birth weight (gene level)

📜 [netViz_BW_gene.R](IV_network_topology/netViz_BW_gene.R)

This R script repeats the module subnetwork visualizations at the gene level, loading the gene-level adjacency matrix (`adjacency_gene.RDS`) and gene-level mediation workspace.

### Network visualization: Gestational age

📜 [netViz_GA.R](IV_network_topology/netViz_GA.R)

This R script repeats the module subnetwork visualizations for the gestational age outcome, loading the GA-specific adjacency matrix (`adjacency_GA.rds`) and GA mediation workspace (filtered to spontaneous labor onset). Includes the same force-directed network layouts and hub-distance density plots.

---
