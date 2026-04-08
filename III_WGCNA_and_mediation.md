---
output:
  html_document: default
  pdf_document: default
---

# Natural variation in transplacental transfer efficiency exposes distinct transcriptional network architectures of PFAS effects on birth weight and gestational age


## III. Weighted gene co-expression network analysis and causal mediation

**Overview**: This section describes the construction of signed co-expression networks using WGCNA, parameter sensitivity analysis via grid-search bootstrapping across 144 parameterizations, and transcript-level causal mediation analysis linking PFAS exposures to birth weight and gestational age through module member transcripts.

**Note**: All analyses described below were performed independently at both the transcript level and the gene level.

### Setup

```bash
################################################################################
# Directory Structure and File Paths
################################################################################

# Primary working directories
DIR_WGCNA=/path/to/WGCNA_analysis        # WGCNA network construction and bootstrapping
DIR_MEDIATION=/path/to/PFAS_mediation     # Causal mediation analysis

# Input data files (shared across analyses)
# input_for_WGCNA.RData:     expr_PFAS (expression matrix), covs_PFAS (covariate/trait data)
# all_module_transcripts_df:  Pre-computed transcript-module-PFAS mapping with batch assignments
```

### WGCNA: Single run and bootstrap setup

**Purpose**: Perform an initial WGCNA network construction to establish baseline module-trait relationships, then prepare a parameter grid for systematic sensitivity analysis across soft-thresholding powers, minimum module sizes, tree-cutting depths, and module merge thresholds.

📜 [WGCNA_single_run_and_bootstrap_setup.R](III_WGCNA_and_mediation/WGCNA_single_run_and_bootstrap_setup.R)

This R script defines core functions for signed WGCNA network construction (`run_optimized_wgcna`), module-trait correlation with BH-adjusted p-values (`correlate_modules_traits`), extraction of significant associations (`extract_significant_associations`), and numeric-to-color module mapping (`create_numeric_module_mapping`). It executes a single baseline run with soft-thresholding power=3 and minModuleSize=50, computes module-trait correlations for p=8 PFAS measures across maternal and cord blood exposure windows, and generates a 144-combination parameter grid (3 powers × 4 minModuleSizes × 3 deepSplits × 4 mergeCutHeights) saved for downstream bootstrapping.

### WGCNA: Bootstrap parameter sensitivity analysis

**Purpose**: Systematically evaluate WGCNA network stability across 144 parameter combinations via an LSF array job, testing sensitivity of module detection and module-trait associations to analytical choices.

📜 [WGCNA_bootstrap.R](III_WGCNA_and_mediation/WGCNA_bootstrap.R)

This R script accepts a single command-line argument (iteration number) indexing into the pre-computed parameter grid. For each parameterization, it constructs a signed co-expression network using the step-by-step approach (adjacency → TOM → hierarchical clustering → dynamic tree cutting → module merging) or blockwise approach for large datasets, computes module-trait correlations, and saves per-iteration results (`WGCNA_output_<iteration>.RData`) containing module assignments, eigengenes, correlation matrices, and significant associations.

```bash
################################################################################
# Submit WGCNA Bootstrap Array Job
################################################################################

# Number of iterations from parameter grid
# 3 powers × 4 minModuleSizes × 3 deepSplits × 4 mergeCutHeights
N_ITERATIONS=144

# Submit array job
bsub -J "distance_analysis[1-${N_ITERATIONS}]" \
     -W 3:00 \
     -n 8 \
     -M 64 \
     -R "rusage[mem=64]" \
     -o "${DIR_WGCNA}/logs/distance_%J_%I.out" \
     -e "${DIR_WGCNA}/logs/distance_%J_%I.err" \
     -q short \
     "cd ${DIR_WGCNA} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript WGCNA_bootstrap.R \${LSB_JOBINDEX}"

# Parameter explanations:
# -J: Array job with indices 1 through N_ITERATIONS
# -W 3:00: 3-hour wall time per iteration
# -n 8: 8 CPU cores
# -M 64 / rusage[mem=64]: 64 GB memory (TOM computation is memory-intensive)
# -q short: Short queue
# LSB_JOBINDEX: Passed as iteration number to index into parameter grid
```

---

### Mediation analysis: Birth weight

**Purpose**: Test whether WGCNA module member transcripts mediate the effect of PFAS exposures on birth weight using bootstrap-based causal mediation analysis, parallelized across 497 batches (~1,000 tests each; 233,152 total).

📜 [PFAS_mediation_BW.R](III_WGCNA_and_mediation/PFAS_mediation_BW.R)

This R script accepts a batch number, loads pre-assigned transcript-module-PFAS mappings (`all_module_transcripts_df`), and runs causal mediation analysis for each transcript in the batch. For each test, the mediator model regresses transcript expression on PFAS exposure adjusting for fetal sex and gestational age, and the outcome model regresses birth weight on PFAS exposure and transcript expression with the same covariates. Mediation is estimated via the `mediation` package with 100 bootstrap resamples. Results (ACME, ADE, proportion mediated, confidence intervals, p-values) are saved per-batch to `Results/mediation_batch_<XXXX>.rds`.

```bash
################################################################################
# Submit Mediation Array Job — Birth Weight
################################################################################

# Create output directories
mkdir -p ${DIR_MEDIATION}/Logs
mkdir -p ${DIR_MEDIATION}/Results

# Number of batches
N_BATCHES=497

# Submit array job
bsub -J "mediation[1-${N_BATCHES}]" \
     -W 1:00 \
     -n 12 \
     -M 48 \
     -R "rusage[mem=48]" \
     -o "${DIR_MEDIATION}/Logs/%J_%I.out" \
     -e "${DIR_MEDIATION}/Logs/%J_%I.err" \
     -q short \
     "cd ${DIR_MEDIATION} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript PFAS_mediation_BW.R \${LSB_JOBINDEX}"

# Each task processes ~1,000 mediation tests
# Total: 233,152 tests across 497 batches
# -n 12: 12 cores for mclapply parallelism within each batch
```

### Mediation analysis: Gestational age

**Purpose**: Repeat causal mediation analysis with gestational age as the outcome, using an adjusted covariate set (fetal sex only, since gestational age is now the outcome rather than a covariate).

📜 [PFAS_mediation_GA.R](III_WGCNA_and_mediation/PFAS_mediation_GA.R)

This R script follows the same structure as the birth weight mediation, but loads a separate expression matrix (`expr_PFAS_GA.rds`) and covariate file (`covs_PFAS_GA.rds`), regresses the outcome on gestational age instead of birth weight, and adjusts for fetal sex only (gestational age is removed from the covariate set since it is the outcome). Results are saved per-batch to `Results_GA/mediation_batch_<XXXX>.rds`.

```bash
################################################################################
# Submit Mediation Array Job — Gestational Age
################################################################################

# Create output directories
mkdir -p ${DIR_MEDIATION}/Logs
mkdir -p ${DIR_MEDIATION}/Results_GA

# Number of batches
N_BATCHES=497

# Submit array job
bsub -J "mediation[1-${N_BATCHES}]" \
     -W 1:00 \
     -n 12 \
     -M 48 \
     -R "rusage[mem=48]" \
     -o "${DIR_MEDIATION}/Logs/%J_%I.out" \
     -e "${DIR_MEDIATION}/Logs/%J_%I.err" \
     -q short \
     "cd ${DIR_MEDIATION} && \
      singularity exec --bind /rsrch5 \
      /risapps/singularity/repo/RStudio/4.3.1/rstudio_4.3.1.sif \
      Rscript PFAS_mediation_GA.R \${LSB_JOBINDEX}"

# Each task processes ~1,000 mediation tests
# Total: 233,152 tests across 497 batches
# -n 12: 12 cores for mclapply parallelism within each batch
```

---
