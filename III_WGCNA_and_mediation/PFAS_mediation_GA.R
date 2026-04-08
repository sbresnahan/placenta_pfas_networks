#!/usr/bin/env Rscript

# Set library path for HPC
.libPaths(c("/home/stbresnahan/R/ubuntu/4.3.1", .libPaths()))

# Get batch number from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No batch number provided. Usage: Rscript PFAS_mediation_GA.R <batch_number>")
}

ITERATION <- as.integer(args[1])
cat("Starting batch:", ITERATION, "\n")
cat("Started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Load required libraries
library(dplyr)
library(mediation)
library(parallel)

# Set working directory
setwd("/rsrch5/home/epi/stbresnahan/scratch/PFAS_mediation")

# Load necessary data
cat("Loading data...\n")
all_module_transcripts_df <- readRDS("all_module_transcripts_df.rds")
expr_PFAS <- readRDS("expr_PFAS_GA.rds")
covs_PFAS <- readRDS("covs_PFAS_GA.rds")

#-------------------------------
# Mediation function
# Outcome: Gestational Age (GA)
#-------------------------------
perform_mediation_all <- function(transcript_id,
                                  pfas_column,
                                  pfas,
                                  exposure_window,
                                  transcript_type,
                                  module = NA) {
  
  # Get expression data
  if (!transcript_id %in% rownames(expr_PFAS)) {
    return(NULL)
  }
  mediator <- expr_PFAS[transcript_id, ]
  
  # Get PFAS exposure
  if (!pfas_column %in% colnames(covs_PFAS)) {
    return(NULL)
  }
  exposure <- covs_PFAS[[pfas_column]]
  
  # Outcome: gestational age
  outcome <- covs_PFAS[["GA"]]
  
  # Covariates
  fetal_sex <- covs_PFAS[["sex"]]
  
  # Length checks
  n_samples <- length(exposure)
  if (length(mediator) != n_samples ||
      length(outcome) != n_samples ||
      length(fetal_sex) != n_samples) {
    return(NULL)
  }
  
  # Assemble analysis frame
  med_data <- data.frame(
    exposure = exposure,
    mediator = mediator,
    outcome = outcome,
    fetal_sex = fetal_sex
  ) %>%
    filter(
      !is.na(exposure) &
        !is.na(mediator) &
        !is.na(outcome) &
        !is.na(fetal_sex)
    )
  
  if (nrow(med_data) < 20) {
    return(NULL)
  }
  
  # Fit models
  tryCatch({
    med_fit <- lm(mediator ~ exposure + fetal_sex, data = med_data)
    out_fit <- lm(outcome ~ exposure + mediator + fetal_sex, data = med_data)
    
    med_result <- mediate(
      med_fit,
      out_fit,
      treat = "exposure",
      mediator = "mediator",
      boot = TRUE,
      sims = 100
    )
    
    data.frame(
      transcript_id = transcript_id,
      PFAS = pfas,
      exposure_window = exposure_window,
      transcript_type = transcript_type,
      module = module,
      ACME = med_result$d0,
      ACME_CI_lower = med_result$d0.ci[1],
      ACME_CI_upper = med_result$d0.ci[2],
      ACME_p = med_result$d0.p,
      ADE = med_result$z0,
      ADE_CI_lower = med_result$z0.ci[1],
      ADE_CI_upper = med_result$z0.ci[2],
      ADE_p = med_result$z0.p,
      prop_mediated = med_result$n0,
      n = nrow(med_data)
    )
  }, error = function(e) {
    return(NULL)
  })
}

#--------------------------------
# Filter to current batch
#--------------------------------
batch_data <- all_module_transcripts_df %>%
  filter(batch_id == ITERATION)

cat("Batch", ITERATION, "contains", nrow(batch_data), "rows\n")

if (nrow(batch_data) == 0) {
  cat("WARNING: Batch", ITERATION, "is empty!\n")
  quit(save = "no", status = 1)
}

# Parallel settings
n_cores <- 12
cat("Running mediation analyses in parallel with", n_cores, "cores\n")

results_list <- mclapply(
  1:nrow(batch_data),
  function(i) {
    tryCatch({
      perform_mediation_all(
        batch_data$transcript_id[i],
        batch_data$PFAS_column[i],
        batch_data$PFAS[i],
        batch_data$exposure_window[i],
        "Module Transcript",
        batch_data$module[i]
      )
    }, error = function(e) NULL)
  },
  mc.cores = n_cores
)

# Combine results
results_list <- results_list[!sapply(results_list, is.null)]

#--------------------------------
# Save results
#--------------------------------
if (!dir.exists("Results_GA")) {
  dir.create("Results_GA", recursive = TRUE)
}

if (length(results_list) > 0) {
  batch_results <- bind_rows(results_list)
  output_file <- sprintf("Results_GA/mediation_batch_%04d.rds", ITERATION)
  saveRDS(batch_results, output_file)
  cat("Saved", nrow(batch_results), "results to", output_file, "\n")
} else {
  cat("WARNING: No results generated for batch", ITERATION, "\n")
}

cat("Batch", ITERATION, "complete at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
