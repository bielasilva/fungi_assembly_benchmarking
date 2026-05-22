# ============================================================================
# 00_orchestrator.R
# ----------------------------------------------------------------------------
# Orchestrates all scripts, runs them and save the figures.
# ============================================================================

# Import libraries
pkgs <- c("tidyverse", "here")

lapply(pkgs, library, character.only=TRUE)

# ----------------------------------------------------------------------------
# Output directory
# ----------------------------------------------------------------------------

results_dir <- here::here("Empirical", "02_Processing", "R", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Run scripts
# ----------------------------------------------------------------------------

# Each script will save its results in the results_dir, so that they can be loaded by the next script without having to re-run the previous ones.

## 01_data_prep.R: Loads the simulated data, processes it and saves the processed data in results_dir.
source(here("Empirical", "02_Processing", "R", "01_data_prep.R"))
# load("results_dir/01_data_prep.RData")

## 02_fit_models.R: Fits the models and saves the results in results_dir.
source(here("Empirical", "02_Processing", "R", "02_fit_models.R"))
# load(file.path(results_dir, "02_fit_models.RData"))

## 03_figures.R: Creates the figures and saves them in results_dir/figures.
source(here("Empirical", "02_Processing", "R", "03_figures.R"))
# load(file.path(results_dir, "03_figures.RData"))

## 04_empirical_simulated.R: Runs the empirical vs simulated analysis.
source(here("Empirical", "02_Processing", "R", "04_empirical_simulated.R"))
