# ============================================================================
# 04_empirical_simulated.R
# ----------------------------------------------------------------------------
# Compares empirical and simulated results, using both raw means and model predictions.
#
# Provides:
#   - 
#   - 
#
# Sources:
#   - 
#   - 
# ============================================================================

# Clean environment and load libraries
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(lme4)
  library(here)
})

# ============================================================================
# 0. SETUP: COLORS, LABELS, AND FIGURE SAVING FUNCTION
# ============================================================================

colors_programs <- c("flye"           = "#2CA02C",
                      "spades_short"  = "#EDC948",
                      "abyss_short"   = "#FF7F0E"
                    )

program_labels <- c("flye"           ="Flye",
                    "spades_short"   ="SPAdes Short",
                    "abyss_short"    ="ABySS Short"
                    )


empirical_dir <- here("Empirical", "02_Processing", "R", "results")
simulated_dir  <- here("Simulated", "02_Processing", "R", "results")

figures_dir <- file.path(empirical_dir, "figures")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(plt, name , size_w, size_h) {
  formats <- c("png", "pdf", "tiff", "jpeg")
  
  for (ext in formats) {
    ggsave(filename = paste0(figures_dir, "/", name, ".", ext),
           plot = plt, device = ext,
           width = size_w, height = size_h,
           units = "in")
  }
}

# ============================================================================
# 1. LOAD AND NORMALIZE DATA
# ============================================================================

emp_stats <- read_csv(file.path(empirical_dir, "empirical_stats_single.csv"), show_col_types = FALSE) %>%
              filter(program %in% c("flye-ovl1000", "spades_short", "abyss_short")) %>%
              mutate(
                program = case_when(
                  program == "flye-ovl1000" ~ "flye",
                  TRUE ~ program
                ),
                dataset = "Empirical"
              )  %>%
              select(dataset, program, depth_eff, n,
                      mean_n50, CI_n50,
                      mean_busco, CI_busco,
                      mean_time_mb_p, CI_time_mb_p,
                      mean_ram_mb_p, CI_ram_mb_p,
                      mean_maxrss_gb, CI_maxrss_gb
                    ) %>%
              rename(mean_time_mb = mean_time_mb_p,
                      CI_time_mb   = CI_time_mb_p,
                      mean_ram_mb  = mean_ram_mb_p,
                      CI_ram_mb    = CI_ram_mb_p)


sim_stats <- read_csv(file.path(simulated_dir, "simulated_stats_single.csv"), show_col_types = FALSE) %>%
              filter(program %in% c("flye", "spades_short", "abyss_short")) %>%
              mutate(dataset = "Simulated") %>%
              select(dataset, program, depth_eff, n,
                      mean_n50, CI_n50,
                      mean_busco, CI_busco,
                      mean_time_mb, CI_time_mb,
                      mean_ram_mb, CI_ram_mb,
                      mean_maxrss_gb, CI_maxrss_gb)

combined <- bind_rows(sim_stats, emp_stats) %>%
  mutate(
    dataset = factor(dataset, levels = c("Simulated", "Empirical")),
    program = factor(program,
                     levels = c("flye", "spades_short", "abyss_short"))
  ) %>%
  mutate(mean_n50 = mean_n50 * 1e-6,  # convert from bp to Mbp for plotting on log scale
          CI_n50   = CI_n50 * 1e-6)

  # ============================================================================
# 2. MAIN FIGURE — RAW MEANS WITH 95% CI RIBBONS
# ============================================================================

# Helper that floors CI lower bound at a tiny positive value when on log scale.
ribbon_floor <- function(mean_val, ci_val, log_scale = FALSE, floor_val = 1) {
  ymin <- mean_val - ci_val
  if (log_scale) ymin <- pmax(ymin, floor_val)
  ymin
}

make_panel <- function(df, y_col, ci_col, y_label, y_log = FALSE, floor_val = 1) {
  df_plot <- df %>%
    mutate(
      .y    = .data[[y_col]],
      .ymin = ribbon_floor(.data[[y_col]], .data[[ci_col]], y_log, floor_val),
      .ymax = .data[[y_col]] + .data[[ci_col]]
    )

  p <- ggplot(df_plot, aes(x = depth_eff,
                           y = .y,
                           color = program,
                           linetype = dataset,
                           fill = program)) +
    geom_ribbon(aes(ymin = .ymin, ymax = .ymax),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(shape = dataset), size = 1.7) +
    scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 75, 100)) +
    scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
    scale_fill_manual(values=colors_programs, label=program_labels, name="Program") +
    scale_linetype_manual(values = c("Simulated" = "solid",
                                     "Empirical" = "dashed")) +
    scale_shape_manual(values = c("Simulated" = 16, "Empirical" = 17)) +
    labs(x = "Long-read depth (X)", y = y_label,
         color = "Assembler", fill = "Assembler",
         linetype = "Dataset", shape = "Dataset") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          plot.tag = element_text(face = "bold"))
  
  if (y_log) {
    p <- p + scale_y_log10(labels = label_number())
  }
  p
}

p_n50 <- make_panel(combined,
                    "mean_n50", "CI_n50",
                    "N50 (Mbp)",
                    y_log = TRUE, floor_val = 0.0005)

p_busco <- make_panel(combined,
                      "mean_busco", "CI_busco",
                      "BUSCO completeness (%)",
                      y_log = FALSE)

p_time <- make_panel(combined,
                     "mean_time_mb", "CI_time_mb",
                     "Runtime per Mbp (seconds)",
                     y_log = TRUE, floor_val = 0.01) +
                     coord_cartesian(ylim = c(5, NA))  # zoom in to better show differences at low depths

fig_main <- ggarrange(p_n50, p_busco, p_time, nrow = 1, common.legend = TRUE, legend = "bottom")

save_plot(fig_main, "Figure_X_sim_vs_empirical_main", 12, 5)

# ============================================================================
# 3. SUPPLEMENTARY FIGURE — MODEL PREDICTIONS WITH BOOTSTRAP CONFIDENCE INTERVALS
# ============================================================================

# Load the fitted models

readRDS(file.path(empirical_dir, "empirical_single_fits.rds")) -> empirical_single_fits
readRDS(file.path(simulated_dir, "simulated_single_fits.rds")) -> simulated_single_fits

# Define computational parameters for bootstrapping predictions
set.seed(123) # For reproducibility of bootstrapping results
n_cpus = 20 # Number of CPU cores to use for parallel processing
n_sim = 1000 # Number of bootstrap simulations

# Function to extract bootstrap predictions for one fitted model
get_model_predictions <- function(fit,
                                  programs_to_keep,
                                  depth_grid = seq(10, 100, by = 5),
                                  is_empirical = FALSE) {
  
  mod       <- fit$model
  response  <- fit$response
  transform <- fit$transform
  
  # Build prediction grid - covariates at their scaled mean (= 0)
  if (is_empirical) {
    pred_grid <- expand_grid(
      program          = programs_to_keep,
      depth_eff        = depth_grid
    ) %>%
      mutate(genome_size_est_s = 0)  # empirical has one covariate
  } else {
    pred_grid <- expand_grid(
      program   = programs_to_keep,
      depth_eff = depth_grid
    ) %>%
      mutate(
        size_mbp_s     = 0,
        gc_content_s   = 0,
        gene_density_s = 0,
        cds_per_gene_s = 0
      )
  }
  
  # Bootstrap predictions (parallelized)
  boot_obj <- bootMer(
    mod,
    nsim     = n_sim,
    parallel = "multicore",
    ncpus    = n_cpus,
    FUN = function(x) {
      predict(x, newdata = pred_grid, re.form = NA, allow.new.levels = TRUE)
    }
  )

  pred_grid$pred_t <- apply(boot_obj$t, 2, median, na.rm = TRUE)
  pred_grid$lwr_t  <- apply(boot_obj$t, 2, quantile,
                            probs = 0.025, na.rm = TRUE)
  pred_grid$upr_t  <- apply(boot_obj$t, 2, quantile,
                            probs = 0.975, na.rm = TRUE)
  
  # Back-transform with biological floor
  inv_xform <- function(x, type) {
    switch(type,
           none  = x,
           log10 = 10^x,
           log1p = expm1(x),
           logit = plogis(x) * 100)
  }
  
  lower_bound_for_response <- function(resp) {
    if (resp %in% c("n50", "busco", "maxrss_gb",
                    "ram_per_mb", "time_per_mb", "ram_mb_p",
                    "time_mb_p", "elapsed_sec", "genome_fraction_percent",
                    "mism_100kbp", "indel_100kbp")) 0 else NA_real_
  }
  
  floor_val <- lower_bound_for_response(response)
  
  pred_grid <- pred_grid %>%
    mutate(
      pred = inv_xform(pred_t, transform),
      lwr  = inv_xform(lwr_t,  transform),
      upr  = inv_xform(upr_t,  transform),
      response = response,
      transform = transform
    )
  
  if (!is.na(floor_val)) {
    pred_grid <- pred_grid %>%
      mutate(
        pred = pmax(pred, floor_val),
        lwr  = pmax(lwr,  floor_val),
        upr  = pmax(upr,  floor_val)
      )
  }
  
  pred_grid
}

# Programs to include in supplementary
shared_programs_sim <- c("flye", "spades_short", "abyss_short")
shared_programs_emp <- c("flye-ovl1000", "spades_short", "abyss_short")

# Responses shared between datasets
shared_responses <- c("n50", "busco", "time_per_mb", "ram_per_mb")
shared_responses_emp <- c("n50", "busco", "time_mb_p", "ram_mb_p")

# Predict for simulated
sim_preds <- map_dfr(shared_responses, function(resp) {
  get_model_predictions(simulated_single_fits[[resp]],
                        programs_to_keep = shared_programs_sim,
                        is_empirical = FALSE) %>%
    mutate(dataset = "Simulated")
})

# Predict for empirical (with response name remapping)
emp_response_map <- c(
  "n50"            = "n50",
  "busco"          = "busco",
  "time_per_mb"    = "time_mb_p",     # use proxy-normalized
  "ram_per_mb"     = "ram_mb_p"
)

emp_preds <- map_dfr(names(emp_response_map), function(canonical) {
  emp_resp <- emp_response_map[[canonical]]
  preds <- get_model_predictions(empirical_single_fits[[emp_resp]],
                                 programs_to_keep = shared_programs_emp,
                                 is_empirical = TRUE) %>%
    mutate(dataset = "Empirical",
           response = canonical,  # rename for joining
           program = case_when(
             program == "flye-ovl1000" ~ "flye",
             TRUE ~ program
           ))
})

# Combine and label
preds_combined <- bind_rows(sim_preds, emp_preds) %>%
  mutate(
    dataset = factor(dataset, levels = c("Simulated", "Empirical")),
    program = factor(program,
                     levels = c("flye", "spades_short", "abyss_short"),
                     labels = c("Flye", "SPAdes (short)", "ABySS (short)")),
    response = factor(response,
                      levels = c("n50", "busco",
                                 "time_per_mb", "ram_per_mb"),
                      labels = c("N50 (bp)",
                                 "BUSCO completeness (%)",
                                 "Runtime per Mbp (sec)",
                                 "RAM per Mbp (GB)"))
  )

# Faceted grid plot
fig_supp <- ggplot(preds_combined,
                   aes(x = depth_eff, y = pred,
                       color = dataset, fill = dataset,
                       linetype = dataset)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_grid(response ~ program, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = c(10, 30, 50, 75, 100)) +
  scale_color_manual(values = c("Simulated" = "#1b9e77",
                                "Empirical" = "#d95f02")) +
  scale_fill_manual(values  = c("Simulated" = "#1b9e77",
                                "Empirical" = "#d95f02")) +
  scale_linetype_manual(values = c("Simulated" = "solid",
                                   "Empirical" = "dashed")) +
  scale_y_log10() +
  labs(x = "Long-read depth (\u00d7)", y = NULL,
       color = "Dataset", fill = "Dataset", linetype = "Dataset") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.placement = "outside",
        strip.background = element_rect(fill = "grey95"),
        panel.spacing = unit(0.6, "lines"))

save_plot(fig_supp, "Figure_SX_sim_vs_empirical_supp", 10, 11)
