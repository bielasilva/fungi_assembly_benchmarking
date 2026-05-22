# ============================================================================
# 02_fit_models.R
# ----------------------------------------------------------------------------
# Fit linear mixed-effects models for empirical data with response transforms,
# DHARMa residual diagnostics, and exported coefficient tables.
#
# Sources:
#   - 01_data_prep.R (must be sourced first to populate `data_empirical`)
#
# Outputs (under Empirical/R/results/):
#   - <model>_global_effects.csv     Type III chi-square per fixed effect
#   - <model>_fixed_effects.csv      Coefficient estimates with CIs
#   - <model>_depth_slopes.csv       Per-program depth trends
#   - <model>_genome_slopes.csv      Per-program genome-size trends
#   - <model>_prog_pairs.csv         Tukey contrasts among programs
#   - <model>_depth_pairs.csv        Tukey contrasts among depths
#   - <model>_dharma_summary.csv     DHARMa formal test results
#   - <model>_dharma_plots/*.png     DHARMa diagnostic plots per response
#   - all_fixed_effects.csv          Stacked coefficient table for forest plots
#   - model_summary_table.csv        One-row-per-model overview
# ============================================================================
mirai::daemons(5)
mirai::everywhere(
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(emmeans)
  library(broom.mixed)
  library(car)
  library(DHARMa)
  library(performance)
  library(here)
  library(mirai)
})
)

# ----------------------------------------------------------------------------
# Significance flagger
# ----------------------------------------------------------------------------

sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    p < 0.1   ~ ".",
    TRUE      ~ ""
  )
}

# ----------------------------------------------------------------------------
# Response transforms
# ----------------------------------------------------------------------------
# Each response gets a transform that makes Gaussian residuals plausible.
# The named list maps response -> transform name; helper functions apply them.

response_transforms <- c(
  n50            = "log10",
  busco_complete = "logit",
  # mism_100kbp    = "log1p",
  # indel_100kbp   = "log1p",
  maxrss_gb      = "log10",
  ram_mb         = "log10",
  time_mb        = "log10",
  ram_mb_p       = "log10",
  time_mb_p      = "log10",
  elapsed_sec    = "log10"
)

apply_transform <- function(x, type, eps = 0.01) {
  switch(
    type,
    none  = x,
    log10 = log10(x + eps),
    log1p = log1p(x),
    logit = qlogis(pmin(pmax(x / 100, eps), 1 - eps)),
    stop("Unknown transform: ", type)
  )
}

inverse_transform <- function(x, type, eps = 0.01) {
  switch(
    type,
    none  = x,
    log10 = 10^x - eps,
    log1p = expm1(x),
    logit = plogis(x) * 100,
    stop("Unknown transform: ", type)
  )
}

# ----------------------------------------------------------------------------
# Model fitting helpers
# ----------------------------------------------------------------------------

mk_form <- function(response, rhs) {
  as.formula(sprintf("%s ~ %s", response, rhs), env = parent.frame())
}

# Fit a single-strategy model (short-only, long-only) with depth_eff as the
# depth variable and genome_size_est_s as the genome covariate.
fit_block_single <- function(df, response, transform_type) {
  
  dff <- df %>%
    drop_na(all_of(response), genome_size_est_s) %>%
    filter(depth_np != 99, depth_il != 99)  # drop OG-depth rows
  
  # Apply transform
  dff[[paste0(response, "_t")]] <- apply_transform(dff[[response]], transform_type)
  response_t <- paste0(response, "_t")
  
  rhs_single <- paste(
    "program * depth_eff + genome_size_est_s",
    "+ program:genome_size_est_s",
    "+ depth_eff:genome_size_est_s",
    "+ (1 | sample)",
    sep = " "
  )
  
  fml <- mk_form(response_t, rhs_single)
  mod <- lme4::lmer(fml, data = dff)
  
  list(
    model         = mod,
    data          = dff,
    response      = response,
    response_t    = response_t,
    transform     = transform_type,
    formula       = format(fml)
  )
}

# Fit a multi-strategy model (hybrid, polished) with separate NP and IL depth.
fit_block_multi <- function(df, response, transform_type) {
  
  dff <- df %>%
    drop_na(all_of(response), genome_size_est_s) %>%
    filter(depth_np != 99, depth_il != 99) %>% 
    select(where(~ !all(is.na(.))))
  
  dff[[paste0(response, "_t")]] <- apply_transform(dff[[response]], transform_type)
  response_t <- paste0(response, "_t")
  
  rhs_multi <- paste(
    "program * depth_np * depth_il + genome_size_est_s",
    "+ program:genome_size_est_s",
    "+ depth_np:genome_size_est_s",
    "+ depth_il:genome_size_est_s",
    "+ (1 | sample)",
    sep = " "
  )
  
  fml <- mk_form(response_t, rhs_multi)
  mod <- lme4::lmer(fml, data = dff)
  
  list(
    model         = mod,
    data          = dff,
    response      = response,
    response_t    = response_t,
    transform     = transform_type,
    formula       = format(fml)
  )
}

# ----------------------------------------------------------------------------
# Post-fit extraction
# ----------------------------------------------------------------------------

extract_global_effects <- function(fit) {
  car::Anova(fit$model, type = 3) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    rename(df = Df, chi_sq = Chisq, p_value = `Pr(>Chisq)`) %>%
    mutate(
      signif    = sig_stars(p_value),
      response  = fit$response,
      transform = fit$transform,
      .before   = 1
    ) %>%
    arrange(p_value)
}

extract_fixed_effects <- function(fit) {
  broom.mixed::tidy(fit$model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      response  = fit$response,
      transform = fit$transform,
      .before   = 1
    )
}

extract_depth_slopes_single <- function(fit) {
  emmeans::emtrends(fit$model, ~ program,
                    var = "depth_eff",
                    pbkrtest.limit = 70000) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    rename(
      slope   = depth_eff.trend,
      t_ratio = t.ratio,
      p_value = p.value
    ) %>%
    mutate(
      signif    = sig_stars(p_value),
      response  = fit$response,
      transform = fit$transform,
      slope_type = "eff",
      .before   = 1
    )
}

extract_depth_slopes_multi <- function(fit) {
  np <- emmeans::emtrends(fit$model, ~ program, 
                          var = "depth_np",
                          pbkrtest.limit = 70000) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    rename(slope = depth_np.trend, t_ratio = t.ratio, p_value = p.value) %>%
    mutate(slope_type = "np", .before = 1)
  
  il <- emmeans::emtrends(fit$model, ~ program,
                          var = "depth_il",
                          pbkrtest.limit = 70000) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    rename(slope = depth_il.trend, t_ratio = t.ratio, p_value = p.value) %>%
    mutate(slope_type = "il", .before = 1)
  
  bind_rows(np, il) %>%
    mutate(
      signif    = sig_stars(p_value),
      response  = fit$response,
      transform = fit$transform,
      .before   = 1
    )
}

extract_genome_slopes <- function(fit) {
  emmeans::emtrends(fit$model, ~ program,
                    var = "genome_size_est_s",
                    pbkrtest.limit = 70000) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    rename(slope = genome_size_est_s.trend, t_ratio = t.ratio, p_value = p.value) %>%
    mutate(
      signif    = sig_stars(p_value),
      covariate = "genome_size_est_s",
      response  = fit$response,
      transform = fit$transform,
      .before   = 1
    )
}


extract_prog_pairs_single <- function(fit) {
  depth_grid <- sort(unique(fit$data$depth_eff))
  emmeans::emmeans(fit$model, ~ program | depth_eff,
                   at = list(depth_eff = depth_grid), 
                   pbkrtest.limit = 70000) %>%
    emmeans::contrast(method = "tukey") %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    rename(t_ratio = t.ratio, p_value = p.value) %>%
    mutate(
      signif          = sig_stars(p_value),
      response        = fit$response,
      transform       = fit$transform,
      comparison_type = "program_at_depth",
      .before         = 1
    )
}

extract_depth_pairs_single <- function(fit) {
  depth_grid <- sort(unique(fit$data$depth_eff))
  emmeans::emmeans(fit$model, ~ depth_eff | program,
                   at = list(depth_eff = depth_grid),
                   pbkrtest.limit = 70000) %>%
    emmeans::contrast(method = "tukey") %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    rename(t_ratio = t.ratio, p_value = p.value) %>%
    mutate(
      signif          = sig_stars(p_value),
      response        = fit$response,
      transform       = fit$transform,
      comparison_type = "depth_within_program",
      .before         = 1
    )
}

# ----------------------------------------------------------------------------
# DHARMa diagnostics
# ----------------------------------------------------------------------------

run_dharma <- function(fit, n_sim = 1000, plot_dir = NULL) {
  sim <- DHARMa::simulateResiduals(fit$model, n = n_sim, plot = FALSE)
  
  # Formal tests
  ks_test    <- DHARMa::testUniformity(sim, plot = FALSE)
  disp_test  <- DHARMa::testDispersion(sim, plot = FALSE)
  out_test   <- DHARMa::testOutliers(sim, plot = FALSE)
  
  # Save diagnostic plot if a directory is provided
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
    png(file.path(plot_dir, paste0(fit$response, "_dharma.png")),
        width = 1200, height = 900, res = 120)
    plot(sim)
    dev.off()
  }
  
  tibble(
    response          = fit$response,
    transform         = fit$transform,
    ks_statistic      = unname(ks_test$statistic),
    ks_p              = ks_test$p.value,
    dispersion_ratio  = unname(disp_test$statistic),
    dispersion_p      = disp_test$p.value,
    outlier_p         = out_test$p.value,
    ks_pass           = ks_test$p.value > 0.05,
    dispersion_pass   = disp_test$p.value > 0.05,
    outlier_pass      = out_test$p.value > 0.05
  )
}

# ----------------------------------------------------------------------------
# Wrapper: fit -> extract -> save
# ----------------------------------------------------------------------------

fit_and_export_single <- function(df, response, transform_type, model_label) {
  message("Fitting single model: ", response, " (transform: ", transform_type, ")")
  
  fit <- fit_block_single(df, response, transform_type)
  
  # Extract tables
  global_effects <- extract_global_effects(fit)
  fixed_effects  <- extract_fixed_effects(fit)
  depth_slopes   <- extract_depth_slopes_single(fit)
  genome_slopes  <- extract_genome_slopes(fit)
  prog_pairs     <- extract_prog_pairs_single(fit)
  depth_pairs    <- extract_depth_pairs_single(fit)
  
  # Diagnostics
  plot_dir <- file.path(results_dir, paste0(model_label, "_dharma_plots"))
  dharma_summary <- run_dharma(fit, plot_dir = plot_dir)
  
  # Save tables
  write_csv(global_effects, file.path(results_dir, paste0(model_label, "_", response, "_global_effects.csv")))
  write_csv(fixed_effects,  file.path(results_dir, paste0(model_label, "_", response, "_fixed_effects.csv")))
  write_csv(depth_slopes,   file.path(results_dir, paste0(model_label, "_", response, "_depth_slopes.csv")))
  write_csv(genome_slopes,  file.path(results_dir, paste0(model_label, "_", response, "_genome_slopes.csv")))
  write_csv(prog_pairs,     file.path(results_dir, paste0(model_label, "_", response, "_prog_pairs.csv")))
  write_csv(depth_pairs,    file.path(results_dir, paste0(model_label, "_", response, "_depth_pairs.csv")))
  write_csv(dharma_summary, file.path(results_dir, paste0(model_label, "_", response, "_dharma_summary.csv")))
  
  list(
    model_label    = model_label,
    fit            = fit,
    global_effects = global_effects,
    fixed_effects  = fixed_effects,
    depth_slopes   = depth_slopes,
    genome_slopes  = genome_slopes,
    prog_pairs     = prog_pairs,
    depth_pairs    = depth_pairs,
    dharma_summary = dharma_summary
  )
}

fit_and_export_multi <- function(df, response, transform_type, model_label) {
  message("Fitting multi model: ", response, " (transform: ", transform_type, ")")
  
  fit <- fit_block_multi(df, response, transform_type)
  
  global_effects <- extract_global_effects(fit)
  fixed_effects  <- extract_fixed_effects(fit)
  depth_slopes   <- extract_depth_slopes_multi(fit)
  genome_slopes  <- extract_genome_slopes(fit)
  
  plot_dir <- file.path(results_dir, paste0(model_label, "_dharma_plots"))
  dharma_summary <- run_dharma(fit, plot_dir = plot_dir)
  
  write_csv(global_effects, file.path(results_dir, paste0(model_label, "_", response, "_global_effects.csv")))
  write_csv(fixed_effects,  file.path(results_dir, paste0(model_label, "_", response, "_fixed_effects.csv")))
  write_csv(depth_slopes,   file.path(results_dir, paste0(model_label, "_", response, "_depth_slopes.csv")))
  write_csv(genome_slopes,  file.path(results_dir, paste0(model_label, "_", response, "_genome_slopes.csv")))
  write_csv(dharma_summary, file.path(results_dir, paste0(model_label, "_", response, "_dharma_summary.csv")))
  
  list(
    model_label    = model_label,
    fit            = fit,
    global_effects = global_effects,
    fixed_effects  = fixed_effects,
    depth_slopes   = depth_slopes,
    genome_slopes  = genome_slopes,
    dharma_summary = dharma_summary
  )
}

# ----------------------------------------------------------------------------
# Run all models
# ----------------------------------------------------------------------------
# Splits the data into single-strategy (short-only and long-only) and
# multi-strategy (hybrid and polished) subsets, fits each response, and
# saves all outputs.

if (!exists("data_empirical")) {
  stop("data_empirical not found. Source 01_data_prep.R first.")
}

data_single <- data_empirical %>%
              filter(strategy %notin% c("hybrid", "polished"))
data_multi  <- data_empirical %>%
              filter(strategy %in% c("hybrid", "polished"),
                    program %notin% c("pilon", "polypolish"))

# --- Single-strategy fits ---
single_fits <- imap(response_transforms, in_parallel(function(transform_type, response) {
  if (response == "mism_100kbp" || response == "indel_100kbp") {
    # Skip because this metric is only in the polished subset
    return(NULL)
  }
  fit_and_export_single(
    df             = data_single,
    response       = response,
    transform_type = transform_type,
    model_label    = "empirical_single"
  )
  }, fit_and_export_single=fit_and_export_single,fit_block_single=fit_block_single,
    data_single=data_single, apply_transform=apply_transform, mk_form=mk_form,
    extract_global_effects=extract_global_effects, sig_stars=sig_stars,
    extract_fixed_effects=extract_fixed_effects, extract_genome_slopes=extract_genome_slopes,
    extract_depth_slopes_single=extract_depth_slopes_single, extract_prog_pairs_single=extract_prog_pairs_single,
    extract_depth_pairs_single=extract_depth_pairs_single, results_dir=results_dir,
    run_dharma=run_dharma))
single_fits <- compact(single_fits)

saveRDS(single_fits, file.path(results_dir, "empirical_single_fits_list.rds"))

# --- Single-strategy fits for each Seq Strategy ---
within_strategy_fits <- list()

for (strat in c("short_reads", "long_reads")) {
  data_strat <- data_single %>% filter(strategy == strat)
  
  fits_this_strat <- imap(response_transforms, in_parallel(function(transform_type, response) {
    
    fit_and_export_single(
      df             = data_strat,
      response       = response,
      transform_type = transform_type,
      model_label    = paste0("empirical_single_", strat)
    )
  }, fit_and_export_single=fit_and_export_single,fit_block_single=fit_block_single,
    data_strat=data_strat, apply_transform=apply_transform, mk_form=mk_form,
    extract_global_effects=extract_global_effects, sig_stars=sig_stars,
    extract_fixed_effects=extract_fixed_effects, extract_genome_slopes=extract_genome_slopes,
    extract_depth_slopes_single=extract_depth_slopes_single, extract_prog_pairs_single=extract_prog_pairs_single,
    extract_depth_pairs_single=extract_depth_pairs_single, results_dir=results_dir,
    run_dharma=run_dharma,strat=strat))
  
  within_strategy_fits[[strat]] <- compact(fits_this_strat)
}

message("Empirical model fitting for single-strategy done.")
save.image(paste(results_dir, "02_fit_models_single.RData", sep = "/"))
saveRDS(within_strategy_fits, file.path(results_dir, "empirical_within_strategy_fits_list.rds"))

# --- Multi-strategy fits ---
multi_fits <- imap(response_transforms, in_parallel(function(transform_type, response) {

  fit_and_export_multi(
    df             = data_multi,
    response       = response,
    transform_type = transform_type,
    model_label    = "empirical_multi"
  )
  }, fit_and_export_multi=fit_and_export_multi,fit_block_multi=fit_block_multi,
    data_multi=data_multi, apply_transform=apply_transform, mk_form=mk_form,
    extract_global_effects=extract_global_effects, sig_stars=sig_stars,
    extract_fixed_effects=extract_fixed_effects, extract_genome_slopes=extract_genome_slopes,
    extract_depth_slopes_multi=extract_depth_slopes_multi, results_dir=results_dir,
    run_dharma=run_dharma))
multi_fits <- compact(multi_fits)

message("Empirical model fitting for multi-strategy done.")
save.image(paste(results_dir, "02_fit_models_multi.RData", sep = "/"))
saveRDS(multi_fits, file.path(results_dir, "empirical_multi_fits_list.rds"))

# ----------------------------------------------------------------------------
# Stacked outputs for forest plots and overview
# ----------------------------------------------------------------------------

all_fixed_effects <- bind_rows(
  map_dfr(single_fits, ~ .x$fixed_effects %>% mutate(model_set = "single")),
  map_dfr(multi_fits,  ~ .x$fixed_effects %>% mutate(model_set = "multi")),
    map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
      f$fixed_effects %>%
        mutate(model_set = strat)
    })
  })
)
write_csv(all_fixed_effects, file.path(results_dir, "empirical_all_fixed_effects.csv"))

all_global_effects <- bind_rows(
  map_dfr(single_fits, ~ .x$global_effects %>% mutate(model_set = "single")),
  map_dfr(multi_fits,  ~ .x$global_effects %>% mutate(model_set = "multi")),
  map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
      f$global_effects %>%
        mutate(model_set = strat)
    })
  })
)
write_csv(all_global_effects, file.path(results_dir, "empirical_all_global_effects.csv"))

all_genome_slopes <- bind_rows(
  map_dfr(single_fits, ~ .x$genome_slopes %>% mutate(model_set = "single")),
  map_dfr(multi_fits,  ~ .x$genome_slopes %>% mutate(model_set = "multi")),
  map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
      f$genome_slopes %>%
        mutate(model_set = strat)
    })
  })
)
write_csv(all_genome_slopes, file.path(results_dir, "empirical_all_genome_slopes.csv"))

all_slopes_effects <- bind_rows(
  map_dfr(single_fits, ~ .x$depth_slopes %>% mutate(model_set = "single", response = .x$fit$response)),
  map_dfr(multi_fits,  ~ .x$depth_slopes %>% mutate(model_set = "multi", response = .x$fit$response)),
  map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
      f$depth_slopes %>%
        mutate(model_set = strat, response = f$fit$response)
    })
  })
)
write_csv(all_slopes_effects, file.path(results_dir, "empirical_all_slopes_effects.csv"))

all_global_effects <- bind_rows(
  map_dfr(single_fits, ~ .x$global_effects %>% mutate(model_set = "single")),
  map_dfr(multi_fits,  ~ .x$global_effects %>% mutate(model_set = "multi")),
  map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
      f$global_effects %>%
        mutate(model_set = strat, response = f$fit$response)
    })
  })
)
write_csv(all_global_effects, file.path(results_dir, "empirical_all_global_effects.csv"))

all_dharma <- bind_rows(
  map_dfr(single_fits, ~ .x$dharma_summary %>% mutate(model_set = "single")),
  map_dfr(multi_fits,  ~ .x$dharma_summary %>% mutate(model_set = "multi")),
  map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
      f$dharma_summary %>%
        mutate(model_set = strat)
    })
  })
)
write_csv(all_dharma, file.path(results_dir, "empirical_all_dharma_diagnostics.csv"))

model_summary_table <- bind_rows(
  map_dfr(single_fits, function(f) {
    tibble(
      model_set = f$model_label,
      response  = f$fit$response,
      transform = f$fit$transform,
      n_obs     = nobs(f$fit$model),
      aic       = AIC(f$fit$model),
      ks_pass   = f$dharma_summary$ks_pass,
      disp_pass = f$dharma_summary$dispersion_pass
    )
  }),
  map_dfr(multi_fits, function(f) {
    tibble(
      model_set = f$model_label,
      response  = f$fit$response,
      transform = f$fit$transform,
      n_obs     = nobs(f$fit$model),
      aic       = AIC(f$fit$model),
      ks_pass   = f$dharma_summary$ks_pass,
      disp_pass = f$dharma_summary$dispersion_pass
    )
  }),
  map_dfr(names(within_strategy_fits), function(strat) {
    map_dfr(within_strategy_fits[[strat]], function(f) {
    tibble(
      model_set = f$model_label,
      response  = f$fit$response,
      transform = f$fit$transform,
      n_obs     = nobs(f$fit$model),
      aic       = AIC(f$fit$model),
      ks_pass   = f$dharma_summary$ks_pass,
      disp_pass = f$dharma_summary$dispersion_pass
    )})
  })
)
write_csv(model_summary_table, file.path(results_dir, "empirical_model_summary_table.csv"))

message("Empirical model fitting complete. ", nrow(model_summary_table), " models fit.")
save.image(file.path(results_dir, "02_fit_models.RData"))
