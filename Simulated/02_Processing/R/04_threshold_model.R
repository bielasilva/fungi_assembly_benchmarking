# ============================================================================
# 04_threshold_model.R
# ----------------------------------------------------------------------------
# Fits models to predict response values across the depth gradient and estimate
# a "threshold" depth for each response and program, where further increases in
# depth yield diminishing returns.
#
# Provides:
#   - Model-based estimates of the depth threshold for each response and program,
#     where the threshold is defined as the depth at which a certain percentage
#     (e.g., 90%) of the total predicted improvement has been achieved.
#   - Saves the threshold estimates in CSV files for use in the paper and further analysis.
#
# Sources:
#   - 02_fit_models.R  (uses fitted model objects in single_fits / multi_fits)
# ============================================================================

library(tidyverse)
library(lme4)

# ---- Inverse transformations ----
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

# ---- Define whether higher or lower values are better ----
metric_direction <- tibble::tribble(
  ~response,                  ~direction,    ~lower_bound,
  "n50",                       "higher",     NA,
  "busco_complete",            "higher",     NA,
  "busco_complete_diff",       "lower_abs",  NA,
  "genome_fraction_percent",   "higher",     NA,
  "length_err_pct",            "lower_abs",  NA,
  "mism_100kbp",               "lower",      0, # non-negative count
  "indel_100kbp",              "lower",      0, # non-negative count
  "ram_per_mb",                "lower",      0,
  "maxrss_gb",                 "lower",      0,
  "time_per_mb",               "lower",      0,
  "elapsed_sec",               "lower",      0
)

# Define computational parameters for bootstrapping predictions
set.seed(123) # For reproducibility of bootstrapping results
n_cpus = 20 # Number of CPU cores to use for parallel processing
n_sim = 1000 # Number of bootstrap simulations

# ---- Find threshold depth for one program ----
## Threshold is the depth where 95% of total predicted improvement has occurred.

## Single
estimate_single <- function(fit,
                            depth_min = 10,
                            depth_max = 100,
                            depth_step = 5,
                            improvement_prop = 0.9) {

  mod       <- fit$model
  response  <- fit$response
  transform <- fit$transform
  
  programs <- levels(model.frame(mod)$program)
  depth_grid <- seq(depth_min, depth_max, by = depth_step)
  
  # Build prediction grid.
  pred_grid <- tidyr::expand_grid(
    program = programs,
    depth_eff = depth_grid
  ) %>%
    mutate(
      size_mbp_s      = 0, # Covariates are centered, so mean = 0.
      gc_content_s    = 0,
      gene_density_s  = 0,
      cds_per_gene_s  = 0
    )

  boot_obj <- lme4::bootMer(mod,
                    nsim = n_sim,
                    parallel = "multicore",
                    ncpus = n_cpus,
                    FUN = function(x) {predict(x,
                                        newdata = pred_grid,
                                        re.form = NA,
                                        allow.new.levels = TRUE)
    }
  )

  pred_grid$pred_t <- apply(boot_obj$t, 2, median,
                            na.rm = TRUE)

  pred_grid$lwr_t  <- apply(boot_obj$t, 2, quantile,
                            probs = 0.025, na.rm = TRUE)

  pred_grid$upr_t  <- apply(boot_obj$t, 2, quantile,
                            probs = 0.975, na.rm = TRUE)

  
  pred_grid <- pred_grid %>%
    mutate(
      pred = inverse_transform(pred_t, transform),
      ci.l = inverse_transform(ci.l_t, transform),
      ci.u = inverse_transform(ci.u_t, transform),
      response = response,
      transform = transform
    )
  
  direction <- metric_direction %>%
    filter(response == !!response) %>%
    pull(direction)
  
  if (length(direction) == 0) {
    stop("No direction defined for response: ", response)
  }

  # Truncate at biological floor for non-negative responses
  this_lower <- metric_direction$lower_bound[metric_direction$response == response]
  if (!is.na(this_lower)) {
    pred_grid <- pred_grid %>%
      mutate(
        pred = pmax(pred, this_lower),
        ci.l = pmax(ci.l, this_lower),
        ci.u = pmax(ci.u, this_lower)
      )
  }
  
  threshold_tbl <- pred_grid %>%
    group_by(response, program) %>%
    group_modify(~ {
      
      df <- .x
      
      if (direction == "higher") {
        start_value <- min(df$pred, na.rm = TRUE)
        best_value <- max(df$pred, na.rm = TRUE)
        # threshold_pred <- start_value + improvement_prop * (best_value - start_value)
        threshold_pred <- best_value * improvement_prop
        
        out <- df %>%
          filter(pred >= threshold_pred) %>%
          slice_head(n = 1)
          
      } else if (direction == "lower") {
        start_value <- max(df$pred, na.rm = TRUE)
        best_value <- min(df$pred, na.rm = TRUE)
        # threshold_pred <- start_value - improvement_prop * (start_value - best_value)
        threshold_pred <- best_value * (2 - improvement_prop)
        
        out <- df %>%
          filter(pred <= threshold_pred) %>%
          slice_head(n = 1)
        
      } else if (direction == "lower_abs") {
        start_value <- max(abs(df$pred), na.rm = TRUE)
        best_value <- min(abs(df$pred), na.rm = TRUE)
        # threshold_pred <- start_value - improvement_prop * (start_value - best_value)
        threshold_pred <- best_value * (2 - improvement_prop)
        
        out <- df %>%
          mutate(abs_response = abs(pred)) %>%
          filter(abs_response <= threshold_pred) %>%
          slice_head(n = 1)
      }
      
      if (nrow(out) == 0) {
        return(tibble(
          threshold_depth = NA_real_,
          threshold_pred = threshold_pred,
          best_value = best_value
        ))
      } else {
        tibble(
          threshold_depth = out$depth_eff[1],
          threshold_pred = threshold_pred,
          best_value = best_value
        )
      }
    }) %>%
    ungroup() %>%
    mutate(
      direction = direction
    )
  
  list(
    predictions = pred_grid,
    thresholds = threshold_tbl
  )
}

single_model_results <- imap(
  single_fits,
  \(x, response_id) {
    estimate_single(x$fit)
  }
) %>% compact()

readr::write_csv(
  map_dfr(single_model_results, ~ .x$thresholds),
  file.path(results_dir, "model_based_single_depth_thresholds.csv")
)

readr::write_csv(
  map_dfr(single_model_results, ~ .x$predictions),
  file.path(results_dir, "model_based_single_depth_predictions.csv")
)

save.image(file.path(results_dir, "04_threshold_model_single.RData"))

# Multi
estimate_multi <- function(fit,
                            depth_np_min = 10,
                            depth_np_max = 100,
                            depth_il_min = 10,
                            depth_il_max = 100,
                            depth_step = 5,
                            improvement_prop = 0.9) {

  mod       <- fit$model
  response  <- fit$response
  transform <- fit$transform
  
  programs <- levels(model.frame(mod)$program)
  
  depth_np_grid <- seq(depth_np_min, depth_np_max, by = depth_step)
  depth_il_grid <- seq(depth_il_min, depth_il_max, by = depth_step)
  
  pred_grid <- tidyr::expand_grid(
    program = programs,
    depth_np = depth_np_grid,
    depth_il = depth_il_grid
  ) %>%
    mutate(
        size_mbp_s      = 0,
        gc_content_s    = 0,
        gene_density_s  = 0,
        cds_per_gene_s  = 0,
        total_depth = depth_np + depth_il
    )
  
  boot_obj <- lme4::bootMer(mod,
                           nsim = n_sim,
                           parallel = "multicore",
                           ncpus = n_cpus,
                           FUN = function(x) {predict(x,
                                                      newdata = pred_grid,
                                                      re.form = NA,
                                                      allow.new.levels = TRUE)
                           })
  
  pred_grid$pred_t <- apply(boot_obj$t, 2, \(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE)))
  pred_grid$ci.l_t <- apply(boot_obj$t, 2, \(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
  pred_grid$ci.u_t <- apply(boot_obj$t, 2, \(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
  
  
  pred_grid <- pred_grid %>%
    mutate(
      pred = inverse_transform(pred_t, transform),
      ci.l = inverse_transform(ci.l_t, transform),
      ci.u = inverse_transform(ci.u_t, transform),
      response = response,
      transform = transform
    )
  
  direction <- metric_direction %>%
    filter(response == !!response) %>%
    pull(direction)
  
  if (length(direction) == 0) {
    stop("No direction defined for response: ", response)
  }
  
  # Truncate at biological floor for non-negative responses
  this_lower <- metric_direction$lower_bound[metric_direction$response == response]
  if (!is.na(this_lower)) {
    pred_grid <- pred_grid %>%
      mutate(
        pred = pmax(pred, this_lower),
        ci.l = pmax(ci.l, this_lower),
        ci.u = pmax(ci.u, this_lower)
      )
  }
  
  threshold_tbl <- pred_grid %>%
    group_by(response, transform, program) %>%
    group_modify(~ {
      
      df <- .x
      
      if (direction == "higher") {
        start_value <- min(df$pred, na.rm = TRUE)
        best_value <- max(df$pred, na.rm = TRUE)
        # threshold_pred <- start_value + improvement_prop * (best_value - start_value)
        threshold_pred <- best_value * improvement_prop

        
        out <- df %>%
          filter(pred >= threshold_pred) %>%
          arrange(depth_np, depth_il) %>%
          slice_head(n = 1)
          
      } else if (direction == "lower") {
        start_value <- max(df$pred, na.rm = TRUE)
        best_value <- min(df$pred, na.rm = TRUE)
        # threshold_pred <- start_value - improvement_prop * (start_value - best_value)
        threshold_pred <- best_value * (2 - improvement_prop)
        
        out <- df %>%
          filter(pred <= threshold_pred) %>%
          arrange(depth_np, depth_il) %>%
          slice_head(n = 1)
        
      } else if (direction == "lower_abs") {
        start_value <- max(abs(df$pred), na.rm = TRUE)
        best_value <- min(abs(df$pred), na.rm = TRUE)
        # threshold_pred <- start_value - improvement_prop * (start_value - best_value)
        threshold_pred <- best_value * (2 - improvement_prop)
        
        out <- df %>%
          mutate(abs_response = abs(pred)) %>%
          filter(abs_response <= threshold_pred) %>%
          arrange(depth_np, depth_il) %>%
          slice_head(n = 1)
      }
      
      if (nrow(out) == 0) {
        tibble(
          threshold_depth_np = NA_real_,
          threshold_depth_il = NA_real_,
          total_depth = NA_real_,
          predicted_at_threshold = NA_real_,
          threshold_pred = threshold_pred,
          best_value = best_value
        )
      } else {
        tibble(
          threshold_depth_np = out$depth_np[1],
          threshold_depth_il = out$depth_il[1],
          total_depth = out$total_depth[1],
          predicted_at_threshold = out$pred[1],
          threshold_pred = threshold_pred,
          best_value = best_value
        )
      }
    }) %>%
    ungroup() %>%
    mutate(
      direction = direction
    )
  
  list(
    predictions = pred_grid,
    thresholds = threshold_tbl
  )
}

multi_model_results <- imap(
  multi_fits,
  \(x, response_id) {
    estimate_multi(x$fit)
  }
) %>% compact()

readr::write_csv(
  map_dfr(multi_model_results, ~ .x$thresholds),
  file.path(results_dir, "model_based_multi_depth_thresholds.csv")
)

readr::write_csv(
  map_dfr(multi_model_results, ~ .x$predictions),
  file.path(results_dir, "model_based_multi_depth_predictions.csv")
)

save.image(file.path(results_dir, "04_threshold_model_multi.RData"))

## Polished
# For polishing we're considering fixed long-reads depth to find the minimum short-reads.
estimate_polished <- function(fit,
                              fixed_lr_depths = c(10, 20, 30, 40),
                              depth_il_min = 10,
                              depth_il_max = 100,
                              depth_step = 1,
                              improvement_prop = 0.9) {
  mod       <- fit$model
  response  <- fit$response
  transform <- fit$transform
  
  programs <- levels(model.frame(mod)$program)
  
  pred_grid <- tidyr::expand_grid(
    program = programs,
    depth_np = fixed_lr_depths,
    depth_il = seq(depth_il_min, depth_il_max, by = depth_step)
  ) %>%
    mutate(
      size_mbp_s      = 0,
      gc_content_s    = 0,
      gene_density_s  = 0,
      cds_per_gene_s  = 0,
      total_depth = depth_np + depth_il
    )
  
  boot_obj <- lme4::bootMer(mod,
                           nsim = n_sim,
                           parallel = "multicore",
                           ncpus = n_cpus,
                           FUN = function(x) {predict(x,
                                                      newdata = pred_grid,
                                                      re.form = NA,
                                                      allow.new.levels = TRUE)
                           })
  
  pred_grid$pred_t <- apply(boot_obj$t, 2, \(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE)))
  pred_grid$ci.l_t <- apply(boot_obj$t, 2, \(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE)))
  pred_grid$ci.u_t <- apply(boot_obj$t, 2, \(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
  
  
  pred_grid <- pred_grid %>%
    mutate(
      pred = inverse_transform(pred_t, transform),
      ci.l = inverse_transform(ci.l_t, transform),
      ci.u = inverse_transform(ci.u_t, transform),
      response = response,
      transform = transform
    )
  
  direction <- metric_direction %>%
    filter(response == !!response) %>%
    pull(direction)

  # Truncate at biological floor for non-negative responses
  this_lower <- metric_direction$lower_bound[metric_direction$response == response]
  if (!is.na(this_lower)) {
    pred_grid <- pred_grid %>%
      mutate(
        pred = pmax(pred, this_lower),
        ci.l = pmax(ci.l, this_lower),
        ci.u = pmax(ci.u, this_lower)
      )
  }
  
  threshold_tbl <- pred_grid %>%
    group_by(response, transform, program, depth_np) %>%
    group_modify(~ {
      
      df <- .x

      if (direction == "higher") {
        start_value <- min(df$pred, na.rm = TRUE)
        best_value <- max(df$pred, na.rm = TRUE)
        threshold_pred <- start_value + improvement_prop * (best_value - start_value)
        
        out <- df %>%
          filter(pred >= threshold_pred) %>%
          arrange(depth_il) %>%
          slice_head(n = 1)
          
      } else if (direction == "lower") {
        start_value <- max(df$pred, na.rm = TRUE)
        best_value <- min(df$pred, na.rm = TRUE)
        threshold_pred <- start_value - improvement_prop * (start_value - best_value)
        
        out <- df %>%
          filter(pred <= threshold_pred) %>%
          arrange(depth_il) %>%
          slice_head(n = 1)
        
      } else if (direction == "lower_abs") {
        start_value <- max(abs(df$pred), na.rm = TRUE)
        best_value <- min(abs(df$pred), na.rm = TRUE)
        threshold_pred <- start_value - improvement_prop * (start_value - best_value)
        
        out <- df %>%
          mutate(abs_response = abs(pred)) %>%
          filter(abs_response <= threshold_pred) %>%
          arrange(depth_il) %>%
          slice_head(n = 1)
      }
      
      if (nrow(out) == 0) {
        tibble(
          threshold_depth_il = NA_real_,
          predicted_at_threshold = NA_real_,
          threshold_pred = threshold_pred,
          best_value = best_value
        )
      } else {
        tibble(
          threshold_depth_il = out$depth_il[1],
          predicted_at_threshold = out$pred[1],
          threshold_pred = threshold_pred,
          best_value = best_value
        )
      }
    }) %>%
    ungroup() %>%
    rename(fixed_depth_np = depth_np) %>%
    mutate(
      direction = direction,
      threshold_rule = case_when(
        direction == "higher" ~ paste0(improvement_prop * 100, "% of predicted maximum"),
        direction == "lower" ~ paste0(improvement_prop * 100, "% of predicted minimum"),
        direction == "lower_abs" ~ paste0(improvement_prop * 100, "% of minimum absolute value")
      )
    )
  
  list(
    predictions = pred_grid,
    thresholds = threshold_tbl
  )
}

polished_model_results <- imap(
  within_strategy_fits$polished,
  \(x, response_id) {
    estimate_polished(
      fit = x$fit
    )
  }
)  %>% compact()

readr::write_csv(
  map_dfr(polished_model_results, ~ .x$thresholds),
  file.path(results_dir, "model_based_polished_depth_thresholds.csv")
)

readr::write_csv(
  map_dfr(polished_model_results, ~ .x$predictions),
  file.path(results_dir, "model_based_polished_depth_predictions.csv")
)

save.image(file.path(results_dir, "04_threshold_model_polished.RData"))

## Arrange

single_threshold_results_for_paper <- 
  map_dfr(single_model_results, ~ .x$thresholds) %>% select(response, program, threshold_depth) %>%
  mutate(response = str_replace_all(response,response_labels),
          program = str_replace_all(program,program_labels)) %>%
  pivot_wider(names_from = program, values_from = threshold_depth)

readr::write_csv(
  single_threshold_results_for_paper, 
  file.path(results_dir, "single_thresholds_for_paper.csv")
)

multi_threshold_results_for_paper <- 
  map_dfr(multi_model_results, ~ .x$thresholds) %>% 
  mutate(threshold_depth = paste0(threshold_depth_np, " | ", threshold_depth_il),
          response = str_replace_all(response,response_labels),
          program = str_replace_all(program,program_labels)) %>%
  select(response, program, threshold_depth) %>%
  pivot_wider(names_from = program, values_from = threshold_depth)

readr::write_csv(
  multi_threshold_results_for_paper, 
  file.path(results_dir, "multi_thresholds_for_paper.csv")
)

polished_threshold_results_for_paper <- 
  map_dfr(polished_model_results, ~ .x$thresholds) %>%
  mutate(threshold_depth = paste0(fixed_depth_np, " | ", threshold_depth_il)) %>%
  select(response, program, threshold_depth) %>%
  pivot_wider(names_from = program, values_from = threshold_depth)

readr::write_csv(
  polished_threshold_results_for_paper, 
  file.path(results_dir, "polished_thresholds_for_paper.csv")
)

# # Plot
# plot_single_threshold_curve <- function(threshold_object,
#                                         response_label = "Predicted response") {
  
  
#   pred <- threshold_object$predictions %>% 
#             mutate(program=factor(program, levels=program_order))
  
#   thr  <- threshold_object$thresholds
  
  
#   ggplot(pred, aes(x = depth_eff, y = pred, color = program)) +
#     geom_ribbon(aes(ymin = ci.l, ymax = ci.u, fill = program), 
#                 show.legend = FALSE, linewidth = 0, alpha = 0.2) +
#     geom_line(linewidth = 1) +
#     geom_vline(
#       data = thr,
#       aes(xintercept = threshold_depth, color = program),
#       linetype = "dashed",
#       alpha = 0.6,
#       show.legend = FALSE
#     ) +
#     scale_color_manual(values=colors_programs, label=program_labels, name="Program") +
#     scale_fill_manual(values=colors_programs, label=program_labels, name="Program") +
#     labs(
#       x = "Sequencing depth (X)",
#       y = response_label,
#       color = "Assembler"
#     ) +
#     theme_minimal(base_size = 13) +
#     theme(
#       panel.grid.minor = element_blank()
#     )
# }

# n50_multi <- estimate_multi(
#   fit = multi_fits$n50$fit
# )

# plot_single_threshold_curve(
#   single_model_results$n50,
#   response_label = "Predicted N50 (Mbp)"
# ) + scale_y_continuous(labels=scales::label_number(scale=1e-6), breaks = 1e6*c(0,1,2,3,4))

# plot_single_threshold_curve(
#   single_model_results$mism_100kbp,
#   response_label = "Mism"
# ) + ggbreak::scale_y_break(c(10, 500), scales = 0.1) 


#   geom_text(data = single_model_results$mism_100kbp$thr,
#             aes(x = 30, y = 7, label = threshold_depth),
#             position = "dodge")

# plot_single_threshold_curve(
#   single_model_results$indel_100kbp,
#   response_label = "Indel"
# ) + ggbreak::scale_y_break(c(10, 450), scales = 0.1)
