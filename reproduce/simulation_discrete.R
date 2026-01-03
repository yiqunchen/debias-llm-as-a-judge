# Discrete Prediction Simulation Runner
#
# DGP: Z ~ Unif({1,...,K}), Y | Z ~ N(mu_Z, sigma^2), Yhat = Z
#
# This simulation demonstrates that discrete predictions (like LLM ratings 1-5)
# require recalibration even when the prediction perfectly identifies the
# underlying category, because Yhat is on the wrong scale.

library(tidyverse)
library(future)
library(furrr)
library(here)
library(mgcv)
library(splines)

source(here("R/sim_discrete.R"))

# ============================================================================
# SIMULATION CONFIGURATION
# ============================================================================

SIM_CONFIG <- list(
  N_total = 2000L,                        # Fixed total sample size
  label_ratios = c(0.05, 0.10, 0.20),     # Labeled data fractions
  B = 500L,
  # Vary mu_3 to control bias magnitude
  # mu = (1, 2, mu_3) where mu_3 varies
  # E[Y] = (1 + 2 + mu_3)/3, E[Yhat] = 2
  # Naive bias = 2 - (3 + mu_3)/3 = (6 - 3 - mu_3)/3 = (3 - mu_3)/3
  mu_1 = 1,
  mu_2 = 2,
  mu_3_values = c(3, 4, 5, 6, 7, 8, 9),   # Varying third component
  sigma = 1.0,
  alpha = 0.10
)

# ============================================================================
# SINGLE SIMULATION REPLICATE
# ============================================================================

#' Run one simulation replicate
#'
#' @param N_total Total sample size
#' @param mu_3 Third component mean (varies)
#' @param label_ratio Fraction of data that's labeled (calibration)
#' @param mu_1 First component mean (fixed)
#' @param mu_2 Second component mean (fixed)
#' @param sigma Noise level
#' @param alpha Miscoverage level
#' @param seed Random seed
sim_one_rep_discrete <- function(N_total, mu_3, label_ratio,
                                  mu_1, mu_2, sigma, alpha = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  mu <- c(mu_1, mu_2, mu_3)

  # Generate data
  dgp <- generate_discrete_dgp(N_total, mu = mu, sigma = sigma)
  Y <- dgp$Y
  Yhat <- dgp$Yhat

  # Split into calibration (labeled) and test (unlabeled)
  m <- max(10L, round(label_ratio * N_total))
  n_test <- N_total - m

  cal_idx <- sample(N_total, m)
  test_idx <- setdiff(seq_len(N_total), cal_idx)

  Y_cal <- Y[cal_idx]
  Y_test <- Y[test_idx]
  Yhat_cal <- Yhat[cal_idx]
  Yhat_test <- Yhat[test_idx]

  # True mean on test set
  theta_true <- mean(Y_test)

  # Collect results
  results <- list()

  # Naive
  naive <- naive_discrete_point_and_ci(Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "Naive",
    theta_hat = naive$theta,
    var_hat = naive$var,
    ci_lower = naive$ci_lower,
    ci_upper = naive$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= naive$ci_lower & theta_true <= naive$ci_upper
  )

  # PPI
  ppi <- ppi_discrete_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "PPI",
    theta_hat = ppi$theta,
    var_hat = ppi$var,
    ci_lower = ppi$ci_lower,
    ci_upper = ppi$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= ppi$ci_lower & theta_true <= ppi$ci_upper
  )

  # PPI++
  ppi_pp <- ppi_pp_discrete_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "PPI++",
    theta_hat = ppi_pp$theta,
    var_hat = ppi_pp$var,
    ci_lower = ppi_pp$ci_lower,
    ci_upper = ppi_pp$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= ppi_pp$ci_lower & theta_true <= ppi_pp$ci_upper,
    lambda = ppi_pp$lambda
  )

  # EIF (per-category calibration - optimal)
  eif <- eif_discrete_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "EIF",
    theta_hat = eif$theta,
    var_hat = eif$var,
    ci_lower = eif$ci_lower,
    ci_upper = eif$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= eif$ci_lower & theta_true <= eif$ci_upper
  )

  # EIF-linear (linear calibration - suboptimal for discrete)
  eif_lin <- eif_linear_discrete_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "EIF-linear",
    theta_hat = eif_lin$theta,
    var_hat = eif_lin$var,
    ci_lower = eif_lin$ci_lower,
    ci_upper = eif_lin$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= eif_lin$ci_lower & theta_true <= eif_lin$ci_upper
  )

  # EIF-gam (GAM calibration)
  eif_gam <- eif_gam_discrete_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "EIF-gam",
    theta_hat = eif_gam$theta,
    var_hat = eif_gam$var,
    ci_lower = eif_gam$ci_lower,
    ci_upper = eif_gam$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= eif_gam$ci_lower & theta_true <= eif_gam$ci_upper
  )

  # EIF-spline (spline calibration)
  eif_spl <- eif_spline_discrete_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "EIF-spline",
    theta_hat = eif_spl$theta,
    var_hat = eif_spl$var,
    ci_lower = eif_spl$ci_lower,
    ci_upper = eif_spl$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= eif_spl$ci_lower & theta_true <= eif_spl$ci_upper
  )

  # Combine
  bind_rows(results) %>%
    mutate(
      N_total = N_total,
      mu_3 = mu_3,
      label_ratio = label_ratio,
      m = m,
      n_test = n_test
    )
}

# ============================================================================
# RUN FULL SIMULATION
# ============================================================================

run_discrete_simulation <- function(config = SIM_CONFIG, parallel = TRUE) {
  # Expand parameters: vary mu_3 and label_ratio
  params <- expand_grid(
    mu_3 = config$mu_3_values,
    label_ratio = config$label_ratios,
    rep = 1:config$B
  )

  cat("\n")
  cat("=======================================================\n")
  cat("       DISCRETE PREDICTION SIMULATION\n")
  cat("=======================================================\n")
  cat(sprintf("Total simulations: %d\n", nrow(params)))
  cat(sprintf("  - Monte Carlo reps: %d\n", config$B))
  cat(sprintf("  - N_total: %d (fixed)\n", config$N_total))
  cat(sprintf("  - mu = (%.0f, %.0f, mu_3) with mu_3 in {%s}\n",
              config$mu_1, config$mu_2, paste(config$mu_3_values, collapse = ", ")))
  cat(sprintf("  - label_ratios: %s\n", paste(config$label_ratios, collapse = ", ")))
  cat(sprintf("  - sigma=%.1f, alpha=%.2f\n", config$sigma, config$alpha))
  cat("\n")
  cat("DGP: Z ~ Unif({1,2,3}), Y|Z ~ N(mu_Z, sigma^2), Yhat = Z\n")
  cat("E[Yhat] = 2, E[Y] = (mu_1 + mu_2 + mu_3)/3\n")
  cat("Naive bias = 2 - E[Y] (increases with mu_3)\n")
  cat("=======================================================\n\n")

  N_total <- config$N_total
  mu_1 <- config$mu_1
  mu_2 <- config$mu_2

  if (parallel) {
    n_workers <- availableCores() - 1
    cat(sprintf("Running in PARALLEL mode with %d workers...\n\n", n_workers))
    plan(multisession, workers = n_workers)

    sim_file <- here("R/sim_discrete.R")

    results <- future_pmap_dfr(
      params,
      function(mu_3, label_ratio, rep) {
        source(sim_file, local = TRUE)

        sim_one_rep_discrete(
          N_total = N_total,
          mu_3 = mu_3,
          label_ratio = label_ratio,
          mu_1 = mu_1,
          mu_2 = mu_2,
          sigma = config$sigma,
          alpha = config$alpha,
          seed = rep * 1000 + mu_3 * 10 + round(label_ratio * 100)
        ) %>%
          mutate(rep = rep)
      },
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )

    plan(sequential)
  } else {
    cat("Running in SEQUENTIAL mode...\n\n")
    results <- pmap_dfr(
      params,
      function(mu_3, label_ratio, rep) {
        if (rep %% 100 == 0) {
          cat(sprintf("  Progress: rep %d / %d\n", rep, config$B))
        }

        sim_one_rep_discrete(
          N_total = N_total,
          mu_3 = mu_3,
          label_ratio = label_ratio,
          mu_1 = mu_1,
          mu_2 = mu_2,
          sigma = config$sigma,
          alpha = config$alpha,
          seed = rep * 1000 + mu_3 * 10 + round(label_ratio * 100)
        ) %>%
          mutate(rep = rep)
      }
    )
  }

  results
}

# ============================================================================
# SUMMARIZE RESULTS
# ============================================================================

summarize_discrete_results <- function(results) {
  results %>%
    filter(!is.na(theta_hat)) %>%
    group_by(method, mu_3, label_ratio) %>%
    summarize(
      n_reps = n(),
      coverage = mean(covers, na.rm = TRUE),
      coverage_se = sqrt(coverage * (1 - coverage) / n_reps),
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      bias = mean(theta_hat - theta_true, na.rm = TRUE),
      rmse = sqrt(mean((theta_hat - theta_true)^2, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(
      # Pretty label for label_ratio
      labeled_pct = paste0("Labeled: ", label_ratio * 100, "%")
    )
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if (sys.nframe() == 0) {
  cat("\n")
  cat("##########################################################\n")
  cat("#  DISCRETE PREDICTION SIMULATION                       #\n")
  cat("##########################################################\n")
  cat("\n")
  cat("DGP: Z ~ Unif({1,2,3}), Y|Z ~ N(mu_Z, sigma^2), Yhat = Z\n")
  cat("\n")
  cat("Key insight: Yhat perfectly identifies the mixture component,\n")
  cat("but E[Yhat] != E[Y] because Yhat is on the wrong scale!\n")
  cat("\n")
  cat("Setup: mu = (1, 2, mu_3) where mu_3 varies from 3 to 9\n")
  cat("  - E[Yhat] = 2 (always)\n")
  cat("  - E[Y] = (1 + 2 + mu_3)/3 = (3 + mu_3)/3\n")
  cat("  - Naive bias = 2 - E[Y] = (3 - mu_3)/3 (negative, grows with mu_3)\n")
  cat("\n")

  start_time <- Sys.time()
  results <- run_discrete_simulation(parallel = TRUE)
  elapsed <- difftime(Sys.time(), start_time, units = "mins")

  cat(sprintf("\nSimulation completed in %.1f minutes\n", as.numeric(elapsed)))

  write_csv(results, "reproduce/discrete_sim_results.csv")
  cat("Raw results saved to: reproduce/discrete_sim_results.csv\n")

  summary_df <- summarize_discrete_results(results)
  write_csv(summary_df, "reproduce/discrete_sim_summary.csv")
  cat("Summary saved to: reproduce/discrete_sim_summary.csv\n")

  cat("\n")
  cat("=======================================================\n")
  cat("  RESULTS PREVIEW (Labeled: 10%)\n")
  cat("=======================================================\n")
  summary_df %>%
    filter(label_ratio == 0.10) %>%
    select(method, mu_3, bias, coverage, rmse) %>%
    arrange(mu_3, method) %>%
    print(n = 40)

  cat("\n")
  cat("##########################################################\n")
  cat("#                 SIMULATION COMPLETE                    #\n")
  cat("##########################################################\n")
}
