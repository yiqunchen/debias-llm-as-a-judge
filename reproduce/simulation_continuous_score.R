# Continuous Score Simulation Runner
# Comparing estimation methods for continuous scores
#
# Structure:
# 1. Generate N observations from DGP: Y = f(X1,...,X5) + noise
# 2. Split: 20% for model training, 80% for evaluation
# 3. Train models (good/medium/poor) on training set
# 4. Predict Yhat on evaluation set
# 5. Split evaluation into calibration (labeled) and test (unlabeled)
# 6. Apply estimators and compare

library(tidyverse)
library(future)
library(furrr)
library(here)
library(mgcv)
library(splines)
library(randomForest)
library(rpart)

# Source functions directly
source(here("R/sim_continuous.R"))

# ============================================================================
# SIMULATION CONFIGURATION
# ============================================================================

SIM_CONFIG <- list(
  N_total = c(500L, 1000L, 2000L, 5000L), # Total sample sizes to try
  train_frac = 0.20,                      # Fraction for model training
  label_ratios = c(0.05, 0.10, 0.20),     # Fraction of eval data that's labeled
  B = 500L,                               # Monte Carlo replicates
  mu_x = 1,                               # X ~ N(mu_x, 1) - same for train & test
  model_qualities = c("good", "medium", "poor-OLS"),
  sigma_y = 2.0,                          # Y noise
  alpha = 0.10
)

# ============================================================================
# SINGLE SIMULATION REPLICATE
# ============================================================================

#' Run one simulation replicate
#'
#' KEY DESIGN: Train and test on SAME distribution X ~ N(mu_x, 1)
#' poor-tree will show bias even without distribution shift (not OLS!)
#'
#' @param N_total Total sample size
#' @param model_quality One of "good", "medium", "poor-OLS", "poor-tree"
#' @param label_ratio Fraction of eval data that's labeled
#' @param train_frac Fraction of data for model training
#' @param mu_x Mean of X for ALL data
#' @param sigma_y Y noise level
#' @param alpha Miscoverage level
#' @param seed Random seed
#' @return Tibble with results for all methods
sim_one_rep_continuous <- function(N_total, model_quality, label_ratio,
                                    train_frac, mu_x,
                                    sigma_y, alpha = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate ALL data from SAME distribution
  dgp <- generate_continuous_dgp(N_total, mu_x = mu_x, sigma_y = sigma_y)
  X <- dgp$X
  Y <- dgp$Y

  # Split: train_frac for training, rest for evaluation
  N_train <- floor(train_frac * N_total)
  train_idx <- sample(N_total, N_train)
  eval_idx <- setdiff(seq_len(N_total), train_idx)

  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  X_eval <- X[eval_idx, , drop = FALSE]
  Y_eval <- Y[eval_idx]

  # Train model on training data
  train_fn <- switch(model_quality,
    "good" = train_model_good,
    "medium" = train_model_medium,
    "poor-OLS" = train_model_poor_ols,
    "poor-tree" = train_model_poor_tree
  )
  model <- train_fn(X_train, Y_train)

  # Predict on evaluation set (SAME distribution as training)
  Yhat_eval <- model(X_eval)

  # Split evaluation into calibration (labeled) and test (unlabeled)
  N_eval <- length(Y_eval)
  m <- max(10L, round(label_ratio * N_eval))
  n_test <- N_eval - m

  cal_idx <- sample(N_eval, m)
  test_idx <- setdiff(seq_len(N_eval), cal_idx)

  Y_cal <- Y_eval[cal_idx]
  Y_test <- Y_eval[test_idx]
  Yhat_cal <- Yhat_eval[cal_idx]
  Yhat_test <- Yhat_eval[test_idx]

  # True mean on test set (target of inference)
  theta_true <- mean(Y_test)

  # Collect results
  results <- list()

  # Naive estimator
  naive <- naive_continuous_point_and_ci(Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "Naive",
    theta_hat = naive$theta,
    var_hat = naive$var,
    ci_lower = naive$ci_lower,
    ci_upper = naive$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= naive$ci_lower & theta_true <= naive$ci_upper
  )

  # PPI estimator
  ppi <- ppi_continuous_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
  results[[length(results) + 1]] <- tibble(
    method = "PPI",
    theta_hat = ppi$theta,
    var_hat = ppi$var,
    ci_lower = ppi$ci_lower,
    ci_upper = ppi$ci_upper,
    theta_true = theta_true,
    covers = theta_true >= ppi$ci_lower & theta_true <= ppi$ci_upper
  )

  # PPI++ estimator
  ppi_pp <- ppi_pp_continuous_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha)
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

  # EIF with different calibration methods
  for (cal_method in c("linear", "gam", "spline")) {
    eif <- tryCatch(
      eif_continuous_point_and_ci(Y_cal, Yhat_cal, Yhat_test, cal_method, alpha),
      error = function(e) list(theta = NA, var = NA, ci_lower = NA, ci_upper = NA)
    )

    method_name <- paste0("EIF-", cal_method)
    results[[length(results) + 1]] <- tibble(
      method = method_name,
      theta_hat = eif$theta,
      var_hat = eif$var,
      ci_lower = eif$ci_lower,
      ci_upper = eif$ci_upper,
      theta_true = theta_true,
      covers = if (!is.na(eif$ci_lower)) {
        theta_true >= eif$ci_lower & theta_true <= eif$ci_upper
      } else NA
    )
  }

  # Combine and add metadata
  bind_rows(results) %>%
    mutate(
      N_total = N_total,
      model_quality = model_quality,
      label_ratio = label_ratio,
      m = m,
      n_test = n_test,
      n_train = N_train
    )
}

# ============================================================================
# RUN FULL SIMULATION
# ============================================================================

#' Run full Monte Carlo simulation
run_continuous_simulation <- function(config = SIM_CONFIG, parallel = TRUE) {
  params <- expand_grid(
    N_total = config$N_total,
    model_quality = config$model_qualities,
    label_ratio = config$label_ratios,
    rep = 1:config$B
  )

  cat("\n")
  cat("=======================================================\n")
  cat("       CONTINUOUS SCORE SIMULATION\n")
  cat("=======================================================\n")
  cat(sprintf("Total simulations: %d\n", nrow(params)))
  cat(sprintf("  - Monte Carlo reps: %d\n", config$B))
  cat(sprintf("  - X ~ N(%.0f, 1) for ALL data (no covariate shift)\n", config$mu_x))
  cat(sprintf("  - N_total in {%s}\n", paste(config$N_total, collapse = ", ")))
  cat(sprintf("  - train_frac: %.0f%%\n", config$train_frac * 100))
  cat(sprintf("  - model_qualities: %s\n", paste(config$model_qualities, collapse = ", ")))
  cat(sprintf("  - label_ratios: %s\n", paste(config$label_ratios, collapse = ", ")))
  cat(sprintf("  - sigma_y=%.1f, alpha=%.2f\n", config$sigma_y, config$alpha))
  cat("=======================================================\n\n")

  if (parallel) {
    n_workers <- availableCores() - 1
    cat(sprintf("Running in PARALLEL mode with %d workers...\n\n", n_workers))
    plan(multisession, workers = n_workers)

    sim_file <- here("R/sim_continuous.R")

    results <- future_pmap_dfr(
      params,
      function(N_total, model_quality, label_ratio, rep) {
        source(sim_file, local = TRUE)

        sim_one_rep_continuous(
          N_total = N_total,
          model_quality = model_quality,
          label_ratio = label_ratio,
          train_frac = config$train_frac,
          mu_x = config$mu_x,
          sigma_y = config$sigma_y,
          alpha = config$alpha,
          seed = rep * 1000 + N_total + round(label_ratio * 100)
        ) %>%
          mutate(rep = rep)
      },
      .options = furrr_options(
        seed = TRUE,
        packages = c("mgcv", "splines", "randomForest", "rpart")
      ),
      .progress = TRUE
    )

    plan(sequential)
  } else {
    cat("Running in SEQUENTIAL mode...\n\n")
    results <- pmap_dfr(
      params,
      function(N_total, model_quality, label_ratio, rep) {
        if (rep %% 50 == 0) {
          cat(sprintf("  Progress: rep %d / %d\n", rep, config$B))
        }
        sim_one_rep_continuous(
          N_total = N_total,
          model_quality = model_quality,
          label_ratio = label_ratio,
          train_frac = config$train_frac,
          mu_x = config$mu_x,
          sigma_y = config$sigma_y,
          alpha = config$alpha,
          seed = rep * 1000 + N_total + round(label_ratio * 100)
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

summarize_continuous_results <- function(results) {
  results %>%
    filter(!is.na(theta_hat)) %>%
    group_by(method, N_total, model_quality, label_ratio) %>%
    summarize(
      n_reps = n(),
      mean_theta_true = mean(theta_true, na.rm = TRUE),  # actual E[Y] for this setting
      coverage = mean(covers, na.rm = TRUE),
      coverage_se = sqrt(coverage * (1 - coverage) / n_reps),
      mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
      sd_ci_width = sd(ci_upper - ci_lower, na.rm = TRUE),
      bias = mean(theta_hat - theta_true, na.rm = TRUE),
      rmse = sqrt(mean((theta_hat - theta_true)^2, na.rm = TRUE)),
      mean_var = mean(var_hat, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(coverage_gap = coverage - 0.90)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if (sys.nframe() == 0) {
  cat("\n")
  cat("##########################################################\n")
  cat("#  CONTINUOUS SCORE SIMULATION - LLM Judge Debiasing    #\n")
  cat("##########################################################\n")
  cat("\n")
  cat("DGP: X ~ N(1, 1), Y = f(X1,...,X5) + noise (no covariate shift!)\n")
  cat("     f has strong nonlinearities (X^2, sin, interactions, exp)\n")
  cat("Models:\n")
  cat("  - good:      GAM on X1-X5 (correct features, flexible)\n")
  cat("  - medium:    RF on X1-X10 (all features, some overfitting)\n")
  cat("  - poor-OLS:  Linear on X1-X10 (OLS guarantees E[Yhat]=E[Y] on train)\n")
  cat("  - poor-tree: Regression tree on X1,X6 (NOT OLS, no mean guarantee!)\n")
  cat("\n")

  start_time <- Sys.time()
  results <- run_continuous_simulation(parallel = TRUE)
  elapsed <- difftime(Sys.time(), start_time, units = "mins")

  cat(sprintf("\nSimulation completed in %.1f minutes\n", as.numeric(elapsed)))

  write_csv(results, "reproduce/continuous_sim_results.csv")
  cat("Raw results saved to: reproduce/continuous_sim_results.csv\n")

  summary_df <- summarize_continuous_results(results)
  write_csv(summary_df, "reproduce/continuous_sim_summary.csv")
  cat("Summary saved to: reproduce/continuous_sim_summary.csv\n")

  cat("\n")
  cat("=======================================================\n")
  cat("  RESULTS PREVIEW (label_ratio=0.10, N_total=2000)\n")
  cat("=======================================================\n")
  summary_df %>%
    filter(label_ratio == 0.10, N_total == 2000) %>%
    select(method, model_quality, bias, coverage, rmse) %>%
    arrange(model_quality, method) %>%
    print(n = 30)

  cat("\n")
  cat("##########################################################\n")
  cat("#                 SIMULATION COMPLETE                    #\n")
  cat("##########################################################\n")
}
