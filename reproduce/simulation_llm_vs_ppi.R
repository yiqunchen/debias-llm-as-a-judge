library(tidyverse)
library(future)
library(furrr)
library(debiasLLMReporting)

# respects Slurm/LSF/containers
plan(multisession, workers = availableCores())
set.seed(123)

# =============================================================================
# SIMULATION CONFIGURATION
# =============================================================================
#
# Sampling design (RANDOM, not stratified):
#   1. Generate N observations with true prevalence theta
#   2. Randomly select m = label_ratio * N for labeling (calibration)
#   3. m0 and m1 are RANDOM - determined by theta and sampling
#   4. Remaining n = N - m are unlabeled (test)
#
# All estimators use the SAME split:
#   - Calibration: m labeled observations (Y, Yhat pairs)
#   - Test: n unlabeled observations (Yhat only)
#
# =============================================================================

SIM_CONFIG <- list(
  full_mle = TRUE,
  N = 2000L,                                    # Total sample size
  label_ratios = c(0.01, 0.05, 0.10, 0.20),     # Fraction labeled
  B = 1000L,                                    # Monte Carlo replicates
  thetas = seq(0.1, 0.9, 0.1),          # True prevalence
  q0_vals = c(0.6, 0.7, 0.8),                        # Specificity values
  q1_vals = c(0.6, 0.7, 0.8),                        # Sensitivity values
  signal = 1,
  nominal_coverage = 0.9
)

# =============================================================================
# SIMULATION FUNCTIONS
# =============================================================================

#' Run one Monte Carlo replicate
#'
#' @param theta True prevalence
#' @param q0 True specificity P(Yhat=0|Y=0)
#' @param q1 True sensitivity P(Yhat=1|Y=1)
#' @param label_ratio Fraction of observations that are labeled
#' @param N Total sample size
#' @param signal DGP signal strength
#' @param rep_id Replicate ID
#' @param nominal Nominal coverage level
sim_one_rep <- function(theta, q0, q1, label_ratio, N, signal, rep_id, nominal) {

  # Empty result for early returns

  empty_row <- tibble(
    method = character(),
    theta_hat = numeric(),
    var_hat = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    lambda = numeric(),
    theta_true = numeric(),
    q0 = numeric(),
    q1 = numeric(),
    label_ratio = numeric(),
    N = integer(),
    n = integer(),
    m = integer(),
    m0 = integer(),
    m1 = integer(),
    rep = integer(),
    q0_hat = numeric(),
    q1_hat = numeric()
  )

  # ---------------------------------------------------------------------------
  # Generate data
  # ---------------------------------------------------------------------------
  dgp <- generate_dgp_data(N, theta, signal)
  Y <- dgp$Y

  # Surrogate predictions with sensitivity q1 and specificity q0
  Yhat <- ifelse(
    Y == 1,
    rbinom(N, 1, q1),        # True positive rate = q1
    rbinom(N, 1, 1 - q0)     # False positive rate = 1 - q0
  )

  # ---------------------------------------------------------------------------
  # Random split into calibration (labeled) and test (unlabeled)
  # ---------------------------------------------------------------------------
  m <- max(4L, round(label_ratio * N))   # Calibration size
  n <- N - m                              # Test size

  idx_cal <- sample(N, m)                 # Random calibration indices

idx_test <- setdiff(seq_len(N), idx_cal)

  # Calibration data (labeled)
  Y_cal <- Y[idx_cal]
  Yhat_cal <- Yhat[idx_cal]

  # Test data (unlabeled - only Yhat observed)
  Yhat_test <- Yhat[idx_test]

  # Observed class counts in calibration (RANDOM, depends on theta)
  m0 <- sum(Y_cal == 0)
  m1 <- sum(Y_cal == 1)

  # Need at least 2 of each class for estimation
  if (m0 < 2 || m1 < 2) {
    return(empty_row)
  }

  # ---------------------------------------------------------------------------
  # Estimate confusion matrix parameters from calibration
  # ---------------------------------------------------------------------------
  q0_hat <- mean(Yhat_cal[Y_cal == 0] == 0)  # Estimated specificity
  q1_hat <- mean(Yhat_cal[Y_cal == 1] == 1)  # Estimated sensitivity
  p_hat <- mean(Yhat_test)                    # Test positive rate

  # ---------------------------------------------------------------------------
  # ESTIMATORS
  # All use same calibration (m) and test (n) split
  # ---------------------------------------------------------------------------
  results <- list()

  # 1. PPI estimator
  ppi_est <- ppi_point_and_ci(Y_L = Y_cal, f_L = Yhat_cal, f_U = Yhat_test,
                               alpha = 1 - nominal)
  results$ppi <- tibble(
    method = "PPI",
    theta_hat = ppi_est$theta,
    var_hat = ppi_est$var,
    ci_lower = ppi_est$ci_lower,
    ci_upper = ppi_est$ci_upper,
    lambda = NA_real_
  )

  # 2. PPI++ estimator
  ppi_pp_est <- ppi_pp_point_and_ci_general(Y_L = Y_cal, f_L = Yhat_cal,
                                             f_U = Yhat_test, alpha = 1 - nominal)
  results$ppi_pp <- tibble(
    method = "PPI++",
    theta_hat = ppi_pp_est$theta,
    var_hat = ppi_pp_est$var,
    ci_lower = ppi_pp_est$ci_lower,
    ci_upper = ppi_pp_est$ci_upper,
    lambda = ppi_pp_est$lambda
  )

  # 3. Rogan-Gladen estimator
  rg_est <- llm_point_and_ci(
    p_hat = p_hat,
    q0_hat = q0_hat,
    q1_hat = q1_hat,
    n = n,
    m0 = m0,
    m1 = m1,
    alpha = 1 - nominal
  )
  results$rg <- tibble(
    method = "Rogan-Gladen",
    theta_hat = rg_est$theta,
    var_hat = rg_est$var,
    ci_lower = rg_est$ci_lower,
    ci_upper = rg_est$ci_upper,
    lambda = NA_real_
  )

  # 4. Joint MLE estimator
  if (isTRUE(SIM_CONFIG$full_mle)) {
    mle_result <- tryCatch({
      mle_fit <- fit_misclass_mle(
        y_cal = Y_cal,
        yhat_cal = Yhat_cal,
        yhat_test = Yhat_test,
        level = nominal
      )
      se_theta <- unname(mle_fit$se_obs["theta"])
      var_hat <- if (is.na(se_theta)) NA_real_ else se_theta^2
      tibble(
        method = "Joint MLE",
        theta_hat = mle_fit$theta_hat,
        var_hat = var_hat,
        ci_lower = mle_fit$ci_theta_obs[1],
        ci_upper = mle_fit$ci_theta_obs[2],
        lambda = NA_real_
      )
    }, error = function(e) {
      tibble(
        method = "Joint MLE",
        theta_hat = NA_real_,
        var_hat = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        lambda = NA_real_
      )
    })
    results$mle <- mle_result
  }

  # 5. Naive estimator (no correction)
  naive_theta <- p_hat
  naive_var <- p_hat * (1 - p_hat) / n
  z_alpha <- qnorm(1 - (1 - nominal) / 2)
  naive_se <- sqrt(naive_var)
  results$naive <- tibble(
    method = "Naive",
    theta_hat = naive_theta,
    var_hat = naive_var,
    ci_lower = pmax(naive_theta - z_alpha * naive_se, 0),
    ci_upper = pmin(naive_theta + z_alpha * naive_se, 1),
    lambda = NA_real_
  )

  # 6. EIF estimator
  eif_est <- eif_point_and_ci(
    Y_cal = Y_cal,
    Yhat_cal = Yhat_cal,
    Yhat_test = Yhat_test,
    alpha = 1 - nominal
  )
  results$eif <- tibble(
    method = "EIF",
    theta_hat = eif_est$theta,
    var_hat = eif_est$var,
    ci_lower = eif_est$ci_lower,
    ci_upper = eif_est$ci_upper,
    lambda = NA_real_
  )

  # Combine all results
  bind_rows(results) %>%
    mutate(
      theta_true = theta,
      q0 = q0,
      q1 = q1,
      label_ratio = label_ratio,
      N = N,
      n = n,
      m = m,
      m0 = m0,
      m1 = m1,
      rep = rep_id,
      q0_hat = q0_hat,
      q1_hat = q1_hat
    )
}

#' Run B replicates for one parameter configuration
run_one_config <- function(theta, q0, q1, label_ratio) {
  map_dfr(
    seq_len(SIM_CONFIG$B),
    ~ sim_one_rep(
      theta = theta,
      q0 = q0,
      q1 = q1,
      label_ratio = label_ratio,
      N = SIM_CONFIG$N,
      signal = SIM_CONFIG$signal,
      rep_id = .x,
      nominal = SIM_CONFIG$nominal_coverage
    )
  )
}

# =============================================================================
# RUN SIMULATION
# =============================================================================

# Parameter grid (no more calibration_balances - m0/m1 are random!)
param_grid <- expand_grid(
  theta = SIM_CONFIG$thetas,
  q0 = SIM_CONFIG$q0_vals,
  q1 = SIM_CONFIG$q1_vals,
  label_ratio = SIM_CONFIG$label_ratios
)

message("Running simulation with ", nrow(param_grid), " configurations x ",
        SIM_CONFIG$B, " replicates")

# Setup output directory
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- file.path("results", timestamp)
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

# Run simulation in parallel
res_all <- future_pmap_dfr(
  list(
    theta = param_grid$theta,
    q0 = param_grid$q0,
    q1 = param_grid$q1,
    label_ratio = param_grid$label_ratio
  ),
  run_one_config,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

# =============================================================================
# SUMMARIZE RESULTS
# =============================================================================

# CI and bias summary
ci_summary <- res_all %>%
  filter(!is.na(theta_hat)) %>%
  mutate(covered = theta_true >= ci_lower & theta_true <= ci_upper) %>%
  group_by(theta_true, q0, q1, label_ratio, N, method) %>%
  summarise(
    n_reps = n(),
    ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    mse = mean((theta_hat - theta_true)^2, na.rm = TRUE),
    bias = mean(theta_hat - theta_true, na.rm = TRUE),
    var_est = var(theta_hat, na.rm = TRUE),
    lambda_mean = mean(lambda, na.rm = TRUE),
    m0_mean = mean(m0, na.rm = TRUE),
    m1_mean = mean(m1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(coverage_gap = coverage - SIM_CONFIG$nominal_coverage)

readr::write_csv(ci_summary, file.path(results_dir, "ci_summary.csv"))

# Calibration parameter estimation summary
q_est_summary <- res_all %>%
  distinct(theta_true, q0, q1, label_ratio, N, rep, q0_hat, q1_hat, m0, m1) %>%
  group_by(theta_true, q0, q1, label_ratio, N) %>%
  summarise(
    q0_bias = mean(q0_hat - q0, na.rm = TRUE),
    q1_bias = mean(q1_hat - q1, na.rm = TRUE),
    q0_rmse = sqrt(mean((q0_hat - q0)^2, na.rm = TRUE)),
    q1_rmse = sqrt(mean((q1_hat - q1)^2, na.rm = TRUE)),
    m0_mean = mean(m0, na.rm = TRUE),
    m1_mean = mean(m1, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(q_est_summary, file.path(results_dir, "q_est_summary.csv"))

# =============================================================================
# PLOTTING
# =============================================================================

save_plot <- function(plot_obj, filename, width = 14, height = 12) {
  ggsave(file.path(plots_dir, filename), plot_obj, width = width, height = height, dpi = 400)
}

# Prepare plotting data
ci_summary_plot <- ci_summary %>%
  mutate(
    bias_pct = if_else(theta_true != 0, 100 * bias / theta_true, NA_real_),
    q0_lab = paste0("q0 = ", q0),
    q1_lab = paste0("q1 = ", q1),
    label_ratio_lab = factor(format_label_ratio(label_ratio),
                              levels = format_label_ratio(sort(unique(label_ratio))))
  )

q_bias_long <- q_est_summary %>%
  pivot_longer(cols = c(q0_bias, q1_bias), names_to = "metric", values_to = "bias") %>%
  mutate(
    metric = recode(metric, q0_bias = "q0", q1_bias = "q1"),
    q0_lab = paste0("q0 = ", q0),
    q1_lab = paste0("q1 = ", q1),
    label_ratio_lab = factor(format_label_ratio(label_ratio),
                              levels = format_label_ratio(sort(unique(label_ratio))))
  )

# -----------------------------------------------------------------------------
# GRID PLOTS: Per label_ratio, faceted by q0 x q1
# -----------------------------------------------------------------------------

for (lr in unique(ci_summary_plot$label_ratio)) {
  suffix <- paste0("_labeled_", sprintf("%02d", as.integer(lr * 100)), ".png")
  subtitle <- paste0(scales::percent(lr, accuracy = 1), " labeled (m0, m1 random)")

  ci_sub <- filter(ci_summary_plot, label_ratio == lr)
  q_sub <- filter(q_bias_long, label_ratio == lr)

  message("Plotting: lr=", lr, " n=", nrow(ci_sub))

  # All methods
  save_plot(plot_coverage(ci_sub, q0_lab ~ q1_lab, "Coverage", subtitle,
                          x_breaks = SIM_CONFIG$thetas), paste0("coverage", suffix))
  save_plot(plot_ci_width(ci_sub, q0_lab ~ q1_lab, "CI Width", subtitle,
                          x_breaks = SIM_CONFIG$thetas), paste0("ci_width", suffix))
  save_plot(plot_bias(ci_sub, q0_lab ~ q1_lab, "Estimator Bias", subtitle,
                      x_breaks = SIM_CONFIG$thetas), paste0("bias", suffix))
  save_plot(plot_bias_pct(ci_sub, q0_lab ~ q1_lab, "Percent Bias", subtitle,
                          x_breaks = SIM_CONFIG$thetas), paste0("bias_pct", suffix))
  save_plot(plot_q_bias(q_sub, q0_lab ~ q1_lab, "Calibration Estimate Bias", subtitle,
                        x_breaks = SIM_CONFIG$thetas), paste0("q_bias", suffix))

  # Excluding Naive for better zoom
  ci_sub_no_naive <- filter(ci_sub, method != "Naive")
  save_plot(plot_coverage(ci_sub_no_naive, q0_lab ~ q1_lab, "Coverage (excl. Naive)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, y_limits = c(0.8, 1)), paste0("coverage_no_naive", suffix))
  save_plot(plot_ci_width(ci_sub_no_naive, q0_lab ~ q1_lab, "CI Width (excl. Naive)", subtitle,
                          x_breaks = SIM_CONFIG$thetas), paste0("ci_width_no_naive", suffix))
  save_plot(plot_bias(ci_sub_no_naive, q0_lab ~ q1_lab, "Estimator Bias (excl. Naive)", subtitle,
                      x_breaks = SIM_CONFIG$thetas), paste0("bias_no_naive", suffix))
  save_plot(plot_bias_pct(ci_sub_no_naive, q0_lab ~ q1_lab, "Percent Bias (excl. Naive)", subtitle,
                          x_breaks = SIM_CONFIG$thetas), paste0("bias_pct_no_naive", suffix))
}

# -----------------------------------------------------------------------------
# DIAGONAL PLOTS: q0 == q1, faceted by q_val x label_ratio
# -----------------------------------------------------------------------------

ci_summary_diag <- ci_summary_plot %>%
  filter(q0 == q1) %>%
  mutate(q_val_lab = paste0("q0 = q1 = ", q0))

q_bias_diag <- q_bias_long %>%
  filter(q0 == q1) %>%
  mutate(q_val_lab = paste0("q0 = q1 = ", q0))

if (nrow(ci_summary_diag) > 0) {
  subtitle <- "Diagonal (q0 = q1), m0/m1 random"

  # All methods
  save_plot(plot_coverage(ci_summary_diag, q_val_lab ~ label_ratio_lab, "Coverage (q0 = q1)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2), "coverage_diag.png", 16, 10)
  save_plot(plot_ci_width(ci_summary_diag, q_val_lab ~ label_ratio_lab, "CI Width (q0 = q1)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2), "ci_width_diag.png", 16, 10)
  save_plot(plot_bias(ci_summary_diag, q_val_lab ~ label_ratio_lab, "Estimator Bias (q0 = q1)", subtitle,
                      x_breaks = SIM_CONFIG$thetas, point_size = 2), "bias_diag.png", 16, 10)
  save_plot(plot_bias_pct(ci_summary_diag, q_val_lab ~ label_ratio_lab, "Percent Bias (q0 = q1)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2), "bias_pct_diag.png", 16, 10)
  if (nrow(q_bias_diag) > 0) {
    save_plot(plot_q_bias(q_bias_diag, q_val_lab ~ label_ratio_lab, "Calibration Estimate Bias (q0 = q1)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2), "q_bias_diag.png", 16, 10)
  }

  # Excluding Naive
  ci_diag_no_naive <- filter(ci_summary_diag, method != "Naive")
  save_plot(plot_coverage(ci_diag_no_naive, q_val_lab ~ label_ratio_lab, "Coverage (q0 = q1, excl. Naive)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2, y_limits = c(0.8, 1)), "coverage_no_naive_diag.png", 16, 10)
  save_plot(plot_ci_width(ci_diag_no_naive, q_val_lab ~ label_ratio_lab, "CI Width (q0 = q1, excl. Naive)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2), "ci_width_no_naive_diag.png", 16, 10)
  save_plot(plot_bias(ci_diag_no_naive, q_val_lab ~ label_ratio_lab, "Estimator Bias (q0 = q1, excl. Naive)", subtitle,
                      x_breaks = SIM_CONFIG$thetas, point_size = 2), "bias_no_naive_diag.png", 16, 10)
  save_plot(plot_bias_pct(ci_diag_no_naive, q_val_lab ~ label_ratio_lab, "Percent Bias (q0 = q1, excl. Naive)", subtitle,
                          x_breaks = SIM_CONFIG$thetas, point_size = 2), "bias_pct_no_naive_diag.png", 16, 10)
}

message("Results saved to ", results_dir)
