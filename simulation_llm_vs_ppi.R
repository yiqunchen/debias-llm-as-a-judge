library(tidyverse)
library(future)
library(furrr)
library(debiasLLMReporting)

source("R/plot_utils.R")

plan(multisession, workers = max(1, parallel::detectCores() - 1))
set.seed(123)

SIM_CONFIG <- list(
  sample_size = 5000L,
  label_ratios = c(0.01, 0.05, 0.10, 0.20, 0.50),
  B = 1000L,
  thetas = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  q0_vals = c(0.6, 0.7, 0.8, 0.9),
  q1_vals = c(0.6, 0.7, 0.8, 0.9),
  signal = 1,
  nominal_coverage = 0.95,
  calibration_balances = tibble(
    split_label = c("20:80", "50:50", "80:20"),
    m0_prop = c(0.2, 0.5, 0.8)
  )
)

format_label_ratio <- function(x) paste0(scales::percent(x, accuracy = 1), " labeled")

g_specs <- list(
  list(label = "PPI", transform = function(scores, context) scores)
)

derive_label_counts <- function(label_ratio, sample_size, m0_prop) {
  total_labels <- max(4L, round(label_ratio * sample_size))
  m0 <- max(2L, round(total_labels * m0_prop))
  m1 <- max(2L, total_labels - m0)
  if (m0 + m1 > total_labels) {
    m1 <- total_labels - m0
  }
  list(
    n_ppi = total_labels,
    N_ppi = sample_size,
    n_llm_test = sample_size,
    m0 = m0,
    m1 = m1
  )
}

sim_one_condition <- function(theta, q0, q1, label_ratio, sample_size,
                              signal, g_specs, rep_id, nominal,
                              m0_prop, split_label) {
  counts <- derive_label_counts(label_ratio, sample_size, m0_prop)
  N_required <- max(counts$n_ppi + counts$N_ppi,
                    counts$n_llm_test + counts$m0 + counts$m1)
  dgp <- generate_dgp_data(N_required, theta, signal)
  Y <- dgp$Y
  Z <- ifelse(
    Y == 1,
    rbinom(N_required, 1, q1),
    rbinom(N_required, 1, 1 - q0)
  )
  f_hat <- Z
  Y_L <- Y[seq_len(counts$n_ppi)]
  f_L <- f_hat[seq_len(counts$n_ppi)]
  f_U <- f_hat[(counts$n_ppi + 1):(counts$n_ppi + counts$N_ppi)]
  Z_test <- f_hat[seq_len(counts$n_llm_test)]
  p_hat <- mean(Z_test)
  idx_neg <- which(Y == 0)
  idx_pos <- which(Y == 1)
  if (length(idx_neg) < counts$m0 || length(idx_pos) < counts$m1) {
    return(tibble())
  }
  Z_cal0 <- Z[sample(idx_neg, counts$m0)]
  Z_cal1 <- Z[sample(idx_pos, counts$m1)]
  q0_hat <- mean(Z_cal0 == 0)
  q1_hat <- mean(Z_cal1 == 1)
  theta_pilot <- mean(Y_L)
  g_context <- list(q0_hat = q0_hat, q1_hat = q1_hat, theta_pilot = theta_pilot)
  ppi_results <- map_dfr(g_specs, function(spec) {
    f_L_g <- spec$transform(f_L, g_context)
    f_U_g <- spec$transform(f_U, g_context)
    est <- ppi_point_and_ci(Y_L, f_L_g, f_U_g)
    tibble(
      method = spec$label,
      theta_hat = est$theta,
      var_hat = est$var,
      ci_lower = est$ci_lower,
      ci_upper = est$ci_upper,
      lambda = NA_real_
    )
  })
  ppi_pp_est <- ppi_pp_point_and_ci_general(Y_L = Y_L, f_L = f_L, f_U = f_U)
  ppi_pp_tbl <- tibble(
    method = "PPI++",
    theta_hat = ppi_pp_est$theta,
    var_hat = ppi_pp_est$var,
    ci_lower = ppi_pp_est$ci_lower,
    ci_upper = ppi_pp_est$ci_upper,
    lambda = ppi_pp_est$lambda
  )
  llm_est <- llm_point_and_ci(
    p_hat = p_hat,
    q0_hat = q0_hat,
    q1_hat = q1_hat,
    n = counts$n_llm_test,
    m0 = counts$m0,
    m1 = counts$m1,
    alpha = 1 - nominal
  )
  llm_tbl <- tibble(
    method = "Rogan-Gladen",
    theta_hat = llm_est$theta,
    var_hat = llm_est$var,
    ci_lower = llm_est$ci_lower,
    ci_upper = llm_est$ci_upper,
    lambda = NA_real_
  )
  # Naive estimator: just use raw LLM predictions with standard Wald CI
  n_naive <- counts$n_llm_test
  naive_theta <- p_hat
  naive_var <- p_hat * (1 - p_hat) / n_naive
  z_alpha <- qnorm(1 - (1 - nominal) / 2)
  naive_se <- sqrt(naive_var)
  naive_ci_lower <- pmax(naive_theta - z_alpha * naive_se, 0)
  naive_ci_upper <- pmin(naive_theta + z_alpha * naive_se, 1)
  naive_tbl <- tibble(
    method = "Naive",
    theta_hat = naive_theta,
    var_hat = naive_var,
    ci_lower = naive_ci_lower,
    ci_upper = naive_ci_upper,
    lambda = NA_real_
  )
  bind_rows(ppi_results, ppi_pp_tbl, llm_tbl, naive_tbl) %>%
    mutate(
      theta_true = theta,
      q0 = q0,
      q1 = q1,
      label_ratio = label_ratio,
      sample_size = sample_size,
      rep = rep_id,
      q0_hat = q0_hat,
      q1_hat = q1_hat,
      n_labeled = counts$n_ppi,
      m0 = counts$m0,
      m1 = counts$m1,
      split_label = split_label,
      m0_prop = m0_prop,
      m1_prop = 1 - m0_prop
    )
}

run_one_grid <- function(theta, q0, q1, label_ratio, m0_prop, split_label) {
  map_dfr(
    1:SIM_CONFIG$B,
    ~ sim_one_condition(
      theta = theta,
      q0 = q0,
      q1 = q1,
      label_ratio = label_ratio,
      sample_size = SIM_CONFIG$sample_size,
      signal = SIM_CONFIG$signal,
      g_specs = g_specs,
      rep_id = .x,
      nominal = SIM_CONFIG$nominal_coverage,
      m0_prop = m0_prop,
      split_label = split_label
    )
  )
}

param_grid <- expand_grid(
  theta = SIM_CONFIG$thetas,
  q0 = SIM_CONFIG$q0_vals,
  q1 = SIM_CONFIG$q1_vals,
  label_ratio = SIM_CONFIG$label_ratios,
  split_label = SIM_CONFIG$calibration_balances$split_label
) %>%
  left_join(SIM_CONFIG$calibration_balances, by = "split_label")

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- file.path("results", timestamp)
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

res_all <- future_pmap_dfr(
  list(
    theta = param_grid$theta,
    q0 = param_grid$q0,
    q1 = param_grid$q1,
    label_ratio = param_grid$label_ratio,
    m0_prop = param_grid$m0_prop,
    split_label = param_grid$split_label
  ),
  run_one_grid,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

readr::write_csv(res_all, file.path(results_dir, "raw_simulation_results.csv"))

ci_summary <- res_all %>%
  filter(!is.na(theta_hat)) %>%
  mutate(
    covered = theta_true >= ci_lower & theta_true <= ci_upper
  ) %>%
  group_by(theta_true, q0, q1, label_ratio, split_label, sample_size, method) %>%
  summarise(
    ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    mse = mean((theta_hat - theta_true)^2, na.rm = TRUE),
    bias = mean(theta_hat - theta_true, na.rm = TRUE),
    var_est = var(theta_hat, na.rm = TRUE),
    lambda_mean = mean(lambda, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_gap = coverage - SIM_CONFIG$nominal_coverage
  )

readr::write_csv(ci_summary, file.path(results_dir, "ci_summary.csv"))

q_est_summary <- res_all %>%
  distinct(theta_true, q0, q1, label_ratio, split_label, sample_size, rep, q0_hat, q1_hat) %>%
  group_by(theta_true, q0, q1, label_ratio, split_label, sample_size) %>%
  summarise(
    q0_bias = mean(q0_hat - q0, na.rm = TRUE),
    q1_bias = mean(q1_hat - q1, na.rm = TRUE),
    q0_rmse = sqrt(mean((q0_hat - q0)^2, na.rm = TRUE)),
    q1_rmse = sqrt(mean((q1_hat - q1)^2, na.rm = TRUE)),
    .groups = "drop"
  )

readr::write_csv(q_est_summary, file.path(results_dir, "q_est_summary.csv"))

# Save plot helper (wraps ggsave with plots_dir)
save_plot <- function(plot_obj, filename, width = 14, height = 12) {
  ggsave(file.path(plots_dir, filename), plot_obj, width = width, height = height, dpi = 400)
}

# =============================================================================
# PLOTTING SECTION
# =============================================================================

# Prepare data for plotting
ci_summary_plot <- ci_summary %>%
  mutate(
    split_label = recode(split_label, "q0_20" = "20:80", "q0_50" = "50:50", "q0_80" = "80:20"),
    bias_pct = dplyr::if_else(theta_true != 0, 100 * bias / theta_true, NA_real_),
    q0_lab = paste0("q0 = ", q0),
    q1_lab = paste0("q1 = ", q1)
  )

q_bias_long <- q_est_summary %>%
  mutate(split_label = recode(split_label, "q0_20" = "20:80", "q0_50" = "50:50", "q0_80" = "80:20")) %>%
  pivot_longer(cols = c(q0_bias, q1_bias), names_to = "metric", values_to = "bias") %>%
  mutate(
    metric = recode(metric, q0_bias = "q0", q1_bias = "q1"),
    q0_lab = paste0("q0 = ", q0),
    q1_lab = paste0("q1 = ", q1)
  )

# -----------------------------------------------------------------------------
# GRID PLOTS: One plot per (split, label_ratio), faceted by q0 x q1
# -----------------------------------------------------------------------------
for (split in unique(ci_summary_plot$split_label)) {
  for (lr in unique(ci_summary_plot$label_ratio)) {
    suffix <- paste0("_", split, "_labeled_", sprintf("%02d", as.integer(lr * 100)), ".png")
    subtitle <- paste0("neg:pos = ", split, " | ", scales::percent(lr, accuracy = 1), " labeled")

    ci_sub <- filter(ci_summary_plot, split_label == split, label_ratio == lr)
    q_sub <- filter(q_bias_long, split_label == split, label_ratio == lr)

    save_plot(plot_coverage(ci_sub, q0_lab ~ q1_lab, "Coverage", subtitle, x_breaks = SIM_CONFIG$thetas), paste0("coverage", suffix))
    save_plot(plot_ci_width(ci_sub, q0_lab ~ q1_lab, "CI Width", subtitle, x_breaks = SIM_CONFIG$thetas), paste0("ci_width", suffix))
    save_plot(plot_bias(ci_sub, q0_lab ~ q1_lab, "Estimator Bias", subtitle, x_breaks = SIM_CONFIG$thetas), paste0("bias", suffix))
    save_plot(plot_bias_pct(ci_sub, q0_lab ~ q1_lab, "Percent Bias", subtitle, x_breaks = SIM_CONFIG$thetas), paste0("bias_pct", suffix))
    save_plot(plot_q_bias(q_sub, q0_lab ~ q1_lab, "Calibration Estimate Bias", subtitle, x_breaks = SIM_CONFIG$thetas), paste0("q_bias", suffix))
  }
}

# -----------------------------------------------------------------------------
# DIAGONAL PLOTS: q0 == q1, faceted by q_val x label_ratio
# -----------------------------------------------------------------------------
ci_summary_diag <- ci_summary_plot %>%
  filter(q0 == q1) %>%
  mutate(
    q_val_lab = paste0("q0 = q1 = ", q0),
    label_ratio_lab = factor(format_label_ratio(label_ratio),
                             levels = format_label_ratio(sort(unique(label_ratio))))
  )

q_bias_diag <- q_bias_long %>%
  filter(q0 == q1) %>%
  mutate(
    q_val_lab = paste0("q0 = q1 = ", q0),
    label_ratio_lab = factor(format_label_ratio(label_ratio),
                             levels = format_label_ratio(sort(unique(label_ratio))))
  )

# Per-split diagonal plots
for (split in unique(ci_summary_diag$split_label)) {
  suffix <- paste0("_diag_", split, ".png")
  subtitle <- paste0("neg:pos = ", split, " | q0 = q1 (diagonal)")

  ci_sub <- filter(ci_summary_diag, split_label == split)
  q_sub <- filter(q_bias_diag, split_label == split)

  save_plot(plot_coverage(ci_sub, q_val_lab ~ label_ratio_lab, "Coverage (q0 = q1)", subtitle, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("coverage", suffix), 16, 10)
  save_plot(plot_ci_width(ci_sub, q_val_lab ~ label_ratio_lab, "CI Width (q0 = q1)", subtitle, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("ci_width", suffix), 16, 10)
  save_plot(plot_bias(ci_sub, q_val_lab ~ label_ratio_lab, "Estimator Bias (q0 = q1)", subtitle, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("bias", suffix), 16, 10)
  save_plot(plot_bias_pct(ci_sub, q_val_lab ~ label_ratio_lab, "Percent Bias (q0 = q1)", subtitle, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("bias_pct", suffix), 16, 10)
  save_plot(plot_q_bias(q_sub, q_val_lab ~ label_ratio_lab, "Calibration Estimate Bias (q0 = q1)", subtitle, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("q_bias", suffix), 16, 10)
}

# -----------------------------------------------------------------------------
# AGGREGATE DIAGONAL PLOTS: All neg:pos splits with linetype
# -----------------------------------------------------------------------------
save_plot(plot_coverage(ci_summary_diag, q_val_lab ~ label_ratio_lab, "Coverage (q0 = q1)", linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), "coverage_diag_aggregate.png", 16, 10)
save_plot(plot_ci_width(ci_summary_diag, q_val_lab ~ label_ratio_lab, "CI Width (q0 = q1)", linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), "ci_width_diag_aggregate.png", 16, 10)
save_plot(plot_bias(ci_summary_diag, q_val_lab ~ label_ratio_lab, "Estimator Bias (q0 = q1)", linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), "bias_diag_aggregate.png", 16, 10)
save_plot(plot_bias_pct(ci_summary_diag, q_val_lab ~ label_ratio_lab, "Percent Bias (q0 = q1)", linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), "bias_pct_diag_aggregate.png", 16, 10)
save_plot(plot_q_bias(q_bias_diag, q_val_lab ~ label_ratio_lab, "Calibration Estimate Bias (q0 = q1)", linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), "q_bias_diag_aggregate.png", 16, 10)

# -----------------------------------------------------------------------------
# OFF-DIAGONAL PLOTS: (q0=0.6, q1=0.9) and (q0=0.9, q1=0.6)
# -----------------------------------------------------------------------------
off_diag_cases <- list(
  list(q0 = 0.6, q1 = 0.9, label = "q0=0.6_q1=0.9"),
  list(q0 = 0.9, q1 = 0.6, label = "q0=0.9_q1=0.6")
)

for (case in off_diag_cases) {
  ci_sub <- ci_summary_plot %>%
    filter(q0 == case$q0, q1 == case$q1) %>%
    mutate(label_ratio_lab = factor(format_label_ratio(label_ratio),
                                    levels = format_label_ratio(sort(unique(label_ratio)))))

  q_sub <- q_bias_long %>%
    filter(q0 == case$q0, q1 == case$q1) %>%
    mutate(label_ratio_lab = factor(format_label_ratio(label_ratio),
                                    levels = format_label_ratio(sort(unique(label_ratio)))))

  suffix <- paste0("_offdiag_", case$label, ".png")
  subtitle <- paste0("q0 = ", case$q0, ", q1 = ", case$q1, " | All neg:pos ratios")
  title_prefix <- paste0(" (q0 = ", case$q0, ", q1 = ", case$q1, ")")

  save_plot(plot_coverage(ci_sub, ~ label_ratio_lab, paste0("Coverage", title_prefix), subtitle, linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("coverage", suffix), 18, 6)
  save_plot(plot_ci_width(ci_sub, ~ label_ratio_lab, paste0("CI Width", title_prefix), subtitle, linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("ci_width", suffix), 18, 6)
  save_plot(plot_bias(ci_sub, ~ label_ratio_lab, paste0("Estimator Bias", title_prefix), subtitle, linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("bias", suffix), 18, 6)
  save_plot(plot_bias_pct(ci_sub, ~ label_ratio_lab, paste0("Percent Bias", title_prefix), subtitle, linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("bias_pct", suffix), 18, 6)
  save_plot(plot_q_bias(q_sub, ~ label_ratio_lab, paste0("Calibration Estimate Bias", title_prefix), subtitle, linetype_var = split_label, use_aggregate = TRUE, x_breaks = SIM_CONFIG$thetas, point_size = 2), paste0("q_bias", suffix), 18, 6)
}

message("Results saved to ", results_dir)
