library(tidyverse)
library(future)
library(furrr)

# Source the multinomial estimators
source("R/multinomial_estimators.R")
source("R/plot_utils.R")

plan(multisession, workers = max(1, parallel::detectCores() - 1))
set.seed(123)

# =============================================================================
# Simulation Configuration for K=5
# =============================================================================

K <- 5L

SIM_CONFIG <- list(
  K = K,
  n_test = 2000L,
  # Calibration ratios (like binary case)
  m_cal_ratios = c(0.01, 0.05, 0.10, 0.20),
  B = 1000L,
  nominal_coverage = 0.95,

  # Different true prevalence vectors (varying "theta" equivalent)
  # We parameterize by pi_1 (first class prevalence) similar to binary theta
  pi1_vals = c(0.1, 0.3, 0.5, 0.7, 0.9),

  # Confusion matrix accuracy (like q0/q1 in binary)
  accuracy_vals = c(0.6, 0.8),

  # Confusion types
  confusion_types = c("diagonal", "asymmetric")
)

format_cal_ratio <- function(x) paste0(scales::percent(x, accuracy = 1), " labeled")

# =============================================================================
# Data Generation Functions
# =============================================================================

#' Generate true prevalence vector given pi_1
#' Remaining mass split between other classes
generate_pi_true <- function(pi1, K) {
  remaining <- 1 - pi1
  # For K=5: class 2 gets 40%, class 3 gets 30%, classes 4-5 get 15% each
  if (K == 5) {
    pi <- c(pi1, remaining * 0.40, remaining * 0.30, remaining * 0.15, remaining * 0.15)
  } else {
    # Generic fallback
    pi <- c(pi1, rep(remaining / (K - 1), K - 1))
  }
  pi / sum(pi)  # ensure sums to 1
}

#' Generate confusion matrix
generate_confusion_matrix <- function(K, type = "diagonal", accuracy = 0.7) {
  if (type == "diagonal") {
    # Diagonal-dominant: correct class with prob `accuracy`, uniform over others
    M <- matrix((1 - accuracy) / (K - 1), nrow = K, ncol = K)
    diag(M) <- accuracy
  } else if (type == "asymmetric") {
    # Adjacent classes get confused more (ordinal-style errors)
    M <- matrix(0.02, nrow = K, ncol = K)
    for (k in 1:K) {
      M[k, k] <- accuracy
      remaining <- 1 - accuracy - 0.02 * (K - 1)
      if (k > 1) M[k-1, k] <- remaining * 0.7
      if (k < K) M[k+1, k] <- remaining * 0.3
      if (k == 1) M[2, k] <- remaining
      if (k == K) M[K-1, k] <- remaining
    }
    # Normalize columns
    M <- sweep(M, 2, colSums(M), "/")
  }
  M
}

#' Generate simulation data
generate_data <- function(n_test, m_cal, pi_true, M_true) {
  K <- length(pi_true)
  N <- n_test + m_cal

  Y <- sample(1:K, N, replace = TRUE, prob = pi_true)
  Yhat <- sapply(Y, function(y) sample(1:K, 1, prob = M_true[, y]))

  # Random split (not stratified!)
  idx_cal <- sample(1:N, m_cal)
  idx_test <- setdiff(1:N, idx_cal)

  list(
    yhat_test = Yhat[idx_test],
    y_cal = Y[idx_cal],
    yhat_cal = Yhat[idx_cal],
    y_test = Y[idx_test]
  )
}

# =============================================================================
# Single Simulation Replicate
# =============================================================================

run_one_sim <- function(pi1, accuracy, m_cal_ratio, confusion_type, rep_id) {
  K <- SIM_CONFIG$K
  n_test <- SIM_CONFIG$n_test
  m_cal <- round(n_test * m_cal_ratio)

  pi_true <- generate_pi_true(pi1, K)
  M_true <- generate_confusion_matrix(K, confusion_type, accuracy)

  data <- generate_data(n_test, m_cal, pi_true, M_true)

  # Prep counts for efficiency
  pre <- prep_counts(
    y_hat_test = data$yhat_test,
    y_true_cal = data$y_cal,
    y_hat_cal = data$yhat_cal,
    K = K
  )

  z <- qnorm(1 - (1 - SIM_CONFIG$nominal_coverage) / 2)

  results <- list()

  # 1. Naive
  naive_est <- tryCatch({
    res <- naive_multinomial(T = pre$T, K = K)
    tibble(
      method = "Naive", class = 1:K, pi_true = pi_true,
      pi_hat = res$pi_hat, se = res$se,
      ci_lower = pmax(res$pi_hat - z * res$se, 0),
      ci_upper = pmin(res$pi_hat + z * res$se, 1)
    )
  }, error = function(e) NULL)
  if (!is.null(naive_est)) results$Naive <- naive_est

  # 2. PPI (onehot)
  ppi_est <- tryCatch({
    res <- multiclass_ppi(T = pre$T, N = pre$N, K = K, g = "onehot", project = TRUE)
    tibble(
      method = "PPI", class = 1:K, pi_true = pi_true,
      pi_hat = res$pi_hat, se = res$se,
      ci_lower = pmax(res$pi_hat - z * res$se, 0),
      ci_upper = pmin(res$pi_hat + z * res$se, 1)
    )
  }, error = function(e) NULL)
  if (!is.null(ppi_est)) results$PPI <- ppi_est

  # 3. EIF (mu)
  eif_est <- tryCatch({
    res <- multiclass_ppi(T = pre$T, N = pre$N, K = K, g = "mu", alpha = 0.5, project = TRUE)
    tibble(
      method = "EIF", class = 1:K, pi_true = pi_true,
      pi_hat = res$pi_hat, se = res$se,
      ci_lower = pmax(res$pi_hat - z * res$se, 0),
      ci_upper = pmin(res$pi_hat + z * res$se, 1)
    )
  }, error = function(e) NULL)
  if (!is.null(eif_est)) results$EIF <- eif_est

  # 4. Rogan-Gladen (plug-in M)
  rg_est <- tryCatch({
    res <- multiclass_rg_plugin(T = pre$T, N = pre$N, K = K, alpha = 0.5, project = TRUE)
    tibble(
      method = "Rogan-Gladen", class = 1:K, pi_true = pi_true,
      pi_hat = res$pi_hat, se = res$se,
      ci_lower = pmax(res$pi_hat - z * res$se, 0),
      ci_upper = pmin(res$pi_hat + z * res$se, 1)
    )
  }, error = function(e) NULL)
  if (!is.null(rg_est)) results$`Rogan-Gladen` <- rg_est

  # 5. MLE (EM)
  mle_est <- tryCatch({
    res <- multiclass_mle_em(T = pre$T, N = pre$N, alpha = 0.5)
    # MLE doesn't return SE directly; leave as NA for now
    tibble(
      method = "MLE", class = 1:K, pi_true = pi_true,
      pi_hat = res$pi_hat, se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_
    )
  }, error = function(e) NULL)
  if (!is.null(mle_est)) results$MLE <- mle_est

  if (length(results) == 0) return(NULL)

  bind_rows(results) %>%
    mutate(
      pi1 = pi1,
      accuracy = accuracy,
      m_cal_ratio = m_cal_ratio,
      m_cal = m_cal,
      confusion_type = confusion_type,
      rep = rep_id
    )
}

# =============================================================================
# Run Simulation (like binary case structure)
# =============================================================================

run_one_grid <- function(pi1, accuracy, m_cal_ratio, confusion_type) {
  map_dfr(
    1:SIM_CONFIG$B,
    ~ run_one_sim(pi1, accuracy, m_cal_ratio, confusion_type, .x)
  )
}

message("Setting up parameter grid for K=", K, "...")

param_grid <- expand_grid(
  pi1 = SIM_CONFIG$pi1_vals,
  accuracy = SIM_CONFIG$accuracy_vals,
  m_cal_ratio = SIM_CONFIG$m_cal_ratios,
  confusion_type = SIM_CONFIG$confusion_types
)

message("Running ", nrow(param_grid), " parameter combinations x ", SIM_CONFIG$B, " reps...")

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- file.path("results", paste0("multinomial_K", K, "_", timestamp))
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

res_all <- future_pmap_dfr(
  list(
    pi1 = param_grid$pi1,
    accuracy = param_grid$accuracy,
    m_cal_ratio = param_grid$m_cal_ratio,
    confusion_type = param_grid$confusion_type
  ),
  run_one_grid,
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

# =============================================================================
# Summarize Results (following binary structure)
# =============================================================================

message("Summarizing results...")

# CI Summary (like binary ci_summary)
ci_summary <- res_all %>%
  filter(!is.na(pi_hat)) %>%
  mutate(
    covered = !is.na(ci_lower) & !is.na(ci_upper) &
              ci_lower <= pi_true & pi_true <= ci_upper
  ) %>%
  group_by(pi1, accuracy, m_cal_ratio, confusion_type, method, class) %>%
  summarise(
    pi_true = first(pi_true),
    mean_pi_hat = mean(pi_hat, na.rm = TRUE),
    bias = mean(pi_hat - pi_true, na.rm = TRUE),
    mse = mean((pi_hat - pi_true)^2, na.rm = TRUE),
    var_est = var(pi_hat, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    n_valid = sum(!is.na(pi_hat)),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_gap = coverage - SIM_CONFIG$nominal_coverage
  )

# Summary focused on class 1 (theta equivalent)
ci_summary_class1 <- ci_summary %>%
  filter(class == 1) %>%
  rename(theta_true = pi_true)

# Summary aggregated across all classes
ci_summary_aggregate <- ci_summary %>%
  group_by(pi1, accuracy, m_cal_ratio, confusion_type, method) %>%
  summarise(
    mean_abs_bias = mean(abs(bias), na.rm = TRUE),
    mean_mse = mean(mse, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(ci_summary, file.path(results_dir, "ci_summary.csv"))
write_csv(ci_summary_class1, file.path(results_dir, "ci_summary_class1.csv"))
write_csv(ci_summary_aggregate, file.path(results_dir, "ci_summary_aggregate.csv"))

# =============================================================================
# Plotting (similar to binary case)
# =============================================================================

message("Creating plots...")

# Save plot helper (wraps ggsave with plots_dir)
save_plot <- function(plot_obj, filename, width = 14, height = 12) {
  ggsave(file.path(plots_dir, filename), plot_obj, width = width, height = height, dpi = 400)
}

# Color palette (matching binary)
method_colors <- c(
  "Naive" = "#D55E00",
  "PPI" = "#E69F00",
  "EIF" = "#56B4E9",
  "Rogan-Gladen" = "#009E73",
  "MLE" = "#0072B2"
)

theme_sim <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.key.size = unit(1.2, "cm"),
      legend.key.width = unit(2, "cm")
    )
}

# Prepare data for plotting (class 1 focus)
ci_summary_plot <- ci_summary_class1 %>%
  mutate(
    accuracy_lab = paste0("Accuracy = ", accuracy),
    m_cal_ratio_lab = factor(format_cal_ratio(m_cal_ratio),
                              levels = format_cal_ratio(sort(unique(m_cal_ratio)))),
    bias_pct = dplyr::if_else(theta_true != 0, 100 * bias / theta_true, NA_real_)
  )

# -----------------------------------------------------------------------------
# GRID PLOTS: One plot per (confusion_type, accuracy), faceted by m_cal_ratio
# -----------------------------------------------------------------------------

for (conf_type in unique(ci_summary_plot$confusion_type)) {
  for (acc in unique(ci_summary_plot$accuracy)) {
    suffix <- paste0("_", conf_type, "_acc", as.integer(acc * 100), ".png")
    subtitle <- paste0("K=", K, " | ", conf_type, " confusion | accuracy=", acc)

    ci_sub <- filter(ci_summary_plot, confusion_type == conf_type, accuracy == acc)

    # Coverage
    p_cov <- ggplot(ci_sub, aes(x = pi1, y = coverage, color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      labs(
        title = "Coverage (Class 1)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "Coverage",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_cov, paste0("coverage", suffix), 16, 8)

    # CI Width
    p_width <- ggplot(ci_sub, aes(x = pi1, y = ci_width, color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      labs(
        title = "CI Width (Class 1)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "Mean CI Width",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_width, paste0("ci_width", suffix), 16, 8)

    # CI Width zoomed
    p_width_zoom <- p_width + coord_cartesian(ylim = c(0, 0.3)) +
      labs(title = "CI Width (Class 1, zoomed)")
    save_plot(p_width_zoom, paste0("ci_width_zoomed", suffix), 16, 8)

    # Bias
    p_bias <- ggplot(ci_sub, aes(x = pi1, y = bias, color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      labs(
        title = "Bias (Class 1)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "Bias",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_bias, paste0("bias", suffix), 16, 8)

    # Percent Bias
    p_bias_pct <- ggplot(ci_sub, aes(x = pi1, y = bias_pct, color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      labs(
        title = "Percent Bias (Class 1)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "Bias (%)",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_bias_pct, paste0("bias_pct", suffix), 16, 8)

    # RMSE (sqrt of mse)
    p_rmse <- ggplot(ci_sub, aes(x = pi1, y = sqrt(mse), color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      labs(
        title = "RMSE (Class 1)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "RMSE",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_rmse, paste0("rmse", suffix), 16, 8)

    # Plots excluding Naive
    ci_sub_no_naive <- filter(ci_sub, method != "Naive")

    p_cov_nn <- ggplot(ci_sub_no_naive, aes(x = pi1, y = coverage, color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      labs(
        title = "Coverage (Class 1, excl. Naive)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "Coverage",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_cov_nn, paste0("coverage_no_naive", suffix), 16, 8)

    p_width_nn <- ggplot(ci_sub_no_naive, aes(x = pi1, y = ci_width, color = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
      coord_cartesian(ylim = c(0, 0.3)) +
      labs(
        title = "CI Width (Class 1, excl. Naive, zoomed)",
        subtitle = subtitle,
        x = expression(pi[1]~"(true)"),
        y = "Mean CI Width",
        color = "Method"
      ) +
      theme_sim()
    save_plot(p_width_nn, paste0("ci_width_no_naive_zoomed", suffix), 16, 8)
  }
}

# -----------------------------------------------------------------------------
# AGGREGATE PLOTS: Across confusion types and accuracies
# -----------------------------------------------------------------------------

agg_plot_data <- ci_summary_aggregate %>%
  mutate(
    accuracy_lab = paste0("Accuracy = ", accuracy),
    m_cal_ratio_lab = factor(format_cal_ratio(m_cal_ratio),
                              levels = format_cal_ratio(sort(unique(m_cal_ratio))))
  )

for (conf_type in unique(agg_plot_data$confusion_type)) {
  agg_sub <- filter(agg_plot_data, confusion_type == conf_type)
  suffix <- paste0("_aggregate_", conf_type, ".png")

  p_cov_agg <- ggplot(agg_sub, aes(x = pi1, y = mean_coverage, color = method, linetype = accuracy_lab)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
    facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
    labs(
      title = paste0("Mean Coverage (All Classes) - K=", K, ", ", conf_type),
      x = expression(pi[1]~"(true)"),
      y = "Mean Coverage",
      color = "Method",
      linetype = "Accuracy"
    ) +
    theme_sim()
  save_plot(p_cov_agg, paste0("coverage", suffix), 18, 8)

  p_rmse_agg <- ggplot(agg_sub, aes(x = pi1, y = sqrt(mean_mse), color = method, linetype = accuracy_lab)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ m_cal_ratio_lab, nrow = 1) +
    scale_color_manual(values = method_colors) +
    scale_x_continuous(breaks = SIM_CONFIG$pi1_vals) +
    labs(
      title = paste0("Mean RMSE (All Classes) - K=", K, ", ", conf_type),
      x = expression(pi[1]~"(true)"),
      y = "Mean RMSE",
      color = "Method",
      linetype = "Accuracy"
    ) +
    theme_sim()
  save_plot(p_rmse_agg, paste0("rmse", suffix), 18, 8)
}

message("Results saved to: ", results_dir)

# Print summary
cat("\n========================================\n")
cat("SIMULATION SUMMARY - K =", K, "\n")
cat("========================================\n\n")

ci_summary_class1 %>%
  filter(confusion_type == "diagonal") %>%
  select(pi1, accuracy, m_cal_ratio, method, bias, coverage, ci_width) %>%
  arrange(pi1, accuracy, m_cal_ratio, method) %>%
  print(n = 50)
