library(tidyverse)
library(future)
library(furrr)

# Source the multinomial estimators
source("R/multinomial_estimators.R")

plan(multisession, workers = max(1, parallel::detectCores() - 1))
set.seed(123)

# =============================================================================
# Simulation Configuration
# =============================================================================

SIM_CONFIG <- list(
  K_vals = c(3, 5),                    # Number of classes
  n_test = 2000L,                       # Test set size
  m_cal_ratios = c(0.05, 0.10, 0.20),   # Calibration set size as fraction of test
  B = 500L,                             # Number of replications
  nominal_coverage = 0.95,
  # Different confusion matrix types
  confusion_types = c("diagonal", "noisy", "asymmetric")
)

# =============================================================================
# Data Generation Functions
# =============================================================================

#' Generate true prevalence vector
#'
#' @param K Number of classes
#' @param type "uniform", "skewed", or "sparse"
generate_pi_true <- function(K, type = "uniform") {
  if (type == "uniform") {
    rep(1/K, K)
  } else if (type == "skewed") {
    # First class dominant
    pi <- c(0.5, rep(0.5/(K-1), K-1))
    pi / sum(pi)
  } else if (type == "sparse") {
    # Only first two classes have mass
    pi <- c(0.6, 0.3, rep(0.1/(K-2), K-2))
    pi / sum(pi)
  } else {
    rep(1/K, K)
  }
}

#' Generate confusion matrix M where M[a,k] = P(Yhat=a | Y=k)
#'
#' @param K Number of classes
#' @param type "diagonal" (accurate), "noisy", or "asymmetric"
#' @param accuracy Base accuracy (diagonal elements)
generate_confusion_matrix <- function(K, type = "diagonal", accuracy = 0.7) {
  if (type == "diagonal") {
    # Diagonal-dominant: correct class with prob `accuracy`, uniform over others
    M <- matrix((1 - accuracy) / (K - 1), nrow = K, ncol = K)
    diag(M) <- accuracy
  } else if (type == "noisy") {
    # More noise, lower accuracy
    acc <- 0.5
    M <- matrix((1 - acc) / (K - 1), nrow = K, ncol = K)
    diag(M) <- acc
  } else if (type == "asymmetric") {
    # Adjacent classes get confused more
    M <- matrix(0.05 / (K - 1), nrow = K, ncol = K)
    for (k in 1:K) {
      M[k, k] <- accuracy
      if (k > 1) M[k-1, k] <- (1 - accuracy) * 0.6  # confuse with previous
      if (k < K) M[k+1, k] <- (1 - accuracy) * 0.3  # confuse with next
    }
    # Normalize columns
    M <- sweep(M, 2, colSums(M), "/")
  }

  # Ensure columns sum to 1
  M <- sweep(M, 2, colSums(M), "/")
  M
}

#' Generate simulation data
#'
#' @param n_test Test set size
#' @param m_cal Calibration set size
#' @param pi_true True prevalence vector
#' @param M_true True confusion matrix
#' @return List with test and calibration data
generate_data <- function(n_test, m_cal, pi_true, M_true) {
  K <- length(pi_true)
  N <- n_test + m_cal

  # Generate true labels Y from pi_true
  Y <- sample(1:K, N, replace = TRUE, prob = pi_true)

  # Generate predicted labels Yhat from M_true
  Yhat <- sapply(Y, function(y) {
    sample(1:K, 1, prob = M_true[, y])
  })

  # Split into test and calibration
  idx_cal <- sample(1:N, m_cal)
  idx_test <- setdiff(1:N, idx_cal)

  list(
    yhat_test = Yhat[idx_test],
    y_cal = Y[idx_cal],
    yhat_cal = Yhat[idx_cal],
    y_test = Y[idx_test]  # for evaluation only
  )
}

# =============================================================================
# Single Simulation Replicate
# =============================================================================

run_one_sim <- function(K, n_test, m_cal_ratio, confusion_type, pi_type, rep_id) {
  m_cal <- round(n_test * m_cal_ratio)
  pi_true <- generate_pi_true(K, pi_type)
  M_true <- generate_confusion_matrix(K, confusion_type)

  # Generate data
  data <- generate_data(n_test, m_cal, pi_true, M_true)

  # Run all estimators
  results <- list()

  # 1. Naive (no correction)
  naive_est <- tryCatch({
    naive_multinomial(data$yhat_test, K)
  }, error = function(e) NULL)

  if (!is.null(naive_est)) {
    results$Naive <- tibble(
      method = "Naive",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = naive_est$pi_hat,
      se = naive_est$se,
      ci_lower = naive_est$ci_lower,
      ci_upper = naive_est$ci_upper
    )
  }

  # 2. PPI (Model A, identity surrogate)
  ppi_est <- tryCatch({
    ppi_multinomial(data$yhat_test, data$y_cal, data$yhat_cal, K)
  }, error = function(e) NULL)

  if (!is.null(ppi_est)) {
    results$PPI <- tibble(
      method = "PPI",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = ppi_est$pi_hat,
      se = ppi_est$se,
      ci_lower = ppi_est$ci_lower,
      ci_upper = ppi_est$ci_upper
    )
  }

  # 3. EIF (Model A, with mu estimation)
  eif_est <- tryCatch({
    eif_multinomial(data$yhat_test, data$y_cal, data$yhat_cal, K, smooth = 0.5)
  }, error = function(e) NULL)

  if (!is.null(eif_est)) {
    results$EIF <- tibble(
      method = "EIF",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = eif_est$pi_hat,
      se = eif_est$se,
      ci_lower = eif_est$ci_lower,
      ci_upper = eif_est$ci_upper
    )
  }

  # 4. Rogan-Gladen (matrix inversion)
  rg_est <- tryCatch({
    rg_multinomial(data$yhat_test, data$y_cal, data$yhat_cal, K, smooth = 0.5)
  }, error = function(e) NULL)

  if (!is.null(rg_est)) {
    results$`Rogan-Gladen` <- tibble(
      method = "Rogan-Gladen",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = rg_est$pi_hat,
      se = rg_est$se,
      ci_lower = rg_est$ci_lower,
      ci_upper = rg_est$ci_upper
    )
  }

  # 5. Ridge-stabilized RG
  ridge_est <- tryCatch({
    rg_ridge_multinomial(data$yhat_test, data$y_cal, data$yhat_cal, K, gamma = 0.01, smooth = 0.5)
  }, error = function(e) NULL)

  if (!is.null(ridge_est)) {
    results$`RG-Ridge` <- tibble(
      method = "RG-Ridge",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = ridge_est$pi_hat,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  }

  # 6. Constrained Least Squares
  cls_est <- tryCatch({
    cls_multinomial(data$yhat_test, data$y_cal, data$yhat_cal, K, smooth = 0.5)
  }, error = function(e) NULL)

  if (!is.null(cls_est)) {
    results$CLS <- tibble(
      method = "CLS",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = cls_est$pi_hat,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  }

  # 7. MLE (Model B)
  mle_est <- tryCatch({
    mle_multinomial(data$y_cal, data$yhat_cal, data$yhat_test, K, smooth = 0.5)
  }, error = function(e) NULL)

  if (!is.null(mle_est)) {
    results$MLE <- tibble(
      method = "MLE",
      class = 1:K,
      pi_true = pi_true,
      pi_hat = mle_est$pi_hat,
      se = mle_est$se,
      ci_lower = mle_est$ci_lower,
      ci_upper = mle_est$ci_upper
    )
  }

  # Combine results
  if (length(results) == 0) return(NULL)

  bind_rows(results) %>%
    mutate(
      K = K,
      n_test = n_test,
      m_cal_ratio = m_cal_ratio,
      m_cal = m_cal,
      confusion_type = confusion_type,
      pi_type = pi_type,
      rep = rep_id
    )
}

# =============================================================================
# Run Simulation
# =============================================================================

message("Setting up parameter grid...")

param_grid <- expand_grid(
  K = SIM_CONFIG$K_vals,
  m_cal_ratio = SIM_CONFIG$m_cal_ratios,
  confusion_type = SIM_CONFIG$confusion_types,
  pi_type = c("uniform", "skewed")
)

message("Running ", nrow(param_grid), " parameter combinations x ", SIM_CONFIG$B, " reps...")

# Run simulations
results_list <- future_pmap(

list(
K = rep(param_grid$K, each = SIM_CONFIG$B),
m_cal_ratio = rep(param_grid$m_cal_ratio, each = SIM_CONFIG$B),
confusion_type = rep(param_grid$confusion_type, each = SIM_CONFIG$B),
pi_type = rep(param_grid$pi_type, each = SIM_CONFIG$B),
rep_id = rep(1:SIM_CONFIG$B, nrow(param_grid))
),
function(K, m_cal_ratio, confusion_type, pi_type, rep_id) {
run_one_sim(K, SIM_CONFIG$n_test, m_cal_ratio, confusion_type, pi_type, rep_id)
},
.progress = TRUE,
.options = furrr_options(seed = TRUE)
)

res_all <- bind_rows(results_list)

# =============================================================================
# Summarize Results
# =============================================================================

message("Summarizing results...")

# Create output directory
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- file.path("results", paste0("multinomial_", timestamp))
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

# Summary statistics
summary_stats <- res_all %>%
  filter(!is.na(pi_hat)) %>%
  mutate(
    bias = pi_hat - pi_true,
    sq_error = (pi_hat - pi_true)^2,
    covered = !is.na(ci_lower) & !is.na(ci_upper) &
              ci_lower <= pi_true & pi_true <= ci_upper,
    ci_width = ci_upper - ci_lower
  ) %>%
  group_by(K, m_cal_ratio, confusion_type, pi_type, method, class) %>%
  summarise(
    pi_true = first(pi_true),
    mean_pi_hat = mean(pi_hat, na.rm = TRUE),
    bias = mean(bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    n_valid = sum(!is.na(pi_hat)),
    .groups = "drop"
  )

# Aggregate across classes
summary_aggregate <- summary_stats %>%
  group_by(K, m_cal_ratio, confusion_type, pi_type, method) %>%
  summarise(
    mean_abs_bias = mean(abs(bias), na.rm = TRUE),
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    mean_ci_width = mean(mean_ci_width, na.rm = TRUE),
    .groups = "drop"
  )

# Save results
write_csv(res_all, file.path(results_dir, "raw_results.csv"))
write_csv(summary_stats, file.path(results_dir, "summary_by_class.csv"))
write_csv(summary_aggregate, file.path(results_dir, "summary_aggregate.csv"))

# =============================================================================
# Plotting
# =============================================================================

message("Creating plots...")

# Color palette for methods
method_colors <- c(
  "Naive" = "#D55E00",
  "PPI" = "#E69F00",
  "EIF" = "#56B4E9",
  "Rogan-Gladen" = "#009E73",
  "RG-Ridge" = "#0072B2",
  "CLS" = "#CC79A7",
  "MLE" = "#F0E442"
)

# Plot theme
theme_sim <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom"
    )
}

# 1. Bias plot by method and K
p_bias <- summary_aggregate %>%
  mutate(m_cal_pct = paste0(100 * m_cal_ratio, "% cal")) %>%
  ggplot(aes(x = method, y = mean_abs_bias, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(K ~ m_cal_pct + confusion_type, scales = "free_y",
             labeller = labeller(K = function(x) paste0("K = ", x))) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Mean Absolute Bias by Method",
    x = "Method",
    y = "Mean |Bias|",
    fill = "Method"
  ) +
  theme_sim() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "bias_by_method.png"), p_bias, width = 14, height = 8, dpi = 300)

# 2. RMSE plot
p_rmse <- summary_aggregate %>%
  mutate(m_cal_pct = paste0(100 * m_cal_ratio, "% cal")) %>%
  ggplot(aes(x = method, y = mean_rmse, fill = method)) +
  geom_col(position = "dodge") +
  facet_grid(K ~ m_cal_pct + confusion_type, scales = "free_y",
             labeller = labeller(K = function(x) paste0("K = ", x))) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Root Mean Squared Error by Method",
    x = "Method",
    y = "RMSE",
    fill = "Method"
  ) +
  theme_sim() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "rmse_by_method.png"), p_rmse, width = 14, height = 8, dpi = 300)

# 3. Coverage plot (only for methods with CIs)
p_coverage <- summary_aggregate %>%
  filter(!is.na(mean_coverage)) %>%
  mutate(m_cal_pct = paste0(100 * m_cal_ratio, "% cal")) %>%
  ggplot(aes(x = method, y = mean_coverage, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_grid(K ~ m_cal_pct + confusion_type,
             labeller = labeller(K = function(x) paste0("K = ", x))) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Coverage of 95% Confidence Intervals",
    x = "Method",
    y = "Coverage",
    fill = "Method"
  ) +
  theme_sim() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plots_dir, "coverage_by_method.png"), p_coverage, width = 14, height = 8, dpi = 300)

# 4. Line plot: RMSE vs calibration ratio
p_rmse_line <- summary_aggregate %>%
  ggplot(aes(x = m_cal_ratio, y = mean_rmse, color = method, linetype = confusion_type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(K ~ pi_type, labeller = labeller(K = function(x) paste0("K = ", x))) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(labels = scales::percent) +
  labs(
    title = "RMSE vs Calibration Set Size",
    x = "Calibration Ratio (% of test set)",
    y = "Mean RMSE",
    color = "Method",
    linetype = "Confusion Type"
  ) +
  theme_sim()

ggsave(file.path(plots_dir, "rmse_vs_cal_ratio.png"), p_rmse_line, width = 12, height = 8, dpi = 300)

# 5. Bias by class (for K=5, uniform pi)
p_bias_class <- summary_stats %>%
  filter(K == 5, pi_type == "uniform", m_cal_ratio == 0.1) %>%
  ggplot(aes(x = factor(class), y = bias, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ confusion_type) +
  scale_fill_manual(values = method_colors) +
  labs(
    title = "Bias by Class (K=5, 10% calibration, uniform pi)",
    x = "Class",
    y = "Bias",
    fill = "Method"
  ) +
  theme_sim()

ggsave(file.path(plots_dir, "bias_by_class_K5.png"), p_bias_class, width = 12, height = 6, dpi = 300)

# 6. Coverage by class
p_coverage_class <- summary_stats %>%
  filter(K == 5, pi_type == "uniform", m_cal_ratio == 0.1, !is.na(coverage)) %>%
  ggplot(aes(x = factor(class), y = coverage, fill = method)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~ confusion_type) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Coverage by Class (K=5, 10% calibration, uniform pi)",
    x = "Class",
    y = "Coverage",
    fill = "Method"
  ) +
  theme_sim()

ggsave(file.path(plots_dir, "coverage_by_class_K5.png"), p_coverage_class, width = 12, height = 6, dpi = 300)

message("Results saved to: ", results_dir)

# Print summary table
cat("\n========================================\n")
cat("SIMULATION SUMMARY (aggregated)\n")
cat("========================================\n\n")

summary_aggregate %>%
  filter(pi_type == "uniform") %>%
  arrange(K, m_cal_ratio, confusion_type, method) %>%
  print(n = 100)
