library(tidyverse)
library(future)
library(furrr)

if (requireNamespace("debiasLLMReporting", quietly = TRUE)) {
  library(debiasLLMReporting)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  stop("Install the debiasLLMReporting package or devtools to load helper functions.")
}

plan(multisession, workers = max(1, parallel::detectCores() - 1))
set.seed(202)

SIM_MULTI <- list(
  sample_size = 5000L,
  label_ratios = c(0.01, 0.05, 0.10, 0.20, 0.50),
  B = 500L,
  nominal_coverage = 0.95
)

# Prevalence patterns: balanced (approx uniform), moderate skew, and extreme skew
# for both K=3 and K=5 settings.
PI_PATTERNS <- tribble(
  ~pi_pattern,       ~K, ~pi_vec,
  "K3-balanced",     3L, list(c(1/3, 1/3, 1/3)),
  "K3-moderate",     3L, list(c(0.2, 0.3, 0.5)),
  "K3-extreme",      3L, list(c(0.05, 0.15, 0.80)),
  "K5-balanced",     5L, list(rep(0.2, 5)),
  "K5-moderate",     5L, list(c(0.4, 0.2, 0.15, 0.15, 0.10)),
  "K5-extreme",      5L, list(c(0.6, 0.15, 0.1, 0.1, 0.05))
)

# Accuracy profiles specify diagonal entries (per-class accuracies) for each K.
ACCURACY_PROFILES <- tribble(
  ~accuracy_id, ~K, ~diag_probs,
  "K3-high",    3L, list(c(0.92, 0.90, 0.94)),
  "K3-medium",  3L, list(c(0.78, 0.74, 0.80)),
  "K3-low",     3L, list(c(0.62, 0.58, 0.60)),
  "K5-high",    5L, list(c(0.90, 0.88, 0.91, 0.87, 0.92)),
  "K5-medium",  5L, list(c(0.78, 0.74, 0.75, 0.73, 0.76)),
  "K5-low",     5L, list(c(0.62, 0.60, 0.58, 0.56, 0.57))
)

format_label_ratio <- function(x) paste0(scales::percent(x, accuracy = 1), " labeled")

build_confusion_matrix <- function(diag_probs) {
  K <- length(diag_probs)
  M <- matrix(0, nrow = K, ncol = K)
  for (k in seq_len(K)) {
    remain <- max(1e-6, 1 - diag_probs[k])
    off_val <- remain / (K - 1)
    M[, k] <- off_val
    M[k, k] <- diag_probs[k]
  }
  M
}

draw_multiclass <- function(pi_vec, n) {
  sample(seq_along(pi_vec), size = n, replace = TRUE, prob = pi_vec)
}

derive_calibration_n <- function(label_ratio, sample_size, K) {
  max(K * 20L, round(label_ratio * sample_size))
}

safe_inverse <- function(M, jitter = 1e-4) {
  attempt <- try(solve(M), silent = TRUE)
  if (inherits(attempt, "try-error")) {
    solve(M + diag(jitter, nrow(M)))
  } else {
    attempt
  }
}

sim_one_condition_multi <- function(pi_vec,
                                    diag_probs,
                                    label_ratio,
                                    sample_size,
                                    nominal,
                                    K,
                                    pi_pattern,
                                    accuracy_id,
                                    rep_id) {
  n_cal <- derive_calibration_n(label_ratio, sample_size, K)
  N_total <- sample_size + n_cal
  Y <- draw_multiclass(pi_vec, N_total)
  conf_true <- build_confusion_matrix(diag_probs)
  S <- vapply(Y, function(y) sample.int(K, 1, prob = conf_true[, y]), integer(1))
  S_test <- S[seq_len(sample_size)]
  S_cal <- S[(sample_size + 1):N_total]
  Y_cal <- Y[(sample_size + 1):N_total]
  p_hat <- tabulate(S_test, nbins = K) / sample_size
  M_hat <- estimate_confusion_matrix(S_cal, Y_cal, K)
  W_hat <- safe_inverse(M_hat)
  context_id <- list(K = K)
  alpha <- 1 - nominal
  ppi_id <- ppi_multiclass(
    S_unlabeled = S_test,
    S_labeled = S_cal,
    Y_labeled = Y_cal,
    g_transform = g_identity_multiclass,
    context = context_id,
    alpha = alpha
  )
  ppi_pp <- ppi_pp_multiclass(
    S_unlabeled = S_test,
    S_labeled = S_cal,
    Y_labeled = Y_cal,
    K = K,
    alpha = alpha
  )
  cov_p <- (diag(p_hat) - tcrossprod(p_hat)) / sample_size
  cov_pi <- W_hat %*% cov_p %*% t(W_hat)
  var_rg <- pmax(diag(cov_pi), 0)
  z <- qnorm(1 - alpha / 2)
  pi_hat_rg_raw <- as.numeric(W_hat %*% p_hat)
  pi_hat_rg_proj <- project_simplex(pi_hat_rg_raw)
  ci_rg_lower <- pmax(pi_hat_rg_proj - z * sqrt(var_rg), 0)
  ci_rg_upper <- pmin(pi_hat_rg_proj + z * sqrt(var_rg), 1)
  pi_hat_rg_cstr <- rg_least_squares_simplex(M_hat, p_hat)
  ci_cstr_lower <- pmax(pi_hat_rg_cstr - z * sqrt(var_rg), 0)
  ci_cstr_upper <- pmin(pi_hat_rg_cstr + z * sqrt(var_rg), 1)
  pi_true_vec <- pi_vec
  assemble_results <- function(method_label, est_obj) {
    tibble(
      method = method_label,
      class_id = seq_len(K),
      pi_hat = est_obj$pi_hat,
      var_hat = est_obj$var,
      ci_lower = est_obj$ci_lower,
      ci_upper = est_obj$ci_upper
    )
  }
  res <- bind_rows(
    assemble_results("ppi_identity", ppi_id),
    assemble_results("ppi_pp", ppi_pp),
    tibble(
      method = "rg_projection",
      class_id = seq_len(K),
      pi_hat = pi_hat_rg_proj,
      var_hat = var_rg,
      ci_lower = ci_rg_lower,
      ci_upper = ci_rg_upper
    ),
    tibble(
      method = "rg_constrained",
      class_id = seq_len(K),
      pi_hat = pi_hat_rg_cstr,
      var_hat = var_rg,
      ci_lower = ci_cstr_lower,
      ci_upper = ci_cstr_upper
    )
  ) %>%
    mutate(
      pi_true = pi_true_vec[class_id],
      pi_pattern = pi_pattern,
      accuracy_id = accuracy_id,
      label_ratio = label_ratio,
      K = K,
      rep = rep_id,
      covered = pi_true >= ci_lower & pi_true <= ci_upper
    )
  res
}

expand_conditions <- function() {
  PI_PATTERNS %>%
    inner_join(ACCURACY_PROFILES,
               by = "K",
               suffix = c("_pi", "_acc")) %>%
    transmute(
      pi_pattern = pi_pattern,
      pi_vec = pi_vec,
      accuracy_id = accuracy_id,
      diag_probs = diag_probs,
      K = K
    )
}

param_grid <- expand_conditions() %>%
  crossing(tibble(label_ratio = SIM_MULTI$label_ratios))

run_param <- function(pi_vec, diag_probs, pi_pattern, accuracy_id, label_ratio, K) {
  pi_vec_num <- as.numeric(pi_vec[[1]])
  diag_vec <- as.numeric(diag_probs[[1]])
  if (length(pi_vec_num) != length(diag_vec)) {
    stop("pi_vec and diag_probs must have the same length (K)")
  }
  K_local <- length(pi_vec_num)
  if (!missing(K) && !is.null(K) && K != K_local) {
    warning("Provided K does not match vector length; using inferred value.")
  }
  K <- K_local
  future_map_dfr(
    1:SIM_MULTI$B,
    ~ sim_one_condition_multi(
      pi_vec = pi_vec_num,
      diag_probs = diag_vec,
      label_ratio = label_ratio,
      sample_size = SIM_MULTI$sample_size,
      nominal = SIM_MULTI$nominal_coverage,
      K = K,
      pi_pattern = pi_pattern,
      accuracy_id = accuracy_id,
      rep_id = .x
    ),
    .options = furrr_options(seed = TRUE)
  ) %>%
    mutate(
      pi_pattern = pi_pattern,
      accuracy_id = accuracy_id,
      K = K
    )
}

results_time <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- file.path("results", paste0("multiclass-", results_time))
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

res_multi <- future_pmap_dfr(
  list(
    pi_vec = param_grid$pi_vec,
    diag_probs = param_grid$diag_probs,
    pi_pattern = param_grid$pi_pattern,
    accuracy_id = param_grid$accuracy_id,
    label_ratio = param_grid$label_ratio,
    K = param_grid$K
  ),
  run_param,
  .progress = TRUE
)

readr::write_csv(res_multi, file.path(results_dir, "raw_multiclass_results.csv"))

label_ratio_levels <- format_label_ratio(sort(unique(SIM_MULTI$label_ratios)))

summary_multi <- res_multi %>%
  group_by(pi_pattern, accuracy_id, label_ratio, method, class_id, pi_true, K) %>%
  summarise(
    coverage = mean(covered),
    ci_width = mean(ci_upper - ci_lower),
    bias = mean(pi_hat - pi_true),
    mse = mean((pi_hat - pi_true)^2),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_gap = coverage - SIM_MULTI$nominal_coverage,
    label_ratio_lab = factor(format_label_ratio(label_ratio), levels = label_ratio_levels),
    class_label = paste0("Class ", class_id)
  )

readr::write_csv(summary_multi, file.path(results_dir, "summary_multiclass.csv"))

plot_theme <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "white"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 12)
    )
}

plot_metric_by_class <- function(metric_col, y_label, title_prefix, filename_prefix) {
  for (K_val in sort(unique(summary_multi$K))) {
    df_K <- summary_multi %>%
      filter(K == K_val)
    for (class_val in sort(unique(df_K$class_id))) {
      df_plot <- df_K %>%
        filter(class_id == class_val)
      if (nrow(df_plot) == 0) next
      plot_obj <- ggplot(df_plot,
                         aes_string(x = "pi_true",
                                    y = metric_col,
                                    color = "method")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey85") +
        geom_point(size = 2) +
        facet_grid(label_ratio_lab ~ accuracy_id) +
        labs(
          title = sprintf("%s (K=%d, Class %d)", title_prefix, K_val, class_val),
          x = "True prevalence (pi_k)",
          y = y_label,
          color = "Method"
        ) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        plot_theme()
      filename <- sprintf("%s_K%d_class%d.png", filename_prefix, K_val, class_val)
      ggsave(file.path(plots_dir, filename), plot_obj, width = 12, height = 9, dpi = 320)
    }
  }
}

plot_metric_by_class("coverage_gap", "Coverage - nominal", "Coverage Gap", "multiclass_coverage_gap")
plot_metric_by_class("ci_width", "CI width", "Confidence Interval Width", "multiclass_ci_width")
plot_metric_by_class("bias", "Bias", "Bias", "multiclass_bias")

message("Multiclass simulation results saved to ", results_dir)
