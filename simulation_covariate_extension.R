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

safe_workers <- tryCatch({
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || cores <= 2) 1 else max(1, cores - 1)
}, error = function(e) 1)
if (safe_workers > 1) {
  plan(multisession, workers = safe_workers)
} else {
  plan(sequential)
}
set.seed(321)

SIM_COV <- list(
  sample_size = 4000L,
  label_ratios = c(0.05, 0.10, 0.20),
  B = 400L,
  covariate_types = c("continuous", "categorical"),
  nominal = 0.95
)

simulate_covariate_data <- function(N, covariate_type) {
  if (covariate_type == "continuous") {
    x_cont <- rnorm(N)
    pi_true <- plogis(-0.2 + 0.8 * x_cont)
    Y <- rbinom(N, 1, pi_true)
    q1 <- plogis(1 + 0.6 * x_cont)
    q0 <- plogis(1 - 0.6 * x_cont)
    prob_hat <- ifelse(Y == 1, q1, 1 - q0)
    S <- rbinom(N, 1, prob_hat)
    tibble(Y = Y, S = S, x_cont = x_cont, x_cat = factor("A"), theta_true = mean(pi_true))
  } else {
    levels_cat <- c("science", "math", "history")
    x_cat <- sample(levels_cat, size = N, replace = TRUE, prob = c(0.3, 0.4, 0.3))
    pi_map <- c(science = 0.25, math = 0.45, history = 0.65)
    Y <- rbinom(N, 1, pi_map[x_cat])
    q1_map <- c(science = 0.9, math = 0.8, history = 0.7)
    q0_map <- c(science = 0.95, math = 0.85, history = 0.75)
    q1 <- q1_map[x_cat]
    q0 <- q0_map[x_cat]
    prob_hat <- ifelse(Y == 1, q1, 1 - q0)
    S <- rbinom(N, 1, prob_hat)
    tibble(Y = Y, S = S, x_cont = NA_real_, x_cat = factor(x_cat), theta_true = sum(pi_map[x_cat]) / N)
  }
}

get_covariate_df <- function(df, covariate_type) {
  if (covariate_type == "continuous") {
    data.frame(x_cont = df$x_cont)
  } else {
    data.frame(x_cat = df$x_cat)
  }
}

run_one_condition <- function(covariate_type, label_ratio) {
  sample_size <- SIM_COV$sample_size
  n_cal <- max(50L, round(label_ratio * sample_size))
  future_map_dfr(
    1:SIM_COV$B,
    ~ {
      data_all <- simulate_covariate_data(sample_size + n_cal, covariate_type)
      test_idx <- seq_len(sample_size)
      cal_idx <- (sample_size + 1):(sample_size + n_cal)
      theta_true <- mean(data_all$Y)
      S_test <- data_all$S[test_idx]
      S_cal <- data_all$S[cal_idx]
      Y_cal <- data_all$Y[cal_idx]
      cov_test <- get_covariate_df(data_all[test_idx, ], covariate_type)
      cov_cal <- get_covariate_df(data_all[cal_idx, ], covariate_type)

      # Baseline identity PPI
      est_ppi <- ppi_point_and_ci(Y_cal, f_L = S_cal, f_U = S_test)

      # Covariate-adjusted RG style g
      mis_models <- fit_misclassification_models(cov_cal, Y_cal, S_cal)
      f_L_cov <- g_covariate_rg(S_cal, cov_cal, mis_models)
      f_U_cov <- g_covariate_rg(S_test, cov_test, mis_models)
      est_cov <- ppi_point_and_ci(Y_cal, f_L = f_L_cov, f_U = f_U_cov)

      # RePPI-style regression of Y on (S, X)
      reppi_mod <- fit_reppi_model(cov_cal, S_cal, Y_cal)
      f_L_reppi <- g_reppi_predict(S_cal, cov_cal, reppi_mod)
      f_U_reppi <- g_reppi_predict(S_test, cov_test, reppi_mod)
      est_reppi <- ppi_point_and_ci(Y_cal, f_L = f_L_reppi, f_U = f_U_reppi)

      tibble(
        method = c("ppi_identity", "ppi_covariate_rg", "ppi_reppi"),
        theta_hat = c(est_ppi$theta, est_cov$theta, est_reppi$theta),
        var_hat = c(est_ppi$var, est_cov$var, est_reppi$var),
        ci_lower = c(est_ppi$ci_lower, est_cov$ci_lower, est_reppi$ci_lower),
        ci_upper = c(est_ppi$ci_upper, est_cov$ci_upper, est_reppi$ci_upper)
      ) %>%
        mutate(
          theta_true = theta_true,
          label_ratio = label_ratio,
          covariate_type = covariate_type,
          covered = theta_true >= ci_lower & theta_true <= ci_upper,
          rep = .x
        )
    },
    .options = furrr_options(seed = TRUE)
  )
}

param_grid <- expand_grid(
  covariate_type = SIM_COV$covariate_types,
  label_ratio = SIM_COV$label_ratios
)

results_time <- format(Sys.time(), "%Y%m%d-%H%M%S")
results_dir <- file.path("results", paste0("covariate-", results_time))
plots_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, showWarnings = FALSE)

cov_results <- future_pmap_dfr(
  list(covariate_type = param_grid$covariate_type,
       label_ratio = param_grid$label_ratio),
  run_one_condition,
  .progress = TRUE
)

readr::write_csv(cov_results, file.path(results_dir, "raw_covariate_results.csv"))

summary_cov <- cov_results %>%
  group_by(covariate_type, label_ratio, method) %>%
  summarise(
    coverage = mean(covered),
    ci_width = mean(ci_upper - ci_lower),
    bias = mean(theta_hat - theta_true),
    mse = mean((theta_hat - theta_true)^2),
    .groups = "drop"
  ) %>%
  mutate(
    coverage_gap = coverage - SIM_COV$nominal,
    label_ratio_lab = factor(paste0(scales::percent(label_ratio, accuracy = 1), " labeled"),
                             levels = paste0(scales::percent(SIM_COV$label_ratios, accuracy = 1), " labeled"))
  )

readr::write_csv(summary_cov, file.path(results_dir, "summary_covariate.csv"))

plot_theme <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.background = element_rect(fill = "white"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 11)
    )
}

coverage_plot <- ggplot(summary_cov,
                        aes(x = covariate_type, y = coverage_gap, color = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  facet_wrap(~ label_ratio_lab) +
  labs(
    title = "Coverage - nominal across covariate-dependent judges",
    x = "Covariate type",
    y = "Coverage gap",
    color = "Estimator"
  ) +
  plot_theme()

ggsave(file.path(plots_dir, "covariate_coverage_gap.png"), coverage_plot, width = 10, height = 6, dpi = 320)

ci_plot <- ggplot(summary_cov,
                  aes(x = covariate_type, y = ci_width, color = method)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  facet_wrap(~ label_ratio_lab) +
  labs(
    title = "Average CI width across covariate types",
    x = "Covariate type",
    y = "CI width",
    color = "Estimator"
  ) +
  plot_theme()

ggsave(file.path(plots_dir, "covariate_ci_width.png"), ci_plot, width = 10, height = 6, dpi = 320)

bias_plot <- ggplot(summary_cov,
                    aes(x = covariate_type, y = bias, color = method)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(position = position_dodge(width = 0.4), size = 2) +
  facet_wrap(~ label_ratio_lab) +
  labs(
    title = "Bias across covariate types",
    x = "Covariate type",
    y = "Bias",
    color = "Estimator"
  ) +
  plot_theme()

ggsave(file.path(plots_dir, "covariate_bias.png"), bias_plot, width = 10, height = 6, dpi = 320)

message("Covariate-dependent simulation results saved to ", results_dir)
