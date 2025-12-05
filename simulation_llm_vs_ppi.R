library(tidyverse)
library(future)
library(furrr)
library(debiasLLMReporting)

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
    split_label = c("q0_20", "q0_50", "q0_80"),
    m0_prop = c(0.2, 0.5, 0.8)
  )
)

format_label_ratio <- function(x) paste0(scales::percent(x, accuracy = 1), " labeled")

g_specs <- list(
  list(label = "ppi", transform = function(scores, context) scores)
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
    method = "ppi_pp",
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
    method = "llm",
    theta_hat = llm_est$theta,
    var_hat = llm_est$var,
    ci_lower = llm_est$ci_lower,
    ci_upper = llm_est$ci_upper,
    lambda = NA_real_
  )
  bind_rows(ppi_results, ppi_pp_tbl, llm_tbl) %>%
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

save_plot <- function(plot_obj, filename, width = 14, height = 12) {
  ggsave(
    filename = file.path(plots_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 320
  )
}

plot_base_theme <- function() {
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "white"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

ci_summary_plot <- ci_summary %>%
  mutate(
    method = recode(method,
                    "llm" = "measurement_error",
                    "ppi" = "PPI",
                    "ppi_pp" = "PPI++"),
    bias_pct = dplyr::if_else(theta_true != 0,
                              100 * bias / theta_true,
                              NA_real_),
    q0_lab = paste0("q0 = ", q0),
    q1_lab = paste0("q1 = ", q1)
  )

q_bias_long <- q_est_summary %>%
  pivot_longer(
    cols = c(q0_bias, q1_bias),
    names_to = "metric",
    values_to = "bias"
  ) %>%
  mutate(
    metric = recode(metric, q0_bias = "q0", q1_bias = "q1"),
    q0_lab = paste0("q0 = ", q0),
    q1_lab = paste0("q1 = ", q1)
  )

for (split in unique(ci_summary_plot$split_label)) {
  for (lr in unique(ci_summary_plot$label_ratio)) {

    lr_pct <- scales::percent(lr, accuracy = 1)
    lr_str <- gsub("%", "", lr_pct)
    filename_suffix <- paste0("_", split, "_labeled_", sprintf("%02d", as.integer(lr * 100)), ".png")
    subtitle_text <- paste0("m0:m1 = ", split, " | ", lr_pct, " labeled")

    df_sub <- ci_summary_plot %>%
      filter(split_label == split, label_ratio == lr)

    # Coverage gap plot
    p_cov <- ggplot(df_sub,
                    aes(x = theta_true, y = coverage_gap, color = method, group = method)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_grid(q0_lab ~ q1_lab) +
      labs(
        title = "Coverage - Nominal",
        subtitle = subtitle_text,
        x = expression(theta[true]),
        y = "Coverage gap",
        color = "Method"
      ) +
      scale_x_continuous(breaks = SIM_CONFIG$thetas) +
      plot_base_theme()
    save_plot(p_cov, paste0("coverage_gap", filename_suffix))

    # CI width plot
    p_width <- ggplot(df_sub,
                      aes(x = theta_true, y = ci_width, color = method, group = method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_grid(q0_lab ~ q1_lab) +
      labs(
        title = "CI Width Across Theta",
        subtitle = subtitle_text,
        x = expression(theta[true]),
        y = "Mean CI width",
        color = "Method"
      ) +
      scale_x_continuous(breaks = SIM_CONFIG$thetas) +
      plot_base_theme()
    save_plot(p_width, paste0("ci_width", filename_suffix))

    # Bias plot
    p_bias <- ggplot(df_sub,
                     aes(x = theta_true, y = bias, color = method, group = method)) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_grid(q0_lab ~ q1_lab) +
      labs(
        title = "Estimator Bias vs Theta",
        subtitle = subtitle_text,
        x = expression(theta[true]),
        y = "Bias",
        color = "Method"
      ) +
      scale_x_continuous(breaks = SIM_CONFIG$thetas) +
      plot_base_theme()
    save_plot(p_bias, paste0("bias", filename_suffix))

    # Bias percent plot
    p_bias_pct <- ggplot(df_sub,
                         aes(x = theta_true, y = bias_pct, color = method, group = method)) +
      geom_hline(yintercept = 0, linetype = "dotted") +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_grid(q0_lab ~ q1_lab) +
      labs(
        title = "Percent Bias vs Theta",
        subtitle = subtitle_text,
        x = expression(theta[true]),
        y = "Bias (%)",
        color = "Method"
      ) +
      scale_x_continuous(breaks = SIM_CONFIG$thetas) +
      plot_base_theme()
    save_plot(p_bias_pct, paste0("bias_pct", filename_suffix))

    # Q bias plot
    df_q_sub <- q_bias_long %>%
      filter(split_label == split, label_ratio == lr)

    p_q_bias <- ggplot(df_q_sub,
                       aes(x = theta_true, y = bias, color = metric, group = metric)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
      geom_line(linewidth = 1) +
      geom_point(size = 1.5) +
      facet_grid(q0_lab ~ q1_lab) +
      labs(
        title = "Calibration Estimate Bias (q0/q1)",
        subtitle = subtitle_text,
        x = expression(theta[true]),
        y = "Bias",
        color = "Parameter"
      ) +
      scale_x_continuous(breaks = SIM_CONFIG$thetas) +
      plot_base_theme()
    save_plot(p_q_bias, paste0("q_bias", filename_suffix))
  }
}

# =============================================================================
# ADDITIONAL PLOTS: q0 == q1 case (diagonal), faceted by q0/q1 and label_ratio
# =============================================================================

ci_summary_diag <- ci_summary_plot %>%
  filter(q0 == q1) %>%
  mutate(
    q_val_lab = paste0("q0 = q1 = ", q0),
    label_ratio_lab = factor(
      format_label_ratio(label_ratio),
      levels = format_label_ratio(sort(unique(label_ratio)))
    )
  )

q_bias_diag <- q_bias_long %>%
  filter(q0 == q1) %>%
  mutate(
    q_val_lab = paste0("q0 = q1 = ", q0),
    label_ratio_lab = factor(
      format_label_ratio(label_ratio),
      levels = format_label_ratio(sort(unique(label_ratio)))
    )
  )

for (split in unique(ci_summary_diag$split_label)) {

  filename_suffix <- paste0("_diag_", split, ".png")
  subtitle_text <- paste0("m0:m1 = ", split, " | q0 = q1 (diagonal)")

  df_sub <- ci_summary_diag %>%
    filter(split_label == split)

  # Coverage gap plot (diagonal)
  p_cov_diag <- ggplot(df_sub,
                       aes(x = theta_true, y = coverage_gap, color = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(q_val_lab ~ label_ratio_lab) +
    labs(
      title = "Coverage - Nominal (q0 = q1)",
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Coverage gap",
      color = "Method"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_cov_diag, paste0("coverage_gap", filename_suffix), width = 16, height = 10)

  # CI width plot (diagonal)
  p_width_diag <- ggplot(df_sub,
                         aes(x = theta_true, y = ci_width, color = method, group = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(q_val_lab ~ label_ratio_lab) +
    labs(
      title = "CI Width (q0 = q1)",
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Mean CI width",
      color = "Method"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_width_diag, paste0("ci_width", filename_suffix), width = 16, height = 10)

  # Bias plot (diagonal)
  p_bias_diag <- ggplot(df_sub,
                        aes(x = theta_true, y = bias, color = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(q_val_lab ~ label_ratio_lab) +
    labs(
      title = "Estimator Bias (q0 = q1)",
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Bias",
      color = "Method"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_bias_diag, paste0("bias", filename_suffix), width = 16, height = 10)

  # Bias percent plot (diagonal)
  p_bias_pct_diag <- ggplot(df_sub,
                            aes(x = theta_true, y = bias_pct, color = method, group = method)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(q_val_lab ~ label_ratio_lab) +
    labs(
      title = "Percent Bias (q0 = q1)",
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Bias (%)",
      color = "Method"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_bias_pct_diag, paste0("bias_pct", filename_suffix), width = 16, height = 10)

  # Q bias plot (diagonal)
  df_q_sub <- q_bias_diag %>%
    filter(split_label == split)

  p_q_bias_diag <- ggplot(df_q_sub,
                          aes(x = theta_true, y = bias, color = metric, group = metric)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_grid(q_val_lab ~ label_ratio_lab) +
    labs(
      title = "Calibration Estimate Bias (q0 = q1)",
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Bias",
      color = "Parameter"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_q_bias_diag, paste0("q_bias", filename_suffix), width = 16, height = 10)
}

# =============================================================================
# AGGREGATE DIAGONAL PLOTS: q0 == q1, all m0:m1 splits with linetype
# =============================================================================

# Coverage gap plot (diagonal, aggregate)
p_cov_diag_agg <- ggplot(ci_summary_diag,
                         aes(x = theta_true, y = coverage_gap,
                             color = method, linetype = split_label,
                             group = interaction(method, split_label))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(q_val_lab ~ label_ratio_lab) +
  labs(
    title = "Coverage - Nominal (q0 = q1)",
    subtitle = "All m0:m1 ratios",
    x = expression(theta[true]),
    y = "Coverage gap",
    color = "Method",
    linetype = "m0:m1"
  ) +
  scale_x_continuous(breaks = SIM_CONFIG$thetas) +
  plot_base_theme()
save_plot(p_cov_diag_agg, "coverage_gap_diag_aggregate.png", width = 16, height = 10)

# CI width plot (diagonal, aggregate)
p_width_diag_agg <- ggplot(ci_summary_diag,
                           aes(x = theta_true, y = ci_width,
                               color = method, linetype = split_label,
                               group = interaction(method, split_label))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(q_val_lab ~ label_ratio_lab) +
  labs(
    title = "CI Width (q0 = q1)",
    subtitle = "All m0:m1 ratios",
    x = expression(theta[true]),
    y = "Mean CI width",
    color = "Method",
    linetype = "m0:m1"
  ) +
  scale_x_continuous(breaks = SIM_CONFIG$thetas) +
  plot_base_theme()
save_plot(p_width_diag_agg, "ci_width_diag_aggregate.png", width = 16, height = 10)

# Bias plot (diagonal, aggregate)
p_bias_diag_agg <- ggplot(ci_summary_diag,
                          aes(x = theta_true, y = bias,
                              color = method, linetype = split_label,
                              group = interaction(method, split_label))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(q_val_lab ~ label_ratio_lab) +
  labs(
    title = "Estimator Bias (q0 = q1)",
    subtitle = "All m0:m1 ratios",
    x = expression(theta[true]),
    y = "Bias",
    color = "Method",
    linetype = "m0:m1"
  ) +
  scale_x_continuous(breaks = SIM_CONFIG$thetas) +
  plot_base_theme()
save_plot(p_bias_diag_agg, "bias_diag_aggregate.png", width = 16, height = 10)

# Bias percent plot (diagonal, aggregate)
p_bias_pct_diag_agg <- ggplot(ci_summary_diag,
                              aes(x = theta_true, y = bias_pct,
                                  color = method, linetype = split_label,
                                  group = interaction(method, split_label))) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(q_val_lab ~ label_ratio_lab) +
  labs(
    title = "Percent Bias (q0 = q1)",
    subtitle = "All m0:m1 ratios",
    x = expression(theta[true]),
    y = "Bias (%)",
    color = "Method",
    linetype = "m0:m1"
  ) +
  scale_x_continuous(breaks = SIM_CONFIG$thetas) +
  plot_base_theme()
save_plot(p_bias_pct_diag_agg, "bias_pct_diag_aggregate.png", width = 16, height = 10)

# Q bias plot (diagonal, aggregate)
p_q_bias_diag_agg <- ggplot(q_bias_diag,
                            aes(x = theta_true, y = bias,
                                color = metric, linetype = split_label,
                                group = interaction(metric, split_label))) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(q_val_lab ~ label_ratio_lab) +
  labs(
    title = "Calibration Estimate Bias (q0 = q1)",
    subtitle = "All m0:m1 ratios",
    x = expression(theta[true]),
    y = "Bias",
    color = "Parameter",
    linetype = "m0:m1"
  ) +
  scale_x_continuous(breaks = SIM_CONFIG$thetas) +
  plot_base_theme()
save_plot(p_q_bias_diag_agg, "q_bias_diag_aggregate.png", width = 16, height = 10)

# =============================================================================
# OFF-DIAGONAL PLOTS: (q0=0.6, q1=0.9) and (q0=0.9, q1=0.6)
# =============================================================================

off_diag_cases <- list(
  list(q0 = 0.6, q1 = 0.9, label = "q0=0.6_q1=0.9"),
  list(q0 = 0.9, q1 = 0.6, label = "q0=0.9_q1=0.6")
)

for (case in off_diag_cases) {

  ci_summary_offdiag <- ci_summary_plot %>%
    filter(q0 == case$q0, q1 == case$q1) %>%
    mutate(
      label_ratio_lab = factor(
        format_label_ratio(label_ratio),
        levels = format_label_ratio(sort(unique(label_ratio)))
      )
    )

  q_bias_offdiag <- q_bias_long %>%
    filter(q0 == case$q0, q1 == case$q1) %>%
    mutate(
      label_ratio_lab = factor(
        format_label_ratio(label_ratio),
        levels = format_label_ratio(sort(unique(label_ratio)))
      )
    )

  filename_suffix <- paste0("_offdiag_", case$label, ".png")
  subtitle_text <- paste0("q0 = ", case$q0, ", q1 = ", case$q1, " | All m0:m1 ratios")

  # Coverage gap plot (off-diagonal)
  p_cov_offdiag <- ggplot(ci_summary_offdiag,
                          aes(x = theta_true, y = coverage_gap,
                              color = method, linetype = split_label,
                              group = interaction(method, split_label))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ label_ratio_lab, nrow = 1) +
    labs(
      title = paste0("Coverage - Nominal (q0 = ", case$q0, ", q1 = ", case$q1, ")"),
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Coverage gap",
      color = "Method",
      linetype = "m0:m1"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_cov_offdiag, paste0("coverage_gap", filename_suffix), width = 18, height = 5)

  # CI width plot (off-diagonal)
  p_width_offdiag <- ggplot(ci_summary_offdiag,
                            aes(x = theta_true, y = ci_width,
                                color = method, linetype = split_label,
                                group = interaction(method, split_label))) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ label_ratio_lab, nrow = 1) +
    labs(
      title = paste0("CI Width (q0 = ", case$q0, ", q1 = ", case$q1, ")"),
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Mean CI width",
      color = "Method",
      linetype = "m0:m1"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_width_offdiag, paste0("ci_width", filename_suffix), width = 18, height = 5)

  # Bias plot (off-diagonal)
  p_bias_offdiag <- ggplot(ci_summary_offdiag,
                           aes(x = theta_true, y = bias,
                               color = method, linetype = split_label,
                               group = interaction(method, split_label))) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ label_ratio_lab, nrow = 1) +
    labs(
      title = paste0("Estimator Bias (q0 = ", case$q0, ", q1 = ", case$q1, ")"),
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Bias",
      color = "Method",
      linetype = "m0:m1"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_bias_offdiag, paste0("bias", filename_suffix), width = 18, height = 5)

  # Bias percent plot (off-diagonal)
  p_bias_pct_offdiag <- ggplot(ci_summary_offdiag,
                               aes(x = theta_true, y = bias_pct,
                                   color = method, linetype = split_label,
                                   group = interaction(method, split_label))) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ label_ratio_lab, nrow = 1) +
    labs(
      title = paste0("Percent Bias (q0 = ", case$q0, ", q1 = ", case$q1, ")"),
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Bias (%)",
      color = "Method",
      linetype = "m0:m1"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_bias_pct_offdiag, paste0("bias_pct", filename_suffix), width = 18, height = 5)

  # Q bias plot (off-diagonal)
  p_q_bias_offdiag <- ggplot(q_bias_offdiag,
                             aes(x = theta_true, y = bias,
                                 color = metric, linetype = split_label,
                                 group = interaction(metric, split_label))) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ label_ratio_lab, nrow = 1) +
    labs(
      title = paste0("Calibration Estimate Bias (q0 = ", case$q0, ", q1 = ", case$q1, ")"),
      subtitle = subtitle_text,
      x = expression(theta[true]),
      y = "Bias",
      color = "Parameter",
      linetype = "m0:m1"
    ) +
    scale_x_continuous(breaks = SIM_CONFIG$thetas) +
    plot_base_theme()
  save_plot(p_q_bias_offdiag, paste0("q_bias", filename_suffix), width = 18, height = 5)
}

message("Results saved to ", results_dir)
