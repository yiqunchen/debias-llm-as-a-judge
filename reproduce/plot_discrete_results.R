# Plotting Functions for Discrete Prediction Simulation
# Line plots with mu_3 on x-axis, faceted by label_ratio

library(tidyverse)

# ============================================================================
# COLOR SCHEMES
# ============================================================================

METHOD_COLORS <- c(
  "Naive" = "#525252",
  "PPI" = "#1B9E77",
  "PPI++" = "#D95F02",
  "EIF" = "#E7298A",
  "EIF-linear" = "#7570B3",
  "EIF-gam" = "#E6AB02",
  "EIF-spline" = "#66A61E"
)

# ============================================================================
# MAIN PLOTTING FUNCTIONS
# ============================================================================

#' Plot coverage vs mu_3
#' Faceted by label_ratio (rows)
plot_coverage_discrete <- function(summary_df) {
  ggplot(summary_df, aes(x = mu_3, y = coverage, color = method)) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray50") +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), alpha = 0.7) +
    facet_wrap(~ labeled_pct, ncol = 1) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = expression(mu[3]),
      y = "Coverage",
      title = "CI Coverage vs Bias Magnitude",
      subtitle = expression("DGP: Z ~ Unif{1,2,3}, Y|Z ~ N(" * mu[Z] * ", 1), " * hat(Y) * " = Z"),
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Plot bias vs mu_3
plot_bias_discrete <- function(summary_df) {
  ggplot(summary_df, aes(x = mu_3, y = bias, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), alpha = 0.7) +
    facet_wrap(~ labeled_pct, ncol = 1) +
    scale_color_manual(values = METHOD_COLORS) +
    labs(
      x = expression(mu[3]),
      y = "Bias",
      title = "Estimation Bias vs Bias Magnitude",
      subtitle = expression("Naive bias = 2 - E[Y] = (3 - " * mu[3] * ")/3"),
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Plot RMSE vs mu_3
plot_rmse_discrete <- function(summary_df) {
  ggplot(summary_df, aes(x = mu_3, y = rmse, color = method)) +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), alpha = 0.7) +
    facet_wrap(~ labeled_pct, ncol = 1) +
    scale_color_manual(values = METHOD_COLORS) +
    labs(
      x = expression(mu[3]),
      y = "RMSE",
      title = "Root Mean Squared Error vs Bias Magnitude",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

#' Plot CI width (excluding Naive)
plot_ci_width_discrete <- function(summary_df) {
  df_no_naive <- summary_df %>% filter(method != "Naive")
  colors_no_naive <- METHOD_COLORS[names(METHOD_COLORS) != "Naive"]

  ggplot(df_no_naive, aes(x = mu_3, y = mean_ci_width, color = method)) +
    geom_point(size = 2, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), alpha = 0.7) +
    facet_wrap(~ labeled_pct, ncol = 1) +
    scale_color_manual(values = colors_no_naive) +
    labs(
      x = expression(mu[3]),
      y = "Mean CI Width",
      title = "CI Width vs Bias Magnitude (excluding Naive)",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# ============================================================================
# Y vs Yhat CALIBRATION PLOT
# ============================================================================

CALIBRATION_COLORS <- c(
  "Per-category" = "#E7298A",
  "Linear" = "#7570B3",
  "GAM" = "#E6AB02",
  "Spline" = "#66A61E"
)

#' Plot Y vs Yhat with calibration fits for discrete predictions
#'
#' Shows scatter of Y vs discrete Yhat (1, 2, 3) with different calibration fits
#' Faceted by mu_3 values to show how bias magnitude affects calibration
#'
#' @param N Sample size for visualization
#' @param mu_3_values Which mu_3 values to show (subset for clarity)
#' @param seed Random seed
plot_calibration_discrete <- function(N = 1000, mu_3_values = c(3, 6, 9), seed = 42) {
  source(here::here("R/sim_discrete.R"))
  set.seed(seed)

  plots <- list()

  for (mu_3 in mu_3_values) {
    mu <- c(1, 2, mu_3)

    # Generate data
    dgp <- generate_discrete_dgp(N, mu = mu, sigma = 1)
    Y <- dgp$Y
    Yhat <- dgp$Yhat

    df <- data.frame(Y = Y, Yhat = Yhat)

    # Fit calibration models
    # Per-category means
    cat_means <- df %>%
      group_by(Yhat) %>%
      summarize(Y_pred = mean(Y), .groups = "drop")

    # Linear
    fit_linear <- lm(Y ~ Yhat, data = df)

    # GAM - use bs() instead of s() for discrete data
    fit_gam <- tryCatch({
      mgcv::gam(Y ~ s(Yhat, k = 3), data = df)
    }, error = function(e) NULL)

    # Spline - use polynomial for 3 points
    fit_spline <- tryCatch({
      lm(Y ~ poly(Yhat, 2), data = df)
    }, error = function(e) NULL)

    # Create prediction data for smooth curves
    # Use fine grid for smooth appearance
    yhat_grid <- seq(0.8, 3.2, length.out = 100)
    pred_df <- data.frame(Yhat = yhat_grid)

    pred_df$Linear <- predict(fit_linear, pred_df)
    pred_df$GAM <- if (!is.null(fit_gam)) predict(fit_gam, pred_df) else pred_df$Linear
    pred_df$Spline <- if (!is.null(fit_spline)) predict(fit_spline, pred_df) else pred_df$Linear

    pred_long <- pred_df %>%
      pivot_longer(cols = c(Linear, GAM, Spline),
                   names_to = "method", values_to = "Y_pred")

    # Calculate naive bias
    naive_bias <- mean(Yhat) - mean(Y)

    # Create plot
    p <- ggplot() +
      # Jittered points
      geom_jitter(data = df, aes(x = Yhat, y = Y),
                  alpha = 0.15, width = 0.15, height = 0, size = 0.8,
                  color = "gray40") +
      # Per-category means as points
      geom_point(data = cat_means, aes(x = Yhat, y = Y_pred, color = "Per-category"),
                 size = 4, shape = 18) +
      # Calibration curves
      geom_line(data = pred_long, aes(x = Yhat, y = Y_pred, color = method),
                linewidth = 1) +
      # Perfect calibration line (Y = Yhat)
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      scale_color_manual(values = CALIBRATION_COLORS) +
      scale_x_continuous(breaks = 1:3, limits = c(0.5, 3.5)) +
      labs(
        x = expression(hat(Y) ~ "(discrete category)"),
        y = "Y",
        title = bquote(mu == "(" * 1 * ", " * 2 * ", " * .(mu_3) * ")"),
        subtitle = sprintf("Naive bias: %.2f", naive_bias),
        color = "Calibration"
      ) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 11),
            plot.subtitle = element_text(size = 9))

    plots[[as.character(mu_3)]] <- p
  }

  # Combine plots
  library(patchwork)

  combined <- wrap_plots(plots, ncol = length(mu_3_values), guides = "collect") +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = expression("Y vs " * hat(Y) * " with Calibration Fits (Discrete Predictions)"),
      subtitle = expression("DGP: Z ~ Unif{1,2,3}, Y|Z ~ N(" * mu[Z] * ", 1), " * hat(Y) * " = Z. Dashed line = perfect calibration.")
    )

  combined
}

# ============================================================================
# COMBINED COVERAGE + CI WIDTH PLOT (for publication)
# ============================================================================

#' Combined coverage and CI width plot
#' Side-by-side with shared legend, bigger fonts, no subtitles
#' Only shows 5% and 20% labeling ratios
plot_combined_coverage_ciwidth <- function(summary_df) {
  library(patchwork)

  # Filter to 5% and 20% only, order 5% before 20%
  summary_df <- summary_df %>%
    filter(labeled_pct %in% c("Labeled: 5%", "Labeled: 20%")) %>%
    mutate(labeled_pct = factor(labeled_pct, levels = c("Labeled: 5%", "Labeled: 20%")))

  # Common theme with bigger fonts
  big_font_theme <- theme_bw(base_size = 18) +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      strip.text = element_text(size = 16, hjust = 0.5),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      plot.title = element_text(size = 20, hjust = 0.5),
      legend.position = "bottom"
    )

  # Coverage plot (WITH legend - includes Naive)
  p_cov <- ggplot(summary_df, aes(x = mu_3, y = coverage, color = method)) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray50") +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), alpha = 0.7, linewidth = 1) +
    facet_wrap(~ labeled_pct, ncol = 1) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      x = expression(mu[3]),
      y = "Coverage",
      title = "Coverage",
      color = "Method"
    ) +
    big_font_theme

  # CI width plot (excluding Naive, NO legend)
  df_no_naive <- summary_df %>% filter(method != "Naive")
  colors_no_naive <- METHOD_COLORS[names(METHOD_COLORS) != "Naive"]

  p_width <- ggplot(df_no_naive, aes(x = mu_3, y = mean_ci_width, color = method)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    geom_line(position = position_dodge(width = 0.3), alpha = 0.7, linewidth = 1) +
    facet_wrap(~ labeled_pct, ncol = 1, scales = "free_y") +
    scale_color_manual(values = colors_no_naive) +
    labs(
      x = expression(mu[3]),
      y = "Mean CI Width",
      title = "CI Width",
      color = "Method"
    ) +
    big_font_theme +
    theme(legend.position = "none")  # Hide legend on CI width plot

  # Combine - legend from coverage plot at bottom
  combined <- p_cov + p_width

  combined
}

# ============================================================================
# GENERATE ALL PLOTS
# ============================================================================

generate_discrete_plot_suite <- function(summary_df, output_dir = "results/discrete_plots") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  cat("Generating discrete simulation plots...\n")

  # Coverage
  p <- plot_coverage_discrete(summary_df)
  ggsave(file.path(output_dir, "coverage.pdf"), p, width = 8, height = 10, dpi = 300)
  cat("  Saved coverage.pdf\n")

  # Bias
  p <- plot_bias_discrete(summary_df)
  ggsave(file.path(output_dir, "bias.pdf"), p, width = 8, height = 10, dpi = 300)
  cat("  Saved bias.pdf\n")

  # RMSE
  p <- plot_rmse_discrete(summary_df)
  ggsave(file.path(output_dir, "rmse.pdf"), p, width = 8, height = 10, dpi = 300)
  cat("  Saved rmse.pdf\n")

  # CI Width
  p <- plot_ci_width_discrete(summary_df)
  ggsave(file.path(output_dir, "ci_width.pdf"), p, width = 8, height = 10, dpi = 300)
  cat("  Saved ci_width.pdf\n")

  # Combined coverage + CI width (for publication)
  p <- plot_combined_coverage_ciwidth(summary_df)
  ggsave(file.path(output_dir, "coverage_ciwidth_combined.pdf"), p, width = 14, height = 8, dpi = 300)
  cat("  Saved coverage_ciwidth_combined.pdf\n")

  # Y vs Yhat calibration plot
  cat("Generating Y vs Yhat calibration plot...\n")
  p <- plot_calibration_discrete()
  ggsave(file.path(output_dir, "calibration_fits.pdf"), p, width = 12, height = 5, dpi = 300)
  cat("  Saved calibration_fits.pdf\n")

  cat(sprintf("\nAll plots saved to %s/\n", output_dir))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if (sys.nframe() == 0) {
  library(here)
  library(mgcv)
  library(splines)
  library(patchwork)

  summary_file <- "reproduce/discrete_sim_summary.csv"

  if (!file.exists(summary_file)) {
    stop("Summary file not found. Run simulation_discrete.R first.")
  }

  summary_df <- read_csv(summary_file, show_col_types = FALSE)
  cat("Loaded summary with", nrow(summary_df), "rows\n")

  generate_discrete_plot_suite(summary_df)

  cat("\nAll plots generated!\n")
}
