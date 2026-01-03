# Plotting Functions for Continuous Score Simulation
# Facet by surrogate quality (columns) and label_ratio (rows)

library(tidyverse)

# ============================================================================
# COLOR SCHEMES (matching binary outcome case)
# ============================================================================

METHOD_COLORS <- c(
  "Naive" = "#525252",
  "PPI" = "#1B9E77",
  "PPI++" = "#D95F02",
  "EIF-linear" = "#E6AB02",
  "EIF-gam" = "#E7298A",
  "EIF-spline" = "#7570B3"
)

CALIBRATION_COLORS <- c(
  "linear" = "#E6AB02",
  "gam" = "#E7298A",
  "spline" = "#7570B3"
)

# ============================================================================
# MAIN PLOTTING FUNCTIONS (faceted by surrogate quality and label_ratio)
# ============================================================================

#' Plot coverage vs N_total
#' Faceted by model_quality (cols) x label_ratio (rows)
plot_coverage <- function(summary_df) {
  ggplot(summary_df, aes(x = N_total, y = coverage, color = method)) +
    geom_hline(yintercept = 0.90, linetype = "dashed", color = "gray50") +
    geom_point(size = 2, position = position_dodge(width = 0.1)) +
    geom_line(position = position_dodge(width = 0.1), alpha = 0.7) +
    facet_grid(label_ratio ~ model_quality,
               labeller = labeller(label_ratio = label_both,
                                   model_quality = label_value)) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_x_continuous(trans = "log10") +
    labs(
      x = "N (Total Sample Size)",
      y = "Coverage",
      title = "CI Coverage (no covariate shift)",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10))
}

#' Plot bias vs N_total
plot_bias <- function(summary_df) {
  ggplot(summary_df, aes(x = N_total, y = bias, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 2, position = position_dodge(width = 0.1)) +
    geom_line(position = position_dodge(width = 0.1), alpha = 0.7) +
    facet_grid(label_ratio ~ model_quality,
               labeller = labeller(label_ratio = label_both,
                                   model_quality = label_value)) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_x_continuous(trans = "log10") +
    labs(
      x = "N (Total Sample Size)",
      y = "Bias",
      title = "Estimation Bias (no covariate shift)",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10))
}

#' Plot RMSE vs N_total
plot_rmse <- function(summary_df) {
  ggplot(summary_df, aes(x = N_total, y = rmse, color = method)) +
    geom_point(size = 2, position = position_dodge(width = 0.1)) +
    geom_line(position = position_dodge(width = 0.1), alpha = 0.7) +
    facet_grid(label_ratio ~ model_quality,
               labeller = labeller(label_ratio = label_both,
                                   model_quality = label_value)) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_x_continuous(trans = "log10") +
    labs(
      x = "N (Total Sample Size)",
      y = "RMSE",
      title = "Root Mean Squared Error (no covariate shift)",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10))
}

#' Plot CI width vs N_total
#' Excludes Naive since it doesn't maintain nominal coverage
plot_ci_width <- function(summary_df) {
  # Filter out Naive - its CI width is not meaningful since coverage is poor
  df_no_naive <- summary_df %>% filter(method != "Naive")

  # Colors without Naive
  colors_no_naive <- METHOD_COLORS[names(METHOD_COLORS) != "Naive"]

  ggplot(df_no_naive, aes(x = N_total, y = mean_ci_width, color = method)) +
    geom_point(size = 2, position = position_dodge(width = 0.1)) +
    geom_line(position = position_dodge(width = 0.1), alpha = 0.7) +
    facet_grid(label_ratio ~ model_quality,
               labeller = labeller(label_ratio = label_both,
                                   model_quality = label_value)) +
    scale_color_manual(values = colors_no_naive) +
    scale_x_continuous(trans = "log10") +
    labs(
      x = "N (Total Sample Size)",
      y = "Mean CI Width",
      title = "CI Width (excluding Naive, no covariate shift)",
      color = "Method"
    ) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10))
}

# ============================================================================
# Y vs Yhat CALIBRATION PLOT
# ============================================================================

#' Generate Y vs Yhat scatter plot with calibration fits
#'
#' Shows one panel per model quality (no covariate shift)
#' All data from same distribution X ~ N(mu_x, 1)
#'
#' @param N Total sample size
#' @param train_frac Fraction for model training
#' @param mu_x Mean of X for all data
#' @param seed Random seed
plot_calibration_fits <- function(N = 2000, train_frac = 0.20, mu_x = 1, seed = 42) {
  # Source the DGP functions
  source(here::here("R/sim_continuous.R"))

  set.seed(seed)

  plots <- list()

  for (qual in c("good", "medium", "poor-OLS")) {
    # Generate ALL data from same distribution
    dgp <- generate_continuous_dgp(N, mu_x = mu_x)
    X <- dgp$X
    Y <- dgp$Y

    # Split into train/eval
    N_train <- floor(train_frac * N)
    train_idx <- sample(N, N_train)
    eval_idx <- setdiff(seq_len(N), train_idx)

    X_train <- X[train_idx, , drop = FALSE]
    Y_train <- Y[train_idx]
    X_eval <- X[eval_idx, , drop = FALSE]
    Y_eval <- Y[eval_idx]

    # Train model
    train_fn <- switch(qual,
      "good" = train_model_good,
      "medium" = train_model_medium,
      "poor-OLS" = train_model_poor_ols
    )
    model <- train_fn(X_train, Y_train)
    Yhat_eval <- model(X_eval)

    # Create data frame for plotting
    df <- data.frame(Y = Y_eval, Yhat = Yhat_eval)

    # Fit calibration models
    fit_linear <- lm(Y ~ Yhat, data = df)
    fit_gam <- tryCatch(
      mgcv::gam(Y ~ s(Yhat, k = 10), data = df),
      error = function(e) lm(Y ~ Yhat, data = df)
    )
    fit_spline <- tryCatch(
      lm(Y ~ splines::ns(Yhat, df = 4), data = df),
      error = function(e) lm(Y ~ Yhat, data = df)
    )

    # Create prediction grid
    yhat_grid <- seq(min(Yhat_eval), max(Yhat_eval), length.out = 100)
    pred_df <- data.frame(Yhat = yhat_grid)

    pred_df$linear <- predict(fit_linear, pred_df)
    pred_df$gam <- predict(fit_gam, pred_df)
    pred_df$spline <- predict(fit_spline, pred_df)

    # Pivot to long format for plotting
    pred_long <- pred_df %>%
      pivot_longer(cols = c(linear, gam, spline),
                   names_to = "method", values_to = "Y_pred")

    # Calculate bias for subtitle
    bias <- mean(Yhat_eval) - mean(Y_eval)

    # Create plot
    p <- ggplot() +
      geom_point(data = df, aes(x = Yhat, y = Y), alpha = 0.3, size = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(data = pred_long, aes(x = Yhat, y = Y_pred, color = method),
                linewidth = 1) +
      scale_color_manual(values = CALIBRATION_COLORS) +
      labs(
        x = expression(hat(Y)),
        y = "Y",
        title = qual,
        subtitle = sprintf("Naive bias: %.2f", bias),
        color = "Calibration"
      ) +
      theme_bw() +
      theme(legend.position = "bottom",
            plot.title = element_text(size = 11),
            plot.subtitle = element_text(size = 9))

    plots[[qual]] <- p
  }

  # Combine plots: 1 row x 3 cols
  library(patchwork)

  combined <- (plots[["good"]] + plots[["medium"]] + plots[["poor-OLS"]]) +
    plot_layout(ncol = 3, guides = "collect")

  combined + plot_annotation(
    title = sprintf("Y vs Yhat with Calibration Fits (N=%d, X ~ N(%d,1))", N, mu_x),
    subtitle = "Dashed line = perfect calibration (Y = Yhat)"
  )
}

# ============================================================================
# GENERATE ALL PLOTS
# ============================================================================

#' Generate all plots and save to disk
#'
#' @param summary_df Summary data frame
#' @param output_dir Output directory for plots
generate_continuous_plot_suite <- function(summary_df, output_dir = "results/continuous_plots") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  cat("Generating main simulation plots...\n")

  # Coverage
  p <- plot_coverage(summary_df)
  ggsave(file.path(output_dir, "coverage.pdf"), p, width = 12, height = 10, dpi = 300)
  cat("  Saved coverage.pdf\n")

  # Bias
  p <- plot_bias(summary_df)
  ggsave(file.path(output_dir, "bias.pdf"), p, width = 12, height = 10, dpi = 300)
  cat("  Saved bias.pdf\n")

  # RMSE
  p <- plot_rmse(summary_df)
  ggsave(file.path(output_dir, "rmse.pdf"), p, width = 12, height = 10, dpi = 300)
  cat("  Saved rmse.pdf\n")

  # CI Width
  p <- plot_ci_width(summary_df)
  ggsave(file.path(output_dir, "ci_width.pdf"), p, width = 12, height = 10, dpi = 300)
  cat("  Saved ci_width.pdf\n")

  # Calibration fits plot (1x3: one per model quality)
  cat("Generating calibration fits plot...\n")
  p <- plot_calibration_fits()
  ggsave(file.path(output_dir, "calibration_fits.pdf"), p, width = 14, height = 5, dpi = 300)
  cat("  Saved calibration_fits.pdf\n")

  cat(sprintf("\nAll plots saved to %s/\n", output_dir))
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if (sys.nframe() == 0) {
  library(here)
  library(patchwork)
  library(mgcv)
  library(splines)
  library(randomForest)
  library(rpart)

  summary_file <- "reproduce/continuous_sim_summary.csv"

  if (!file.exists(summary_file)) {
    stop("Summary file not found. Run simulation_continuous_score.R first.")
  }

  summary_df <- read_csv(summary_file, show_col_types = FALSE)
  cat("Loaded summary with", nrow(summary_df), "rows\n")

  generate_continuous_plot_suite(summary_df)

  cat("\nAll plots generated!\n")
}
