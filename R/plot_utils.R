# Plot utilities for simulation results
# Reusable plotting functions for coverage, CI width, bias, and calibration plots

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "theta_true",
    "bias",
    "bias_pct",
    "ci_width",
    "coverage",
    "method",
    "metric"
  ))
}

# Colorblind-friendly palette (fixed mapping for methods)
# Uses ColorBrewer "Dark2" palette - highly distinguishable and publication-safe
METHOD_COLORS <- c(
  "PPI" = "#1B9E77",
  "PPI++" = "#D95F02",
  "Rogan-Gladen" = "#7570B3",
  "Joint MLE" = "#E7298A",
  "Naive" = "#525252",
  "EIF" = "#E6AB02",
  "q0" = "#1B9E77",
  "q1" = "#D95F02"
)

#' Get the fixed color scale for methods
#' @param methods Character vector of method names to include
#' @return ggplot2 scale_color_manual
#' @export
get_method_color_scale <- function(methods = NULL) {
  if (is.null(methods)) {
    scale_color_manual(values = METHOD_COLORS, drop = FALSE)
  } else {
    scale_color_manual(values = METHOD_COLORS[methods], drop = FALSE)
  }
}

#' Base theme for standard plots
#' @return A ggplot2 theme object
#' @export
plot_base_theme <- function() {
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

#' Larger legend theme for aggregate plots (stacked in two rows)
#' @return A ggplot2 theme object
#' @export
plot_base_theme_aggregate <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 16),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 14),
      axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.size = unit(1.5, "cm"),
      legend.key.width = unit(3, "cm"),
      legend.spacing.y = unit(0.2, "cm")
    )
}

#' Format a labeling ratio as a percentage string
#'
#' @param x Numeric vector of label fractions (between 0 and 1).
#' @return A character vector like "5% labeled".
#' @export
format_label_ratio <- function(x) paste0(scales::percent(x, accuracy = 1), " labeled")


#' Generic line plot builder
#'
#' @param data Data frame
#' @param x_var X variable (unquoted)
#' @param y_var Y variable (unquoted)
#' @param color_var Color variable (unquoted)
#' @param linetype_var Optional linetype variable (unquoted), NULL for none
#' @param facet_formula Facet formula (e.g., \code{q0_lab ~ q1_lab})
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param y_lab Y-axis label
#' @param color_lab Color legend label
#' @param linetype_lab Linetype legend label (if applicable)
#' @param hline Optional horizontal line y-intercept
#' @param hline_style "dashed" or "dotted"
#' @param x_breaks X-axis breaks
#' @param y_limits Optional y-axis limits (c(min, max))
#' @param use_aggregate_theme Use larger theme for aggregate plots
#' @param point_size Point size
#' @param line_width Line width
#' @return A ggplot object
#' @export
build_line_plot <- function(data,
                            x_var,
                            y_var,
                            color_var,
                            linetype_var = NULL,
                            facet_formula,
                            title,
                            subtitle = "",
                            y_lab,
                            color_lab = "Method",
                            linetype_lab = "neg:pos",
                            hline = NULL,
                            hline_style = "dotted",
                            x_breaks = NULL,
                            y_limits = NULL,
                            use_aggregate_theme = FALSE,
                            point_size = 1.5,
                            line_width = 1) {

  x_var <- rlang::enquo(x_var)
  y_var <- rlang::enquo(y_var)
  color_var <- rlang::enquo(color_var)
  # linetype_var is passed as a string (column name) or NULL
  use_linetype <- !is.null(linetype_var) && is.character(linetype_var)
  if (use_linetype) {
    linetype_sym <- rlang::sym(linetype_var)
  }

  # Build aesthetic mapping
  if (!use_linetype) {
    p <- ggplot(data, aes(x = !!x_var, y = !!y_var,
                          color = !!color_var, group = !!color_var))
  } else {
    p <- ggplot(data, aes(x = !!x_var, y = !!y_var,
                          color = !!color_var, linetype = !!linetype_sym,
                          group = interaction(!!color_var, !!linetype_sym)))
  }

  # Add horizontal line if specified
  if (!is.null(hline)) {
    hline_lt <- if (hline_style == "dashed") "dashed" else "dotted"
    p <- p + geom_hline(yintercept = hline, linetype = hline_lt, color = "grey40")
  }

  # Add geoms
  p <- p +
    geom_line(linewidth = line_width) +
    geom_point(size = point_size)

  # Add faceting
  if (inherits(facet_formula, "formula")) {
    # Check if it's a one-sided formula (facet_wrap style)
    if (length(facet_formula) == 2) {
      p <- p + facet_wrap(facet_formula, nrow = 1)
    } else {
      p <- p + facet_grid(facet_formula)
    }
  }

  # Build labels
  labs_args <- list(
    title = title,
    subtitle = subtitle,
    x = expression(theta[true]),
    y = y_lab,
    color = color_lab
  )
  if (use_linetype) {
    labs_args$linetype <- linetype_lab
  }
  p <- p + do.call(labs, labs_args)

  # Add scales
  if (!is.null(x_breaks)) {
    p <- p + scale_x_continuous(breaks = x_breaks)
  }
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits)
  }

  # Add fixed colorblind-friendly color scale
  # This ensures consistent colors even when methods are filtered out
  p <- p + scale_color_manual(values = METHOD_COLORS, drop = FALSE)

  # Add theme
  p <- p + if (use_aggregate_theme) plot_base_theme_aggregate() else plot_base_theme()

  p
}

#' Plot empirical coverage across simulation settings
#'
#' @param data A data frame containing summary statistics per configuration.
#'   Must include columns such as \code{theta_true}, \code{coverage},
#'   and \code{method}.
#' @param facet_formula A faceting formula (e.g. \code{q0_lab ~ q1_lab})
#'   specifying how to arrange simulation panels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param linetype_var Optional name of a column in \code{data} that controls
#'   line type. If \code{NULL}, all lines are solid.
#' @param use_aggregate Logical; whether to use a larger theme suitable for
#'   aggregate plots with many legend entries.
#' @param x_breaks Optional numeric vector of breaks for the x-axis (in
#'   \code{theta_true} space).
#' @param point_size Point size for the scatter overlay.
#' @param y_limits Optional y-axis limits (c(min, max)).
#' @param alpha Nominal coverage level for the horizontal reference line.
#'   Defaults to 0.9.
#'
#' @return A \code{ggplot} object showing empirical coverage as a function of
#'   \code{theta_true}, with a horizontal reference line at \code{alpha}.
#' @export
plot_coverage <- function(data, facet_formula, title = "Coverage", subtitle = "",
                          linetype_var = NULL, use_aggregate = FALSE,
                          x_breaks = NULL, point_size = 1.5, y_limits = NULL,
                          alpha = 0.9) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = coverage,
    color_var = method,
    linetype_var = linetype_var,
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Coverage",
    hline = alpha,
    hline_style = "dashed",
    x_breaks = x_breaks,
    y_limits = y_limits,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Plot mean confidence interval width across simulation settings
#'
#' @param data A data frame containing summary statistics per configuration.
#'   Must include columns such as \code{theta_true}, \code{ci_width},
#'   and \code{method}.
#' @param facet_formula A faceting formula (e.g. \code{q0_lab ~ q1_lab})
#'   specifying how to arrange simulation panels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param linetype_var Optional name of a column in \code{data} that controls
#'   line type. If \code{NULL}, all lines are solid.
#' @param use_aggregate Logical; whether to use a larger theme suitable for
#'   aggregate plots with many legend entries.
#' @param x_breaks Optional numeric vector of breaks for the x-axis (in
#'   \code{theta_true} space).
#' @param point_size Point size for the scatter overlay.
#' @param y_limits Optional numeric vector of length 2 giving y-axis limits.
#'
#' @return A \code{ggplot} object showing mean confidence interval width
#'   as a function of \code{theta_true}.
#' @export
plot_ci_width <- function(data, facet_formula, title = "CI Width", subtitle = "",
                          linetype_var = NULL, use_aggregate = FALSE,
                          x_breaks = NULL, point_size = 1.5, y_limits = NULL) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = ci_width,
    color_var = method,
    linetype_var = linetype_var,
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Mean CI width",
    x_breaks = x_breaks,
    y_limits = y_limits,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Plot estimator bias across simulation settings
#'
#' @param data A data frame containing summary statistics per configuration.
#' @param facet_formula A faceting formula (e.g. \code{q0_lab ~ q1_lab}).
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param linetype_var Optional name of a column in \code{data} that controls line type.
#'   If \code{NULL}, all lines are solid.
#' @param use_aggregate Logical; whether to use a larger theme suitable for
#'   aggregate plots with many legend entries.
#' @param x_breaks Optional numeric vector of breaks for the x-axis.
#' @param point_size Point size for the scatter overlay.
#' @return A \code{ggplot} object showing bias curves.
#' @export
plot_bias <- function(data, facet_formula, title = "Estimator Bias", subtitle = "",
                      linetype_var = NULL, use_aggregate = FALSE,
                      x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = bias,
    color_var = method,
    linetype_var = linetype_var,
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Bias",
    hline = 0,
    hline_style = "dotted",
    x_breaks = x_breaks,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Plot percent bias across simulation settings
#'
#' @param data A data frame containing summary statistics per configuration.
#'   Must include columns such as \code{theta_true}, \code{bias_pct},
#'   and \code{method}.
#' @param facet_formula A faceting formula (e.g. \code{q0_lab ~ q1_lab})
#'   specifying how to arrange simulation panels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param linetype_var Optional name of a column in \code{data} that controls
#'   line type. If \code{NULL}, all lines are solid.
#' @param use_aggregate Logical; whether to use a larger theme suitable for
#'   aggregate plots with many legend entries.
#' @param x_breaks Optional numeric vector of breaks for the x-axis (in
#'   \code{theta_true} space).
#' @param point_size Point size for the scatter overlay.
#'
#' @return A \code{ggplot} object showing percent bias curves.
#' @export
plot_bias_pct <- function(data, facet_formula, title = "Percent Bias", subtitle = "",
                          linetype_var = NULL, use_aggregate = FALSE,
                          x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = bias_pct,
    color_var = method,
    linetype_var = linetype_var,
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Bias (%)",
    hline = 0,
    hline_style = "dotted",
    x_breaks = x_breaks,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Plot calibration bias for \eqn{\hat q_0} and \eqn{\hat q_1}
#'
#' @param data A data frame containing calibration summary statistics,
#'   typically with columns \code{theta_true}, \code{bias},
#'   \code{metric} (e.g. "q0" or "q1"), and any faceting variables.
#' @param facet_formula A faceting formula (e.g. \code{q0_lab ~ q1_lab}
#'   or \code{q_val_lab ~ label_ratio_lab}) specifying how panels are arranged.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param linetype_var Optional name of a column in \code{data} that controls
#'   line type. If \code{NULL}, all lines are solid.
#' @param use_aggregate Logical; whether to use a larger theme suitable for
#'   aggregate plots with many legend entries.
#' @param x_breaks Optional numeric vector of breaks for the x-axis (in
#'   \code{theta_true} space).
#' @param point_size Point size for the scatter overlay.
#'
#' @return A \code{ggplot} object showing calibration bias for
#'   \eqn{\hat q_0} and \eqn{\hat q_1}.
#' @export
plot_q_bias <- function(data, facet_formula, title = "Calibration Estimate Bias", subtitle = "",
                        linetype_var = NULL, use_aggregate = FALSE,
                        x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = bias,
    color_var = metric,
    linetype_var = linetype_var,
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Bias",
    color_lab = "Parameter",
    hline = 0,
    hline_style = "dotted",
    x_breaks = x_breaks,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Generate all 5 standard plots for a given data subset
#'
#' @param ci_data CI summary data (filtered)
#' @param q_data Q bias data (filtered)
#' @param facet_formula Facet formula
#' @param filename_suffix Suffix for filenames
#' @param subtitle Plot subtitle
#' @param save_fn Function to save plots (takes plot and filename)
#' @param x_breaks X-axis breaks
#' @param linetype_var Optional linetype variable for aggregate plots
#' @param use_aggregate Use aggregate theme
#' @param width Plot width
#' @param height Plot height
#' @param point_size Point size
#' @return NULL (saves plots as side effect)
#' @export
generate_plot_suite <- function(ci_data,
                                q_data,
                                facet_formula,
                                filename_suffix,
                                subtitle = "",
                                save_fn,
                                x_breaks = NULL,
                                linetype_var = NULL,
                                use_aggregate = FALSE,
                                width = 14,
                                height = 12,
                                point_size = 1.5) {

  title_suffix <- if (use_aggregate) "" else ""

  # Coverage
  p_cov <- plot_coverage(ci_data, facet_formula, "Coverage", subtitle,
                         linetype_var, use_aggregate, x_breaks, point_size)
  save_fn(p_cov, paste0("coverage", filename_suffix), width, height)

  # CI Width
  p_width <- plot_ci_width(ci_data, facet_formula, "CI Width", subtitle,
                           linetype_var, use_aggregate, x_breaks, point_size)
  save_fn(p_width, paste0("ci_width", filename_suffix), width, height)

  # Bias
  p_bias <- plot_bias(ci_data, facet_formula, "Estimator Bias", subtitle,
                      linetype_var, use_aggregate, x_breaks, point_size)
  save_fn(p_bias, paste0("bias", filename_suffix), width, height)

  # Percent Bias
  p_bias_pct <- plot_bias_pct(ci_data, facet_formula, "Percent Bias", subtitle,
                              linetype_var, use_aggregate, x_breaks, point_size)
  save_fn(p_bias_pct, paste0("bias_pct", filename_suffix), width, height)

  # Q Bias
  p_q_bias <- plot_q_bias(q_data, facet_formula, "Calibration Estimate Bias", subtitle,
                          linetype_var, use_aggregate, x_breaks, point_size)
  save_fn(p_q_bias, paste0("q_bias", filename_suffix), width, height)
}
