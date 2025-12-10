# Plot utilities for simulation results
# Reusable plotting functions for coverage, CI width, bias, and calibration plots

#' Base theme for standard plots
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

#' Format label ratio as percentage string
format_label_ratio <- function(x) {
  paste0(as.integer(x * 100), "% labeled")
}

#' Generic line plot builder
#'
#' @param data Data frame
#' @param x_var X variable (unquoted)
#' @param y_var Y variable (unquoted)
#' @param color_var Color variable (unquoted)
#' @param linetype_var Optional linetype variable (unquoted), NULL for none
#' @param facet_formula Facet formula (e.g., q0_lab ~ q1_lab)
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

  # Build aesthetic mapping
  if (is.null(linetype_var)) {
    p <- ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }},
                          color = {{ color_var }}, group = {{ color_var }}))
  } else {
    p <- ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }},
                          color = {{ color_var }}, linetype = {{ linetype_var }},
                          group = interaction({{ color_var }}, {{ linetype_var }})))
  }

  # Add horizontal line if specified
  if (!is.null(hline)) {
    hline_lt <- if (hline_style == "dashed") "dashed" else "dotted"
    hline_col <- if (hline == 0.95) "grey40" else "grey50"
    p <- p + geom_hline(yintercept = hline, linetype = hline_lt, color = hline_col)
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
  if (!is.null(linetype_var)) {
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

  # Add theme
  p <- p + if (use_aggregate_theme) plot_base_theme_aggregate() else plot_base_theme()

  p
}

#' Create coverage plot
plot_coverage <- function(data, facet_formula, title = "Coverage", subtitle = "",
                          linetype_var = NULL, use_aggregate = FALSE,
                          x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = coverage,
    color_var = method,
    linetype_var = {{ linetype_var }},
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Coverage",
    hline = 0.95,
    hline_style = "dashed",
    x_breaks = x_breaks,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Create CI width plot
plot_ci_width <- function(data, facet_formula, title = "CI Width", subtitle = "",
                          linetype_var = NULL, use_aggregate = FALSE,
                          x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = ci_width,
    color_var = method,
    linetype_var = {{ linetype_var }},
    facet_formula = facet_formula,
    title = title,
    subtitle = subtitle,
    y_lab = "Mean CI width",
    x_breaks = x_breaks,
    use_aggregate_theme = use_aggregate,
    point_size = point_size
  )
}

#' Create bias plot
plot_bias <- function(data, facet_formula, title = "Estimator Bias", subtitle = "",
                      linetype_var = NULL, use_aggregate = FALSE,
                      x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = bias,
    color_var = method,
    linetype_var = {{ linetype_var }},
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

#' Create percent bias plot
plot_bias_pct <- function(data, facet_formula, title = "Percent Bias", subtitle = "",
                          linetype_var = NULL, use_aggregate = FALSE,
                          x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = bias_pct,
    color_var = method,
    linetype_var = {{ linetype_var }},
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

#' Create calibration (q) bias plot
plot_q_bias <- function(data, facet_formula, title = "Calibration Estimate Bias", subtitle = "",
                        linetype_var = NULL, use_aggregate = FALSE,
                        x_breaks = NULL, point_size = 1.5) {
  build_line_plot(
    data = data,
    x_var = theta_true,
    y_var = bias,
    color_var = metric,
    linetype_var = {{ linetype_var }},
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
                         {{ linetype_var }}, use_aggregate, x_breaks, point_size)
  save_fn(p_cov, paste0("coverage", filename_suffix), width, height)

  # CI Width
  p_width <- plot_ci_width(ci_data, facet_formula, "CI Width", subtitle,
                           {{ linetype_var }}, use_aggregate, x_breaks, point_size)
  save_fn(p_width, paste0("ci_width", filename_suffix), width, height)

  # Bias
  p_bias <- plot_bias(ci_data, facet_formula, "Estimator Bias", subtitle,
                      {{ linetype_var }}, use_aggregate, x_breaks, point_size)
  save_fn(p_bias, paste0("bias", filename_suffix), width, height)

  # Percent Bias
  p_bias_pct <- plot_bias_pct(ci_data, facet_formula, "Percent Bias", subtitle,
                              {{ linetype_var }}, use_aggregate, x_breaks, point_size)
  save_fn(p_bias_pct, paste0("bias_pct", filename_suffix), width, height)

  # Q Bias
  p_q_bias <- plot_q_bias(q_data, facet_formula, "Calibration Estimate Bias", subtitle,
                          {{ linetype_var }}, use_aggregate, x_breaks, point_size)
  save_fn(p_q_bias, paste0("q_bias", filename_suffix), width, height)
}
