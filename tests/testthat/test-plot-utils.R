if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
}
if (requireNamespace("grid", quietly = TRUE)) {
  library(grid)
}

test_that("build_line_plot returns a ggplot object without linetype", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    theta_true = rep(c(0.1, 0.2, 0.3), 2),
    coverage = runif(6, 0.8, 1.0),
    method = rep(c("Naive", "PPI"), each = 3),
    q0_lab = "q0=0.8",
    q1_lab = "q1=0.8",
    stringsAsFactors = FALSE
  )

  p <- build_line_plot(
    data = df,
    x_var = theta_true,
    y_var = coverage,
    color_var = method,
    linetype_var = NULL,
    facet_formula = q0_lab ~ q1_lab,
    title = "Coverage",
    y_lab = "Coverage",
    hline = 0.95,
    hline_style = "dashed"
  )

  expect_s3_class(p, "ggplot")
})

test_that("build_line_plot returns a ggplot object with linetype string column", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    theta_true = rep(c(0.1, 0.2, 0.3), 4),
    bias = rnorm(12, 0, 0.01),
    method = rep(c("Naive", "PPI"), each = 6),
    split_label = rep(c("20:80", "50:50"), times = 6),
    q_val_lab = "q0=q1=0.8",
    label_ratio_lab = "10% labeled",
    stringsAsFactors = FALSE
  )

  p <- build_line_plot(
    data = df,
    x_var = theta_true,
    y_var = bias,
    color_var = method,
    linetype_var = "split_label",
    facet_formula = q_val_lab ~ label_ratio_lab,
    title = "Bias",
    y_lab = "Bias",
    hline = 0,
    hline_style = "dotted",
    use_aggregate_theme = TRUE
  )

  expect_s3_class(p, "ggplot")
})

test_that("plot_coverage wrapper works and returns ggplot", {
  skip_if_not_installed("ggplot2")

  df <- data.frame(
    theta_true = rep(c(0.1, 0.2, 0.3), 2),
    coverage = runif(6, 0.8, 1.0),
    method = rep(c("Naive", "PPI"), each = 3),
    q0_lab = "q0=0.8",
    q1_lab = "q1=0.8",
    stringsAsFactors = FALSE
  )

  p <- plot_coverage(df, q0_lab ~ q1_lab, x_breaks = c(0.1, 0.2, 0.3))
  expect_s3_class(p, "ggplot")
})
