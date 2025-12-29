# Joint misclassification MLE utilities ---------------------------------------

joint_mle_counts <- function(y_cal, yhat_cal, yhat_test) {
  y_cal <- as.integer(y_cal)
  yhat_cal <- as.integer(yhat_cal)
  yhat_test <- as.integer(yhat_test)
  if (length(y_cal) != length(yhat_cal)) {
    stop("y_cal and yhat_cal must have the same length.")
  }
  valid_vals <- c(0L, 1L)
  if (!all(y_cal %in% valid_vals) ||
      !all(yhat_cal %in% valid_vals) ||
      !all(yhat_test %in% valid_vals)) {
    stop("All inputs must be in {0,1}.")
  }
  N11 <- sum(y_cal == 1L & yhat_cal == 1L)
  N10 <- sum(y_cal == 1L & yhat_cal == 0L)
  N01 <- sum(y_cal == 0L & yhat_cal == 1L)
  N00 <- sum(y_cal == 0L & yhat_cal == 0L)
  T1 <- sum(yhat_test == 1L)
  T0 <- length(yhat_test) - T1
  list(
    N11 = N11, N10 = N10,
    N01 = N01, N00 = N00,
    T1 = T1, T0 = T0
  )
}

loglik_misclass <- function(par, counts) {
  a <- par[1]
  b <- par[2]
  c <- par[3]
  theta <- stats::plogis(a)
  q0 <- stats::plogis(b)
  q1 <- stats::plogis(c)
  with(counts, {
    p11 <- theta * q1
    p10 <- theta * (1 - q1)
    p01 <- (1 - theta) * (1 - q0)
    p00 <- (1 - theta) * q0
    p <- theta * q1 + (1 - theta) * (1 - q0)
    J <- q0 + q1 - 1
    if (p11 <= 0 || p10 <= 0 || p01 <= 0 || p00 <= 0 ||
        p <= 0 || p >= 1 || J <= 0) {
      return(-1e20)
    }
    ll_cal <- N11 * log(p11) + N10 * log(p10) +
      N01 * log(p01) + N00 * log(p00)
    ll_test <- T1 * log(p) + T0 * log(1 - p)
    ll_cal + ll_test
  })
}

expected_info_misclass <- function(par, n_test, m_cal) {
  theta <- par[1]
  q0 <- par[2]
  q1 <- par[3]
  if (theta <= 0 || theta >= 1 ||
      q0 <= 0 || q0 >= 1 ||
      q1 <= 0 || q1 >= 1) {
    stop("Parameters must be in (0,1).")
  }
  pi11 <- theta * q1
  pi10 <- theta * (1 - q1)
  pi01 <- (1 - theta) * (1 - q0)
  pi00 <- (1 - theta) * q0
  s11 <- c(1 / theta, 0, 1 / q1)
  s10 <- c(1 / theta, 0, -1 / (1 - q1))
  s01 <- c(-1 / (1 - theta), -1 / (1 - q0), 0)
  s00 <- c(-1 / (1 - theta), 1 / q0, 0)
  I_cal_one <- matrix(0, 3, 3)
  I_cal_one <- I_cal_one + pi11 * (s11 %*% t(s11))
  I_cal_one <- I_cal_one + pi10 * (s10 %*% t(s10))
  I_cal_one <- I_cal_one + pi01 * (s01 %*% t(s01))
  I_cal_one <- I_cal_one + pi00 * (s00 %*% t(s00))
  p <- theta * q1 + (1 - theta) * (1 - q0)
  J <- q0 + q1 - 1
  if (p <= 0 || p >= 1 || J <= 0) {
    stop("Invalid parameters: p must be in (0,1) and J>0.")
  }
  s1 <- c(J / p, -(1 - theta) / p, theta / p)
  s0 <- c(-J / (1 - p), (1 - theta) / (1 - p), -theta / (1 - p))
  I_test_one <- matrix(0, 3, 3)
  I_test_one <- I_test_one + p * (s1 %*% t(s1))
  I_test_one <- I_test_one + (1 - p) * (s0 %*% t(s0))
  I_cal <- m_cal * I_cal_one
  I_test <- n_test * I_test_one
  I_tot <- I_cal + I_test
  colnames(I_tot) <- rownames(I_tot) <- c("theta", "q0", "q1")
  I_tot
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

logit_safe <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)

wald_ci_prob <- function(est, se, level = 0.90, clip = TRUE) {
  if (is.na(se)) {
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- est + c(-1, 1) * z * se
  if (clip) {
    ci <- clip01(ci)
  }
  ci
}

#' Logit-transformed Wald CI for probability parameters.
#'
#' Constructs a confidence interval on the logit scale and transforms back,
#' which provides better coverage for probabilities near 0 or 1.
#'
#' @param est Point estimate on probability scale.
#' @param se Standard error of the estimate.
#' @param level Confidence level.
#'
#' @return Numeric vector of length 2 with (lower, upper) bounds.
#' @export
logit_ci_prob <- function(est, se, level = 0.90) {
  if (is.na(se) || is.na(est)) {
    return(c(NA_real_, NA_real_))
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  # Clamp estimate away from boundaries for logit transform
  th <- logit_safe(est)
  # Delta method: Var(logit(theta)) ≈ Var(theta) / (theta * (1-theta))^2
  gprime <- 1 / (th * (1 - th))
  se_logit <- gprime * se
  # CI on logit scale
  logit_lo <- stats::qlogis(th) - z * se_logit
  logit_hi <- stats::qlogis(th) + z * se_logit
  # Transform back
  c(stats::plogis(logit_lo), stats::plogis(logit_hi))
}

#' Joint MLE under the binary misclassification model.
#'
#' Estimates \\theta, q0, and q1 by maximizing the joint likelihood built from
#' labeled calibration data and unlabeled surrogate outputs, then obtains
#' standard errors via the inverse Fisher information (observed and expected).
#'
#' @param y_cal Human labels on the calibration split (0/1).
#' @param yhat_cal Surrogate predictions on the calibration split (0/1).
#' @param yhat_test Surrogate predictions on the unlabeled/test split (0/1).
#' @param level Confidence level for the Wald intervals. Defaults to 0.90.
#'
#' @return A list containing the MLEs, variance-covariance matrices from the
#'   observed and expected Fisher information, and Wald intervals for \\theta.
#' @export
fit_misclass_mle <- function(y_cal, yhat_cal, yhat_test, level = 0.90) {
  counts <- joint_mle_counts(y_cal, yhat_cal, yhat_test)
  n_test <- counts$T1 + counts$T0
  m_cal <- counts$N11 + counts$N10 + counts$N01 + counts$N00
  m1 <- counts$N11 + counts$N10
  m0 <- counts$N01 + counts$N00
  m <- m0 + m1
  theta0 <- if (m > 0) m1 / m else 0.5
  q10 <- if (m1 > 0) counts$N11 / m1 else 0.9
  q00 <- if (m0 > 0) counts$N00 / m0 else 0.9
  p_hat <- if (n_test > 0) counts$T1 / n_test else theta0
  J_init <- q00 + q10 - 1
  theta_mom <- if (J_init > 1e-4) clip01((p_hat + q00 - 1) / J_init) else theta0
  theta_init <- if (is.finite(theta_mom)) theta_mom else theta0
  par0 <- stats::qlogis(logit_safe(c(theta_init, q00, q10)))
  negloglik <- function(par) -loglik_misclass(par, counts)
  opt <- stats::optim(
    par = par0,
    fn = negloglik,
    method = "BFGS",
    control = list(maxit = 2000, reltol = 1e-10)
  )
  if (opt$convergence != 0) {
    warning("Optimization may not have converged: code = ", opt$convergence)
  }
  par_hat <- opt$par
  a_hat <- par_hat[1]
  b_hat <- par_hat[2]
  c_hat <- par_hat[3]
  theta_hat <- stats::plogis(a_hat)
  q0_hat <- stats::plogis(b_hat)
  q1_hat <- stats::plogis(c_hat)
  dtheta_da <- theta_hat * (1 - theta_hat)
  dq0_db <- q0_hat * (1 - q0_hat)
  dq1_dc <- q1_hat * (1 - q1_hat)
  J_jac <- diag(c(dtheta_da, dq0_db, dq1_dc))
  rownames(J_jac) <- c("theta", "q0", "q1")
  colnames(J_jac) <- c("a", "b", "c")
  H_abc <- numDeriv::hessian(
    func = loglik_misclass,
    x = par_hat,
    counts = counts
  )
  info_abc_obs <- -H_abc
  vcov_abc_obs <- tryCatch(
    solve(info_abc_obs),
    error = function(e) {
      warning("Observed information nearly singular; returning NA vcov.")
      matrix(NA_real_, nrow = 3, ncol = 3)
    }
  )
  vcov_tq0q1_obs <- J_jac %*% vcov_abc_obs %*% t(J_jac)
  colnames(vcov_tq0q1_obs) <- rownames(vcov_tq0q1_obs) <- c("theta", "q0", "q1")
  exp_info <- tryCatch(
    expected_info_misclass(
      par = c(theta_hat, q0_hat, q1_hat),
      n_test = n_test,
      m_cal = m_cal
    ),
    error = function(e) {
      warning("Expected information unavailable: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(exp_info)) {
    vcov_tq0q1_exp <- matrix(NA_real_, nrow = 3, ncol = 3)
    colnames(vcov_tq0q1_exp) <- rownames(vcov_tq0q1_exp) <- c("theta", "q0", "q1")
  } else {
    vcov_tq0q1_exp <- tryCatch(
      solve(exp_info),
      error = function(e) {
        warning("Expected information nearly singular; returning NA vcov.")
        matrix(NA_real_, nrow = 3, ncol = 3)
      }
    )
  }
  se_obs <- sqrt(diag(vcov_tq0q1_obs))
  se_exp <- sqrt(diag(vcov_tq0q1_exp))
  # Use logit-transformed CI for better coverage near boundaries
  ci_theta_obs <- logit_ci_prob(theta_hat, se_obs["theta"], level = level)
  ci_theta_exp <- logit_ci_prob(theta_hat, se_exp["theta"], level = level)
  list(
    theta_hat = theta_hat,
    q0_hat = q0_hat,
    q1_hat = q1_hat,
    par_hat = c(theta = theta_hat, q0 = q0_hat, q1 = q1_hat),
    vcov_obs = vcov_tq0q1_obs,
    vcov_exp = vcov_tq0q1_exp,
    se_obs = se_obs,
    se_exp = se_exp,
    ci_theta_obs = ci_theta_obs,
    ci_theta_exp = ci_theta_exp,
    logLik = loglik_misclass(par_hat, counts),
    counts = counts,
    opt = opt
  )
}

#' Grab the theta-specific summary from a joint MLE fit.
#'
#' @param fit Output from [fit_misclass_mle()].
#' @param level Confidence level for the Wald interval.
#'
#' @return List with theta_hat, observed/expected SEs, and CIs.
#' @export
grab_theta <- function(fit, level = 0.90) {
  list(
    theta_hat = fit$theta_hat,
    se_obs = unname(fit$se_obs["theta"]),
    se_exp = unname(fit$se_exp["theta"]),
    ci_obs = logit_ci_prob(fit$theta_hat, fit$se_obs["theta"], level = level),
    ci_exp = logit_ci_prob(fit$theta_hat, fit$se_exp["theta"], level = level)
  )
}

summarize_parameter <- function(name, estimates, se_obs, se_exp, truth) {
  data.frame(
    parameter = name,
    true_value = truth,
    mean_hat = mean(estimates),
    bias = mean(estimates) - truth,
    empirical_sd = stats::sd(estimates),
    mean_se_obs = mean(se_obs, na.rm = TRUE),
    mean_se_exp = mean(se_exp, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_sim_plots <- function(sim_data, true_vals, coverage_tbl, level) {
  hist_plot <- function(values, truth, label) {
    df <- data.frame(est = values)
    ggplot2::ggplot(df, ggplot2::aes(x = est)) +
      ggplot2::geom_histogram(bins = 30, color = "white", fill = "#3182bd", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = truth, color = "#d62728", linetype = "dashed", linewidth = 0.8) +
      ggplot2::labs(title = paste0(label, " MLEs"), x = paste0(label, " hat"), y = "Count") +
      ggplot2::theme_minimal()
  }
  list(
    theta_hist = hist_plot(sim_data$theta_hat, true_vals["theta"], "theta"),
    q0_hist = hist_plot(sim_data$q0_hat, true_vals["q0"], "q0"),
    q1_hist = hist_plot(sim_data$q1_hat, true_vals["q1"], "q1"),
    coverage = ggplot2::ggplot(coverage_tbl, ggplot2::aes(x = parameter, y = coverage, fill = method)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.6), width = 0.5, alpha = 0.8) +
      ggplot2::geom_hline(yintercept = level, linetype = "dotted", color = "#636363") +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(y = "Coverage", x = NULL, title = "Empirical coverage by parameter") +
      ggplot2::theme_minimal()
  )
}

#' Monte Carlo diagnostics for the joint MLE.
#'
#' Simulates repeated calibration/test splits under the binary misclassification
#' model to study estimator bias, standard errors, and interval coverage.
#'
#' @param B Number of Monte Carlo repetitions.
#' @param N Total sample size to simulate.
#' @param theta_true,q0_true,q1_true Ground-truth prevalence/specificity/sensitivity.
#' @param m_cal Size of the calibration set.
#' @param level Confidence level for coverage diagnostics.
#'
#' @return List containing aggregate summaries, coverage, raw simulation draws,
#'   and ggplot2 objects for quick inspection.
#' @export
sim_compare_info <- function(B = 1000, N = 5000,
                             theta_true = 0.3,
                             q0_true = 0.8,
                             q1_true = 0.85,
                             m_cal = 500,
                             level = 0.90) {
  if (m_cal >= N) {
    stop("m_cal must be less than N.")
  }
  sim_data <- vector("list", B)
  true_vals <- c(theta = theta_true, q0 = q0_true, q1 = q1_true)
  for (b in seq_len(B)) {
    Y <- stats::rbinom(N, 1, theta_true)
    Yhat <- ifelse(
      Y == 1,
      stats::rbinom(N, 1, q1_true),
      stats::rbinom(N, 1, 1 - q0_true)
    )
    idx_cal <- sample(seq_len(N), m_cal)
    idx_test <- setdiff(seq_len(N), idx_cal)
    y_cal <- Y[idx_cal]
    yhat_cal <- Yhat[idx_cal]
    yhat_test <- Yhat[idx_test]
    fit <- fit_misclass_mle(y_cal, yhat_cal, yhat_test, level = level)
    ci_theta_obs <- logit_ci_prob(fit$theta_hat, fit$se_obs["theta"], level)
    ci_theta_exp <- logit_ci_prob(fit$theta_hat, fit$se_exp["theta"], level)
    ci_q0_obs <- logit_ci_prob(fit$q0_hat, fit$se_obs["q0"], level)
    ci_q0_exp <- logit_ci_prob(fit$q0_hat, fit$se_exp["q0"], level)
    ci_q1_obs <- logit_ci_prob(fit$q1_hat, fit$se_obs["q1"], level)
    ci_q1_exp <- logit_ci_prob(fit$q1_hat, fit$se_exp["q1"], level)
    sim_data[[b]] <- data.frame(
      iter = b,
      theta_hat = fit$theta_hat,
      q0_hat = fit$q0_hat,
      q1_hat = fit$q1_hat,
      se_theta_obs = fit$se_obs["theta"],
      se_theta_exp = fit$se_exp["theta"],
      se_q0_obs = fit$se_obs["q0"],
      se_q0_exp = fit$se_exp["q0"],
      se_q1_obs = fit$se_obs["q1"],
      se_q1_exp = fit$se_exp["q1"],
      theta_lwr_obs = ci_theta_obs[1],
      theta_upr_obs = ci_theta_obs[2],
      theta_lwr_exp = ci_theta_exp[1],
      theta_upr_exp = ci_theta_exp[2],
      q0_lwr_obs = ci_q0_obs[1],
      q0_upr_obs = ci_q0_obs[2],
      q0_lwr_exp = ci_q0_exp[1],
      q0_upr_exp = ci_q0_exp[2],
      q1_lwr_obs = ci_q1_obs[1],
      q1_upr_obs = ci_q1_obs[2],
      q1_lwr_exp = ci_q1_exp[1],
      q1_upr_exp = ci_q1_exp[2]
    )
  }
  sim_df <- do.call(rbind, sim_data)
  summary_tbl <- rbind(
    summarize_parameter("theta", sim_df$theta_hat, sim_df$se_theta_obs, sim_df$se_theta_exp, theta_true),
    summarize_parameter("q0", sim_df$q0_hat, sim_df$se_q0_obs, sim_df$se_q0_exp, q0_true),
    summarize_parameter("q1", sim_df$q1_hat, sim_df$se_q1_obs, sim_df$se_q1_exp, q1_true)
  )
  coverage_tbl <- data.frame(
    parameter = rep(c("theta", "q0", "q1"), each = 2),
    method = rep(c("observed", "expected"), times = 3),
    coverage = c(
      mean(sim_df$theta_lwr_obs <= theta_true & theta_true <= sim_df$theta_upr_obs, na.rm = TRUE),
      mean(sim_df$theta_lwr_exp <= theta_true & theta_true <= sim_df$theta_upr_exp, na.rm = TRUE),
      mean(sim_df$q0_lwr_obs <= q0_true & q0_true <= sim_df$q0_upr_obs, na.rm = TRUE),
      mean(sim_df$q0_lwr_exp <= q0_true & q0_true <= sim_df$q0_upr_exp, na.rm = TRUE),
      mean(sim_df$q1_lwr_obs <= q1_true & q1_true <= sim_df$q1_upr_obs, na.rm = TRUE),
      mean(sim_df$q1_lwr_exp <= q1_true & q1_true <= sim_df$q1_upr_exp, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  coverage_tbl$method <- factor(coverage_tbl$method, levels = c("observed", "expected"))
  plots <- make_sim_plots(sim_df, true_vals, coverage_tbl, level)
  list(
    summary = summary_tbl,
    coverage = coverage_tbl,
    theta_summary = summary_tbl[summary_tbl$parameter == "theta", , drop = FALSE],
    sim_data = sim_df,
    plots = plots
  )
}
