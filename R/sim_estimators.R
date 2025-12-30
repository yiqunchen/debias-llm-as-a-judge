# Core estimators for LLM-as-a-judge debiasing

clamp01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

#' Generate Bernoulli responses with controllable prevalence.
#'
#' Draws i.i.d. binary outcomes with the specified prevalence. The simpler
#' exchangeable design is appropriate when estimators do not condition on
#' covariates.
#'
#' @param N Integer; number of observations to simulate.
#' @param theta Target marginal prevalence for the Bernoulli outcome.
#' @param signal Unused; retained for API compatibility.
#'
#' @return A list containing outcome vector `Y`.
#' @keywords internal
#' @export
generate_dgp_data <- function(N, theta, signal = 1) {
  Y <- rbinom(N, 1, theta)
  list(Y = Y)
}

#' Wald-style confidence interval via logit transform.
#'
#' @param theta_hat Point estimate on probability scale.
#' @param var_hat Estimated variance of `theta_hat`.
#' @param alpha Miscoverage level for the interval.
#'
#' @return A list with `lower` and `upper` bounds.
#' @export
logit_ci <- function(theta_hat, var_hat, alpha = 0.10) {
  th <- clamp01(theta_hat)
  z <- qnorm(1 - alpha / 2)
  se <- sqrt(max(var_hat, 0))
  gprime <- 1 / (th * (1 - th))
  se_logit <- gprime * se
  lo <- plogis(qlogis(th) - z * se_logit)
  hi <- plogis(qlogis(th) + z * se_logit)
  list(lower = lo, upper = hi)
}

#' Prediction-powered mean estimator with logit CI.
#'
#' @param Y_L Vector of human labels on the calibration set.
#' @param f_L Surrogate predictions on the calibration set.
#' @param f_U Surrogate predictions on the large unlabeled/test set.
#' @param alpha Miscoverage level for the confidence interval.
#'
#' @return List with `theta`, `var`, and CI endpoints.
#' @export
ppi_point_and_ci <- function(Y_L, f_L, f_U, alpha = 0.10) {
  n <- length(Y_L)
  N <- length(f_U)
  theta_hat <- mean(f_U) + mean(Y_L - f_L)
  eps_L <- Y_L - f_L
  var_eps <- if (n > 1) var(eps_L) else 0
  var_fU <- if (N > 1) var(f_U) else 0
  var_hat <- var_eps / n + var_fU / N
  ci <- logit_ci(theta_hat, var_hat, alpha)
  list(theta = theta_hat,
       var = var_hat,
       ci_lower = ci$lower,
       ci_upper = ci$upper)
}

#' LLM-as-a-judge corrected estimator with finite-sample adjustment.
#'
#' @param p_hat Average LLM score on the test set.
#' @param q0_hat Calibration specificity estimate.
#' @param q1_hat Calibration sensitivity estimate.
#' @param n Test set size.
#' @param m0,m1 Calibration sample sizes for negatives/positives.
#' @param alpha Miscoverage level.
#' @param eps_J Small threshold to guard against near-random judges.
#'
#' @return List with `theta`, `var`, and CI endpoints.
#' @export
llm_point_and_ci <- function(p_hat, q0_hat, q1_hat, n, m0, m1,
                             alpha = 0.10, eps_J = 1e-4) {
  z <- qnorm(1 - alpha / 2)
  p_tilde <- (n * p_hat + z^2 / 2) / (n + z^2)
  q0_tilde <- (m0 * q0_hat + 1) / (m0 + 2)
  q1_tilde <- (m1 * q1_hat + 1) / (m1 + 2)
  n_tilde <- n + z^2
  m0_tilde <- m0 + 2
  m1_tilde <- m1 + 2
  J_tilde <- q0_tilde + q1_tilde - 1
  if (is.na(J_tilde) || abs(J_tilde) < eps_J) {
    return(list(theta = NA_real_, var = NA_real_,
                ci_lower = NA_real_, ci_upper = NA_real_))
  }
  theta_tilde <- (p_tilde + q0_tilde - 1) / J_tilde
  theta_tilde <- pmin(pmax(theta_tilde, 0), 1)
  var_hat <- (
    p_tilde * (1 - p_tilde) / n_tilde +
      (1 - theta_tilde)^2 * q0_tilde * (1 - q0_tilde) / m0_tilde +
      theta_tilde^2 * q1_tilde * (1 - q1_tilde) / m1_tilde
  ) / (J_tilde^2)
  var_hat <- max(var_hat, 0)
  se_hat <- sqrt(var_hat)
  th <- clamp01(theta_tilde)
  gprime <- 1 / (th * (1 - th))
  se_logit <- gprime * se_hat
  logit_lo <- qlogis(th) - z * se_logit
  logit_hi <- qlogis(th) + z * se_logit
  ci_lower <- plogis(logit_lo)
  ci_upper <- plogis(logit_hi)
  list(theta = theta_tilde,
       var = var_hat,
       ci_lower = ci_lower,
       ci_upper = ci_upper)
}

ppi_pp_var_hat <- function(lambda, Y_L, f_L, f_U) {
  n <- length(Y_L)
  N <- length(f_U)
  Z_L <- Y_L - lambda * f_L
  var_Z_L <- if (n > 1) var(Z_L) else 0
  var_f_U <- if (N > 1) var(f_U) else 0
  var_Z_L / n + lambda^2 * var_f_U / N
}

lambda_hat_ppi_pp_numeric <- function(Y_L, f_L, f_U,
                                      lambda_range = c(0, 1)) {
  optimize(
    f = ppi_pp_var_hat,
    interval = lambda_range,
    Y_L = Y_L,
    f_L = f_L,
    f_U = f_U
  )
}

#' Tuned PPI++ estimator using numeric lambda search.
#'
#' @param Y_L,f_L,f_U Observed labels and surrogate scores.
#' @param alpha Miscoverage level.
#' @param lambda_range Interval to constrain lambda.
#'
#' @return List containing `theta`, `lambda`, variance, and CI limits.
#' @export
ppi_pp_point_and_ci_general <- function(Y_L, f_L, f_U,
                                        alpha = 0.10,
                                        lambda_range = c(0, 1)) {
  lam_fit <- lambda_hat_ppi_pp_numeric(Y_L, f_L, f_U, lambda_range)
  lambda_hat <- lam_fit$minimum
  theta_hat <- mean(Y_L) + lambda_hat * (mean(f_U) - mean(f_L))
  var_hat <- ppi_pp_var_hat(lambda_hat, Y_L, f_L, f_U)
  var_hat <- max(var_hat, 0)
  theta_tilde <- clamp01(theta_hat)
  var_logit <- var_hat / (theta_tilde^2 * (1 - theta_tilde)^2)
  se_logit <- sqrt(max(var_logit, 0))
  z <- qnorm(1 - alpha / 2)
  logit_lo <- qlogis(theta_tilde) - z * se_logit
  logit_hi <- qlogis(theta_tilde) + z * se_logit
  ci_lower <- plogis(logit_lo)
  ci_upper <- plogis(logit_hi)
  theta_hat <- pmin(pmax(theta_hat, 0), 1)
  list(theta = theta_hat,
       lambda = lambda_hat,
       var = var_hat,
       ci_lower = ci_lower,
       ci_upper = ci_upper)
}

#' EIF (Efficient Influence Function) estimator for binary surrogates.
#'
#' Implements the post-stratified efficient estimator from the EIF theory.
#' For binary surrogates, this estimates \eqn{\theta = E[Y]} using:
#' \deqn{\hat\theta_{EIF} = (1-\hat p)\hat\mu_0 + \hat p \hat\mu_1}
#' where \eqn{\hat\mu_0, \hat\mu_1} are the conditional means of Y given
#' \eqn{\hat Y = 0, 1} estimated from calibration data, and \eqn{\hat p}
#' is the proportion of \eqn{\hat Y = 1} in the test sample.
#'
#' This estimator is asymptotically efficient in Model A (surrogates are
#' conditionally independent of outcomes given covariates).
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions (0/1) on the calibration set.
#' @param Yhat_test Vector of surrogate predictions (0/1) on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#'
#' @return A list containing:
#'   \describe{
#'     \item{theta}{Point estimate of the prevalence.}
#'     \item{var}{Estimated variance of the estimator.}
#'     \item{ci_lower}{Lower bound of the confidence interval.}
#'     \item{ci_upper}{Upper bound of the confidence interval.}
#'     \item{mu0_hat}{Estimated \eqn{E[Y | \hat Y = 0]}.}
#'     \item{mu1_hat}{Estimated \eqn{E[Y | \hat Y = 1]}.}
#'     \item{p_hat}{Proportion of \eqn{\hat Y = 1} in test set.}
#'   }
#' @export
eif_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  # Return NA results if inputs are problematic
  na_result <- list(
    theta = NA_real_,
    var = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_,
    mu0_hat = NA_real_,
    mu1_hat = NA_real_,
    p_hat = NA_real_
  )

  # Check for sufficient valid data
  if (length(Y_cal) != length(Yhat_cal) ||
      length(Y_cal) == 0 ||
      length(Yhat_test) == 0) {
    return(na_result)
  }

  # Filter to valid calibration observations
  valid_cal <- !is.na(Y_cal) & !is.na(Yhat_cal) & Yhat_cal %in% c(0, 1)
  Y_cal <- Y_cal[valid_cal]
  Yhat_cal <- Yhat_cal[valid_cal]

  # Filter to valid test observations
  valid_test <- !is.na(Yhat_test) & Yhat_test %in% c(0, 1)
  Yhat_test <- Yhat_test[valid_test]

  # Check we still have enough data
  if (length(Y_cal) == 0 || length(Yhat_test) == 0) {
    return(na_result)
  }

  # Calibration sample sizes
  idx0 <- which(Yhat_cal == 0)
  idx1 <- which(Yhat_cal == 1)
  m0 <- length(idx0)
  m1 <- length(idx1)
  m <- m0 + m1

  # Test sample size
  n <- length(Yhat_test)

  # Handle edge cases
  if (m0 == 0 || m1 == 0 || n == 0) {
    return(list(
      theta = NA_real_,
      var = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      mu0_hat = NA_real_,
      mu1_hat = NA_real_,
      p_hat = NA_real_
    ))
  }

  # Estimate mu0 and mu1 from calibration data
  mu0_hat <- mean(Y_cal[idx0])
  mu1_hat <- mean(Y_cal[idx1])

  # Estimate p from test data
  p_hat <- mean(Yhat_test)

  # Point estimate: post-stratified estimator
  theta_hat <- (1 - p_hat) * mu0_hat + p_hat * mu1_hat

  # Variance estimation using delta method
  var_p <- p_hat * (1 - p_hat) / n
  var_mu0 <- mu0_hat * (1 - mu0_hat) / m0
  var_mu1 <- mu1_hat * (1 - mu1_hat) / m1

  var_hat <- (mu1_hat - mu0_hat)^2 * var_p +
    (1 - p_hat)^2 * var_mu0 +
    p_hat^2 * var_mu1

  var_hat <- max(var_hat, 0)

  # Confidence interval using logit transform for better coverage
  ci <- logit_ci(theta_hat, var_hat, alpha)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = ci$lower,
    ci_upper = ci$upper,
    mu0_hat = mu0_hat,
    mu1_hat = mu1_hat,
    p_hat = p_hat
  )
}

