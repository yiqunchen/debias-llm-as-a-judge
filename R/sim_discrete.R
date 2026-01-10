# Discrete Prediction Simulation: Core Functions
#
# DGP: Z ~ Unif({1, 2, 3}), Y | Z ~ N(mu_Z, sigma^2), Yhat = Z
#
# Key insight: Yhat is the mixture component label (discrete: 1, 2, 3)
# Even though Yhat perfectly identifies which component Y came from,
# naive mean(Yhat) is biased for E[Y] because Yhat is on wrong scale!
#
# Example: mu = c(1, 5, 9), sigma = 1
#   E[Y] = (1 + 5 + 9) / 3 = 5
#   E[Yhat] = (1 + 2 + 3) / 3 = 2
#   Naive bias = 2 - 5 = -3

# ============================================================================
# 1. DATA GENERATING PROCESS
# ============================================================================

#' Generate discrete prediction data
#'
#' \eqn{Z \sim \mathrm{Unif}(\{1,\ldots,K\}), Y \mid Z \sim N(\mu_Z, \sigma^2), \hat{Y} = Z}
#'
#' @param N Sample size
#' @param mu Vector of means for each mixture component
#' @param sigma Standard deviation (same for all components)
#' @param seed Random seed
#' @return List with Z, Y, Yhat
#' @export
generate_discrete_dgp <- function(N, mu = c(1, 5, 9), sigma = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  K <- length(mu)

  # Z ~ Unif({1, ..., K})
  Z <- sample(1:K, N, replace = TRUE)

  # Y | Z ~ N(mu[Z], sigma^2)
  Y <- rnorm(N, mean = mu[Z], sd = sigma)

  # Yhat = Z (discrete label prediction)
  Yhat <- Z

  list(
    Z = Z,
    Y = Y,
    Yhat = Yhat,
    mu = mu,
    sigma = sigma,
    K = K,
    true_mean = mean(mu)  # E[Y] = mean of component means (uniform mixture)
  )
}

# ============================================================================
# 2. ESTIMATORS
# ============================================================================

#' Naive estimator: just use mean of discrete predictions
#'
#' @param Yhat_test Vector of discrete surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
naive_discrete_point_and_ci <- function(Yhat_test, alpha = 0.10) {
  n <- length(Yhat_test)
  theta_hat <- mean(Yhat_test)
  var_hat <- var(Yhat_test) / n

  z <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z * se,
    ci_upper = theta_hat + z * se
  )
}

#' PPI estimator for discrete predictions
#'
#' \eqn{\hat{\theta} = \mathrm{mean}(\hat{Y}_{test}) + (\mathrm{mean}(Y_{cal}) - \mathrm{mean}(\hat{Y}_{cal}))}
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions on the calibration set.
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
ppi_discrete_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  # Bias correction
  bias_hat <- mean(Y_cal) - mean(Yhat_cal)
  theta_hat <- mean(Yhat_test) + bias_hat

  # Variance
  eps_cal <- Y_cal - Yhat_cal
  var_eps <- var(eps_cal)
  var_test <- var(Yhat_test)
  var_hat <- var_eps / m + var_test / n

  z <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z * se,
    ci_upper = theta_hat + z * se
  )
}

#' PPI++ estimator for discrete predictions
#'
#' Uses closed-form optimal lambda (no constraints):
#' \eqn{\lambda^* = \frac{\mathrm{Cov}(Y,\hat{Y})}{\mathrm{Var}(\hat{Y}_{cal}) + (m/n)\mathrm{Var}(\hat{Y}_{test})}}
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions on the calibration set.
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
ppi_pp_discrete_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  # Closed-form optimal lambda
  cov_hat <- cov(Y_cal, Yhat_cal)
  var_cal <- var(Yhat_cal)
  var_test <- var(Yhat_test)

  lambda_hat <- cov_hat / (var_cal + (m / n) * var_test)

  # Point estimate
  theta_hat <- mean(Y_cal) + lambda_hat * (mean(Yhat_test) - mean(Yhat_cal))

  # Variance at optimal lambda
  # V(lambda) = Var(Y - lambda*Yhat) / m + lambda^2 * Var(Yhat_test) / n
  var_resid <- var(Y_cal - lambda_hat * Yhat_cal)
  var_hat <- var_resid / m + lambda_hat^2 * var_test / n
  var_hat <- max(var_hat, 1e-10)

  z <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    lambda = lambda_hat,
    ci_lower = theta_hat - z * se,
    ci_upper = theta_hat + z * se
  )
}

#' EIF estimator with per-category calibration (optimal for discrete)
#'
#' Learns \eqn{g: \{1,\ldots,K\} \to \mathbb{R}} such that
#' \eqn{g(z) = E[Y \mid \hat{Y} = z]}.
#' For discrete Yhat, this is just the conditional mean per category.
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions on the calibration set.
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
eif_discrete_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  # Learn g(z) = E[Y | Yhat = z] for each discrete value
  categories <- sort(unique(c(Yhat_cal, Yhat_test)))

  g_map <- sapply(categories, function(z) {
    idx <- Yhat_cal == z
    if (sum(idx) > 0) {
      mean(Y_cal[idx])
    } else {
      mean(Y_cal)  # fallback if no observations
    }
  })
  names(g_map) <- as.character(categories)

  # Apply g to calibration and test
  g_cal <- g_map[as.character(Yhat_cal)]
  g_test <- g_map[as.character(Yhat_test)]

  # Estimate
  theta_hat <- mean(g_test)

  # Variance: var(g(Yhat_test))/n + var(Y - g(Yhat))/m
  var_g_test <- var(g_test) / n
  resid_cal <- Y_cal - g_cal
  var_cal <- var(resid_cal) / m

  var_hat <- var_g_test + var_cal
  var_hat <- max(var_hat, 1e-10)

  z_crit <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z_crit * se,
    ci_upper = theta_hat + z_crit * se,
    g_map = g_map
  )
}

#' EIF estimator with linear calibration (suboptimal for discrete)
#'
#' Fits g(Yhat) = a + b*Yhat using linear regression.
#' This is suboptimal when Yhat is discrete but shows importance of proper calibration.
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions on the calibration set.
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
eif_linear_discrete_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  # Fit linear calibration: Y ~ Yhat
  fit <- lm(Y_cal ~ Yhat_cal)

  # Apply g to calibration and test
  g_cal <- predict(fit, newdata = data.frame(Yhat_cal = Yhat_cal))
  g_test <- predict(fit, newdata = data.frame(Yhat_cal = Yhat_test))

  # Estimate
  theta_hat <- mean(g_test)

  # Variance
  var_g_test <- var(g_test) / n
  resid_cal <- Y_cal - g_cal
  var_cal <- var(resid_cal) / m

  var_hat <- var_g_test + var_cal
  var_hat <- max(var_hat, 1e-10)

  z_crit <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z_crit * se,
    ci_upper = theta_hat + z_crit * se,
    coef = coef(fit)
  )
}

#' EIF estimator with GAM calibration
#'
#' Fits g(Yhat) using a GAM with smooth term.
#' For discrete Yhat with few categories, this behaves similarly to per-category.
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions on the calibration set.
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
eif_gam_discrete_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  df_cal <- data.frame(Y = Y_cal, Yhat = Yhat_cal)
  df_test <- data.frame(Yhat = Yhat_test)

  # Fit GAM - use k = number of unique values (at most 3 for our case)
  k_val <- min(length(unique(Yhat_cal)), 10)

  fit <- tryCatch(
    mgcv::gam(Y ~ s(Yhat, k = k_val), data = df_cal),
    error = function(e) lm(Y ~ Yhat, data = df_cal)
  )

  # Apply g to calibration and test
  g_cal <- predict(fit, newdata = df_cal)
  g_test <- predict(fit, newdata = df_test)

  # Estimate
  theta_hat <- mean(g_test)

  # Variance
  var_g_test <- var(g_test) / n
  resid_cal <- Y_cal - g_cal
  var_cal <- var(resid_cal) / m

  var_hat <- var_g_test + var_cal
  var_hat <- max(var_hat, 1e-10)

  z_crit <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z_crit * se,
    ci_upper = theta_hat + z_crit * se
  )
}

#' EIF estimator with spline calibration
#'
#' Fits g(Yhat) using natural splines.
#' For discrete Yhat with few categories, this behaves similarly to per-category.
#'
#' @param Y_cal Vector of true labels on the calibration set.
#' @param Yhat_cal Vector of surrogate predictions on the calibration set.
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
eif_spline_discrete_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  # For 3 discrete values, use df = 2 (linear interpolation through 3 points)
  df_spline <- min(length(unique(Yhat_cal)) - 1, 4)
  df_spline <- max(df_spline, 1)

  fit <- tryCatch(
    lm(Y_cal ~ splines::ns(Yhat_cal, df = df_spline)),
    error = function(e) lm(Y_cal ~ Yhat_cal)
  )

  # Apply g to calibration and test
  g_cal <- predict(fit, newdata = data.frame(Yhat_cal = Yhat_cal))
  g_test <- predict(fit, newdata = data.frame(Yhat_cal = Yhat_test))

  # Estimate
  theta_hat <- mean(g_test)

  # Variance
  var_g_test <- var(g_test) / n
  resid_cal <- Y_cal - g_cal
  var_cal <- var(resid_cal) / m

  var_hat <- var_g_test + var_cal
  var_hat <- max(var_hat, 1e-10)

  z_crit <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z_crit * se,
    ci_upper = theta_hat + z_crit * se
  )
}

#' Oracle estimator (knows true mu values)
#'
#' @param Yhat_test Vector of surrogate predictions on the test set.
#' @param mu Vector of true per-class means.
#' @param alpha Miscoverage level for the confidence interval.
#' @export
oracle_discrete_point_and_ci <- function(Yhat_test, mu, alpha = 0.10) {
  n <- length(Yhat_test)

  # Apply true g(z) = mu[z]
  g_test <- mu[Yhat_test]

  theta_hat <- mean(g_test)
  var_hat <- var(g_test) / n

  z <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z * se,
    ci_upper = theta_hat + z * se
  )
}
