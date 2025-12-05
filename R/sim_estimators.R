# Helper utilities and estimators used across the simulation driver.

clamp01 <- function(x, eps = 1e-6) {
  pmin(pmax(x, eps), 1 - eps)
}

safe_var <- function(x) {
  if (length(x) <= 1) {
    0
  } else {
    stats::var(x)
  }
}

#' Euclidean projection onto the probability simplex.
#'
#' Implements the standard sorting/threshold algorithm used for sparsemax and
#' simplex projection: finds the closest probability vector to an unconstrained
#' input under the L2 norm.
#'
#' @param v Numeric vector.
#'
#' @return Numeric vector on the simplex (nonnegative, sums to 1).
#' @export
project_simplex <- function(v) {
  if (length(v) == 0) {
    return(v)
  }
  if (all(is.na(v))) {
    return(rep(1 / length(v), length(v)))
  }
  v_clean <- v
  v_clean[is.na(v_clean)] <- 0
  u <- sort(v_clean, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u > (cssv - 1) / seq_along(u)))
  theta <- (cssv[rho] - 1) / rho
  w <- pmax(v_clean - theta, 0)
  if (sum(w) <= 0) {
    rep(1 / length(v), length(v))
  } else {
    w / sum(w)
  }
}

class_one_hot <- function(y, K) {
  y <- as.integer(y)
  n <- length(y)
  mat <- matrix(0, nrow = n, ncol = K)
  ok <- !is.na(y) & y >= 1 & y <= K
  if (any(ok)) {
    mat[cbind(which(ok), y[ok])] <- 1
  }
  mat
}

#' Generate Bernoulli responses with controllable prevalence.
#'
#' @param N Integer; number of observations to simulate.
#' @param theta Target marginal prevalence for the Bernoulli outcome.
#' @param signal Numeric; strength of the latent signal.
#'
#' @return A list containing feature matrix `X`, outcome vector `Y`, and
#'   per-sample probabilities `p_true`.
#' @export
generate_dgp_data <- function(N, theta, signal = 1) {
  p <- 5L
  X <- matrix(rnorm(N * p), N, p)
  eta <- signal * (1.5 * X[, 1] - 1 * X[, 2]^2 + 0.5 * sin(X[, 3]))
  fn <- function(c) mean(plogis(eta + c)) - theta
  c_hat <- uniroot(fn, interval = c(-15, 15))$root
  prob <- plogis(eta + c_hat)
  Y <- rbinom(N, 1, prob)
  list(X = X, Y = Y, p_true = prob)
}

#' Wald-style confidence interval via logit transform.
#'
#' @param theta_hat Point estimate on probability scale.
#' @param var_hat Estimated variance of `theta_hat`.
#' @param alpha Miscoverage level for the interval.
#'
#' @return A list with `lower` and `upper` bounds.
#' @export
logit_ci <- function(theta_hat, var_hat, alpha = 0.05) {
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
ppi_point_and_ci <- function(Y_L, f_L, f_U, alpha = 0.05) {
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
                             alpha = 0.05, eps_J = 1e-4) {
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
                                        alpha = 0.05,
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

#' Rogan–Gladen style transform usable inside PPI.
#'
#' @param z Numeric vector of surrogate scores.
#' @param q0_pilot,q1_pilot Pilot specificity and sensitivity inputs.
#' @param eps Tolerance to avoid division by zero when the judge is random.
#'
#' @return Scores mapped onto the corrected probability scale.
#' @export
g_rogan_gladen <- function(z, q0_pilot, q1_pilot, eps = 1e-3) {
  J <- q0_pilot + q1_pilot - 1
  ifelse(abs(J) < eps, z, pmin(pmax((z + q0_pilot - 1) / J, 0), 1))
}

#' Posterior P(Y=1 | \\hat Z=z) under plug-in confusion matrix.
#'
#' This transform corresponds to the Bayes correction described in the
#' general PPI-with-g framework: use pilot prevalence and judge sensitivity/
#' specificity to predict the posterior probability that the underlying
#' human label is positive.
#'
#' @param z Binary surrogate outputs (0/1) or probabilities in [0,1].
#' @param q0_pilot,q1_pilot Pilot specificity/sensitivity.
#' @param theta_pilot Pilot prevalence (e.g., mean of labeled outcomes).
#' @param eps Small positive value preventing division by zero.
#'
#' @return Posterior probabilities usable as a g-transform inside PPI.
#' @export
g_posterior_prob <- function(z, q0_pilot, q1_pilot, theta_pilot, eps = 1e-6) {
  pi1 <- clamp01(theta_pilot, eps)
  pi0 <- 1 - pi1
  p_z1 <- q1_pilot * pi1 + (1 - q0_pilot) * pi0
  p_z0 <- (1 - q1_pilot) * pi1 + q0_pilot * pi0
  post_z1 <- ifelse(p_z1 < eps, pi1, (q1_pilot * pi1) / p_z1)
  post_z0 <- ifelse(p_z0 < eps, pi1, ((1 - q1_pilot) * pi1) / p_z0)
  ifelse(z >= 0.5, post_z1, post_z0)
}

#' Constrained least-squares Rogan–Gladen estimator.
#'
#' Solves \\min_{\\pi \\in \\Delta_{K-1}} ||\\hat p - \\hat M \\pi||_2^2 to enforce
#' simplex structure on the prevalence vector.
#'
#' @param M_hat Estimated confusion matrix (K x K).
#' @param p_hat Observed surrogate class probabilities (length K).
#' @param jitter Small ridge term to stabilize the quadratic program.
#'
#' @return Numeric vector of constrained prevalence estimates.
#' @export
rg_least_squares_simplex <- function(M_hat, p_hat, jitter = 1e-6) {
  K <- length(p_hat)
  D <- crossprod(M_hat) + diag(jitter, K)
  d <- crossprod(M_hat, p_hat)
  Amat <- cbind(rep(1, K), diag(K))
  bvec <- c(1, rep(0, K))
  qp <- quadprog::solve.QP(Dmat = D, dvec = d, Amat = Amat, bvec = bvec, meq = 1)
  project_simplex(qp$solution)
}

#' Identity g-transform for multiclass surrogates.
#'
#' @param S Integer vector of surrogate class labels in {1,...,K}.
#' @param context List containing at least `K`.
#'
#' @return Matrix with one-hot rows representing S.
#' @export
g_identity_multiclass <- function(S, context) {
  K <- context$K
  class_one_hot(S, K)
}

#' Confusion-matrix inverse g-transform (multiclass RG).
#'
#' @param S Integer vector of surrogate labels.
#' @param context List containing inverse confusion matrix `W_inv`.
#'
#' @return Matrix whose rows are `W_inv %*% e_S`.
#' @export
g_confusion_inverse <- function(S, context) {
  W <- context$W_inv
  K <- nrow(W)
  out <- matrix(0, nrow = length(S), ncol = K)
  for (i in seq_along(S)) {
    s_val <- as.integer(S[i])
    if (!is.na(s_val) && s_val >= 1 && s_val <= K) {
      out[i, ] <- W[, s_val]
    }
  }
  out
}

#' Estimate a multiclass confusion matrix using calibration data.
#'
#' @param S Observed surrogate labels on the calibration set.
#' @param Y True labels on the calibration set.
#' @param K Number of classes.
#' @param laplace Small positive smoothing constant added to each cell.
#'
#' @return K x K column-stochastic matrix \hat M with entries
#'   P(S=a | Y=b).
#' @export
estimate_confusion_matrix <- function(S, Y, K, laplace = 1e-3) {
  M <- matrix(laplace, nrow = K, ncol = K)
  for (b in seq_len(K)) {
    idx <- which(Y == b)
    if (length(idx) > 0) {
      counts <- tabulate(S[idx], nbins = K)
      M[, b] <- M[, b] + counts
    }
  }
  col_sums <- colSums(M)
  col_sums[col_sums == 0] <- 1
  sweep(M, 2, col_sums, "/")
}

#' General multiclass PPI estimator with vector-valued g-transform.
#'
#' @param S_unlabeled Integer vector of surrogate labels (test set).
#' @param S_labeled Integer vector of surrogate labels on calibration set.
#' @param Y_labeled Integer vector of true labels on calibration set.
#' @param g_transform Function mapping (S, context) to an n x K matrix.
#' @param context List passed into `g_transform` (must include `K`).
#' @param alpha Miscoverage level for per-class Wald intervals.
#'
#' @return A list containing per-class prevalence estimates, variance
#'   estimates, and Wald confidence limits.
#' @export
ppi_multiclass <- function(S_unlabeled,
                           S_labeled,
                           Y_labeled,
                           g_transform,
                           context,
                           alpha = 0.05,
                           project = TRUE) {
  stopifnot(length(S_labeled) == length(Y_labeled))
  K <- context$K
  g_unlabeled <- g_transform(S_unlabeled, context)
  g_labeled <- g_transform(S_labeled, context)
  if (!is.matrix(g_unlabeled)) {
    g_unlabeled <- as.matrix(g_unlabeled)
  }
  if (!is.matrix(g_labeled)) {
    g_labeled <- as.matrix(g_labeled)
  }
  Y_onehot <- class_one_hot(Y_labeled, K)
  n <- nrow(g_unlabeled)
  m <- nrow(g_labeled)
  plug_in <- colMeans(g_unlabeled)
  correction <- colMeans(Y_onehot - g_labeled)
  pi_hat <- plug_in + correction
  if (project) {
    pi_hat <- project_simplex(pi_hat)
  }
  var_unlabeled <- apply(g_unlabeled, 2, safe_var)
  var_correction <- apply(Y_onehot - g_labeled, 2, safe_var)
  var_hat <- var_unlabeled / max(n, 1) + var_correction / max(m, 1)
  z <- qnorm(1 - alpha / 2)
  se <- sqrt(pmax(var_hat, 0))
  ci_lower <- pmax(pi_hat - z * se, 0)
  ci_upper <- pmin(pi_hat + z * se, 1)
  list(
    pi_hat = pi_hat,
    var = var_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    plug_in = plug_in,
    correction = correction
  )
}

#' Multiclass PPI++ via per-class binary reducers.
#'
#' Runs the scalar PPI++ estimator separately for each class indicator
#' (1{Y=k}) using the one-hot representation of surrogate labels.
#'
#' @param S_unlabeled Surrogate labels for the large sample.
#' @param S_labeled Surrogate labels for the calibration sample.
#' @param Y_labeled True labels for the calibration sample.
#' @param K Number of classes.
#' @param alpha Miscoverage level for each class.
#'
#' @return List with per-class prevalence estimates, variance, CI limits,
#'   and the tuned lambda values.
#' @export
ppi_pp_multiclass <- function(S_unlabeled,
                              S_labeled,
                              Y_labeled,
                              K,
                              alpha = 0.05) {
  f_U <- class_one_hot(S_unlabeled, K)
  f_L <- class_one_hot(S_labeled, K)
  Y_mat <- class_one_hot(Y_labeled, K)
  pi_hat <- var_hat <- ci_lower <- ci_upper <- lambda_hat <- numeric(K)
  for (k in seq_len(K)) {
    est <- ppi_pp_point_and_ci_general(
      Y_L = Y_mat[, k],
      f_L = f_L[, k],
      f_U = f_U[, k],
      alpha = alpha
    )
    pi_hat[k] <- est$theta
    var_hat[k] <- est$var
    ci_lower[k] <- est$ci_lower
    ci_upper[k] <- est$ci_upper
    lambda_hat[k] <- est$lambda
  }
  pi_hat <- project_simplex(pi_hat)
  list(
    pi_hat = pi_hat,
    var = var_hat,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    lambda = lambda_hat
  )
}

covariate_formula <- function(response, covariates) {
  covariates <- as.data.frame(covariates)
  if (ncol(covariates) == 0) {
    stats::as.formula(paste(response, "~ 1"))
  } else {
    stats::as.formula(paste(response, "~", paste(names(covariates), collapse = " + ")))
  }
}

#' Fit covariate-dependent misclassification models.
#'
#' Estimates P(S=1 | Y=1, X) and P(S=1 | Y=0, X) via logistic regression on the
#' calibration data.
#'
#' @param covariates Data frame of covariates observed on the calibration set.
#' @param Y Calibration human labels.
#' @param S Calibration surrogate outputs (0/1).
#'
#' @return A list with two glm objects: one for the positive class (tp_model)
#'   and one for the negative class (fp_model).
#' @export
fit_misclassification_models <- function(covariates, Y, S) {
  covariates <- as.data.frame(covariates)
  data_full <- data.frame(S = S, covariates)
  pos_idx <- which(Y == 1)
  neg_idx <- which(Y == 0)
  if (length(pos_idx) < 2) {
    tp_model <- stats::glm(S ~ 1, data = data_full[pos_idx, , drop = FALSE], family = stats::binomial())
  } else {
    tp_model <- stats::glm(covariate_formula("S", covariates), data = data_full[pos_idx, , drop = FALSE], family = stats::binomial())
  }
  if (length(neg_idx) < 2) {
    fp_model <- stats::glm(S ~ 1, data = data_full[neg_idx, , drop = FALSE], family = stats::binomial())
  } else {
    fp_model <- stats::glm(covariate_formula("S", covariates), data = data_full[neg_idx, , drop = FALSE], family = stats::binomial())
  }
  list(tp_model = tp_model, fp_model = fp_model)
}

#' Covariate-adjusted Rogan--Gladen transform.
#'
#' Uses covariate-specific estimates of q0(X) and q1(X) to form an adjusted
#' surrogate via the RG formula.
#'
#' @param S Surrogate outputs (0/1) on either labeled or unlabeled data.
#' @param covariates Data frame of covariates matching `S`.
#' @param models Output from `fit_misclassification_models`.
#' @param eps Small positive number to avoid division by zero.
#'
#' @return Numeric vector of adjusted surrogates.
#' @export
g_covariate_rg <- function(S, covariates, models, eps = 1e-4) {
  covariates <- as.data.frame(covariates)
  q1_hat <- stats::predict(models$tp_model, newdata = covariates, type = "response")
  p_s1_given_y0 <- stats::predict(models$fp_model, newdata = covariates, type = "response")
  q0_hat <- 1 - p_s1_given_y0
  J_hat <- q0_hat + q1_hat - 1
  J_hat <- pmax(pmin(J_hat, 1 - eps), eps)
  g_val <- (S + q0_hat - 1) / J_hat
  clamp01(g_val, eps = eps)
}

#' Fit a RePPI-style logistic regression for g(S,X).
#'
#' @param covariates Data frame of covariates on the calibration set.
#' @param S Calibration surrogate outputs (0/1).
#' @param Y Calibration human labels (0/1).
#'
#' @return Fitted glm predicting Y from S and covariates (with interactions).
#' @export
fit_reppi_model <- function(covariates, S, Y) {
  covariates <- as.data.frame(covariates)
  if (ncol(covariates) == 0) {
    stats::glm(Y ~ S, family = stats::binomial(), data = data.frame(Y = Y, S = S))
  } else {
    terms <- paste(names(covariates), collapse = " + ")
    formula_str <- paste("Y ~ S * (", terms, ")")
    stats::glm(stats::as.formula(formula_str), family = stats::binomial(),
               data = data.frame(Y = Y, S = S, covariates))
  }
}

#' Predict using the RePPI-style logistic regression.
#'
#' @param S Surrogate outputs (0/1).
#' @param covariates Data frame of covariates.
#' @param model Fitted glm from `fit_reppi_model`.
#'
#' @return Predicted probabilities of Y=1 given (S,X).
#' @export
g_reppi_predict <- function(S, covariates, model) {
  covariates <- as.data.frame(covariates)
  preds <- stats::predict(model, newdata = data.frame(S = S, covariates), type = "response")
  clamp01(preds)
}
