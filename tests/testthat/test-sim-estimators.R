test_that("clamp01 clamps into (eps, 1-eps)", {
  x <- c(-1, 0, 0.2, 1, 2, NA_real_)
  out <- clamp01(x, eps = 1e-3)
  expect_true(all(out[!is.na(out)] >= 1e-3))
  expect_true(all(out[!is.na(out)] <= 1 - 1e-3))
  expect_true(is.na(out[is.na(x)]))
})

test_that("safe_var returns 0 for length <= 1", {
  expect_equal(safe_var(numeric(0)), 0)
  expect_equal(safe_var(1), 0)
  expect_equal(safe_var(c(1, 2)), stats::var(c(1, 2)))
})

test_that("generate_dgp_data produces correct shapes and marginal prevalence ~ theta", {
  set.seed(1)
  N <- 20000
  theta <- 0.3
  dgp <- generate_dgp_data(N, theta, signal = 1)

  expect_true(is.matrix(dgp$X))
  expect_equal(nrow(dgp$X), N)
  expect_equal(ncol(dgp$X), 5L)
  expect_equal(length(dgp$Y), N)
  expect_equal(length(dgp$p_true), N)

  # marginal mean should be close to target (Monte Carlo tolerance)
  expect_lt(abs(mean(dgp$Y) - theta), 0.01)
  expect_true(all(dgp$p_true >= 0 & dgp$p_true <= 1))
})

test_that("logit_ci returns bounds in [0,1] and lower<=upper", {
  ci <- logit_ci(theta_hat = 0.2, var_hat = 0.01, alpha = 0.05)
  expect_true(ci$lower >= 0 && ci$lower <= 1)
  expect_true(ci$upper >= 0 && ci$upper <= 1)
  expect_true(ci$lower <= ci$upper)
})

test_that("project_simplex returns nonnegative vector summing to 1", {
  set.seed(1)
  v <- rnorm(10)
  w <- project_simplex(v)

  expect_equal(length(w), length(v))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-10)
})

test_that("project_simplex handles all-NA input", {
  v <- rep(NA_real_, 4)
  w <- project_simplex(v)
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1, tolerance = 1e-10)
  expect_equal(w, rep(1/4, 4))
})

test_that("class_one_hot produces correct one-hot encoding", {
  y <- c(1, 2, 2, NA, 3, 0, 4)
  K <- 3
  mat <- class_one_hot(y, K)

  expect_equal(dim(mat), c(length(y), K))
  expect_equal(mat[1, ], c(1, 0, 0))
  expect_equal(mat[2, ], c(0, 1, 0))
  expect_equal(mat[3, ], c(0, 1, 0))
  expect_equal(mat[4, ], c(0, 0, 0)) # NA -> all zeros
  expect_equal(mat[5, ], c(0, 0, 1))
  expect_equal(mat[6, ], c(0, 0, 0)) # out of range -> all zeros
  expect_equal(mat[7, ], c(0, 0, 0)) # out of range -> all zeros
})

test_that("estimate_confusion_matrix is column-stochastic", {
  set.seed(1)
  K <- 3
  Y <- sample(1:K, 100, replace = TRUE)
  # noisy surrogate labels
  S <- pmin(pmax(Y + sample(c(-1,0,1), 100, replace = TRUE), 1), K)

  Mhat <- estimate_confusion_matrix(S, Y, K, laplace = 1e-3)
  expect_equal(dim(Mhat), c(K, K))
  expect_true(all(Mhat >= 0))
  expect_equal(colSums(Mhat), rep(1, K), tolerance = 1e-10)
})

test_that("ppi_point_and_ci basic sanity", {
  set.seed(1)
  n <- 200
  N <- 2000
  Y_L <- rbinom(n, 1, 0.3)
  # Surrogate with some noise
  f_L <- pmin(pmax(Y_L + rnorm(n, sd = 0.3), 0), 1)
  f_U <- pmin(pmax(rnorm(N, mean = 0.3, sd = 0.2), 0), 1)

  est <- ppi_point_and_ci(Y_L, f_L, f_U, alpha = 0.05)

  expect_true(is.finite(est$theta))
  expect_true(est$theta >= 0 && est$theta <= 1)
  expect_true(est$var >= 0)
  expect_true(est$ci_lower <= est$ci_upper)
  expect_true(est$ci_lower >= 0 && est$ci_upper <= 1)
})

test_that("rg_point_and_ci returns NA when judge nearly random (J ~ 0)", {
  # q0 + q1 - 1 ~ 0
  out <- rg_point_and_ci(
    p_hat = 0.5, q0_hat = 0.5, q1_hat = 0.5,
    n = 1000, m0 = 50, m1 = 50,
    alpha = 0.05, eps_J = 1e-4
  )
  expect_true(is.na(out$theta))
  expect_true(is.na(out$var))
})

test_that("rg_point_and_ci returns bounded estimate and CI when J is safe", {
  out <- rg_point_and_ci(
    p_hat = 0.4, q0_hat = 0.9, q1_hat = 0.8,
    n = 1000, m0 = 200, m1 = 200,
    alpha = 0.05
  )
  expect_true(out$theta >= 0 && out$theta <= 1)
  expect_true(out$var >= 0)
  expect_true(out$ci_lower <= out$ci_upper)
  expect_true(out$ci_lower >= 0 && out$ci_upper <= 1)
})

test_that("ppi_pp_point_and_ci returns lambda within range and bounded theta", {
  set.seed(2)
  n <- 300
  N <- 3000
  Y_L <- rbinom(n, 1, 0.4)
  f_L <- rbinom(n, 1, 0.7)
  f_U <- rbinom(N, 1, 0.6)

  est <- ppi_pp_point_and_ci(Y_L, f_L, f_U, alpha = 0.05)

  # Lambda is now unconstrained (closed-form optimal)
  expect_true(is.numeric(est$lambda))
  expect_true(est$theta >= 0 && est$theta <= 1)
  expect_true(est$var >= 0)
  expect_true(est$ci_lower <= est$ci_upper)
  expect_true(est$ci_lower >= 0 && est$ci_upper <= 1)
})

test_that("lambda_hat_ppi_pp_optimal finds variance-minimizing lambda", {
  set.seed(3)
  n <- 200
  N <- 2000
  Y_L <- rbinom(n, 1, 0.3)
  f_L <- rbinom(n, 1, 0.6)
  f_U <- rbinom(N, 1, 0.5)

  # Get optimal lambda from closed-form solution
  lam_hat <- lambda_hat_ppi_pp_optimal(Y_L, f_L, f_U)

  # Evaluate variance on a wide grid (lambda is now unconstrained)
  grid <- seq(-1, 3, length.out = 201)
  vars <- vapply(
    grid,
    function(l) ppi_pp_var_hat(l, Y_L, f_L, f_U),
    numeric(1)
  )

  v_hat <- ppi_pp_var_hat(lam_hat, Y_L, f_L, f_U)

 # λ_hat should be near the global minimum
  expect_true(
    v_hat <= min(vars) + 1e-6
  )
})

test_that("ppi_multiclass returns simplex prevalence and valid CIs", {
  set.seed(1)
  K <- 3
  nU <- 5000
  nL <- 300

  # True prevalence
  pi <- c(0.2, 0.5, 0.3)

  Y_L <- sample(1:K, nL, replace = TRUE, prob = pi)
  S_L <- sample(1:K, nL, replace = TRUE, prob = pi)  # (identity surrogate for test)
  S_U <- sample(1:K, nU, replace = TRUE, prob = pi)

  ctx <- list(K = K)
  out <- ppi_multiclass(
    S_unlabeled = S_U,
    S_labeled = S_L,
    Y_labeled = Y_L,
    g_transform = g_identity_multiclass,
    context = ctx,
    alpha = 0.05,
    project = TRUE
  )

  expect_equal(length(out$pi_hat), K)
  expect_true(all(out$pi_hat >= 0))
  expect_equal(sum(out$pi_hat), 1, tolerance = 1e-10)

  expect_equal(length(out$ci_lower), K)
  expect_equal(length(out$ci_upper), K)
  expect_true(all(out$ci_lower <= out$ci_upper))
  expect_true(all(out$ci_lower >= 0))
  expect_true(all(out$ci_upper <= 1))
})

test_that("ppi_pp_multiclass returns simplex prevalence and lambda per class", {
  set.seed(2)
  K <- 4
  nU <- 3000
  nL <- 400
  pi <- rep(1/K, K)

  Y_L <- sample(1:K, nL, replace = TRUE, prob = pi)
  S_L <- sample(1:K, nL, replace = TRUE, prob = pi)
  S_U <- sample(1:K, nU, replace = TRUE, prob = pi)

  out <- ppi_pp_multiclass(S_U, S_L, Y_L, K, alpha = 0.05)

  expect_equal(length(out$pi_hat), K)
  expect_true(all(out$pi_hat >= 0))
  expect_equal(sum(out$pi_hat), 1, tolerance = 1e-10)

  expect_equal(length(out$lambda), K)
  expect_true(all(out$lambda >= 0 & out$lambda <= 1))
})

test_that("rg_least_squares_simplex returns simplex vector", {
  skip_if_not_installed("quadprog")

  set.seed(1)
  K <- 3
  # confusion matrix (column-stochastic): P(S=a | Y=b)
  M_hat <- matrix(c(
    0.8, 0.1, 0.1,
    0.1, 0.8, 0.1,
    0.1, 0.1, 0.8
  ), nrow = K, byrow = TRUE)

  # observed surrogate distribution
  p_hat <- c(0.3, 0.5, 0.2)

  pi_hat <- rg_least_squares_simplex(M_hat, p_hat)

  expect_equal(length(pi_hat), K)
  expect_true(all(pi_hat >= 0))
  expect_equal(sum(pi_hat), 1, tolerance = 1e-10)
})

test_that("rg_least_squares_simplex approximately solves least squares objective", {
  skip_if_not_installed("quadprog")

  K <- 3
  M_hat <- diag(K)  # identity s.t. solution should be close to p_hat projected
  p_hat <- c(0.2, 0.3, 0.5)

  pi_hat <- rg_least_squares_simplex(M_hat, p_hat)
  expect_equal(pi_hat, p_hat, tolerance = 1e-6)
})

test_that("eif_point_and_ci returns bounded estimate and CI for binary surrogate", {
  set.seed(123)

  # Generate test data
  n_test <- 1000
  m_cal <- 100
  theta_true <- 0.3
  q0 <- 0.8  # specificity
  q1 <- 0.9  # sensitivity

  # True labels
  Y_cal <- rbinom(m_cal, 1, theta_true)

  # Surrogate predictions
  Yhat_test <- rbinom(n_test, 1, theta_true * q1 + (1 - theta_true) * (1 - q0))
  Yhat_cal <- ifelse(Y_cal == 1, rbinom(m_cal, 1, q1), rbinom(m_cal, 1, 1 - q0))

  # Run EIF estimator
  result <- eif_point_and_ci(Y_cal, Yhat_cal, Yhat_test)

  expect_true(result$theta >= 0 && result$theta <= 1)
  expect_true(result$var >= 0)
  expect_true(result$ci_lower <= result$ci_upper)
  expect_true(result$ci_lower >= 0 && result$ci_upper <= 1)
  expect_true(result$mu0_hat >= 0 && result$mu0_hat <= 1)
  expect_true(result$mu1_hat >= 0 && result$mu1_hat <= 1)
  expect_true(result$p_hat >= 0 && result$p_hat <= 1)
})

test_that("eif_point_and_ci returns NA when no calibration data for a stratum", {
  # All calibration samples have Yhat = 1 (no Yhat = 0)
  Y_cal <- c(1, 1, 0, 1)
  Yhat_cal <- c(1, 1, 1, 1)  # all ones
  Yhat_test <- c(0, 1, 1, 0, 1)

  result <- eif_point_and_ci(Y_cal, Yhat_cal, Yhat_test)

  expect_true(is.na(result$theta))
  expect_true(is.na(result$var))
})

test_that("eif_point_and_ci achieves approximately correct coverage", {
  set.seed(42)
  n_reps <- 200
  n_test <- 500
  m_cal <- 50
  theta_true <- 0.4
  q0 <- 0.85
  q1 <- 0.80

  covered <- 0
  for (rep in 1:n_reps) {
    Y_cal <- rbinom(m_cal, 1, theta_true)
    Yhat_cal <- ifelse(Y_cal == 1, rbinom(m_cal, 1, q1), rbinom(m_cal, 1, 1 - q0))
    Yhat_test <- rbinom(n_test, 1, theta_true * q1 + (1 - theta_true) * (1 - q0))

    result <- eif_point_and_ci(Y_cal, Yhat_cal, Yhat_test, alpha = 0.05)

    if (!is.na(result$ci_lower) && !is.na(result$ci_upper)) {
      if (theta_true >= result$ci_lower && theta_true <= result$ci_upper) {
        covered <- covered + 1
      }
    }
  }

  coverage_rate <- covered / n_reps
  # Coverage should be approximately 0.95 (allow margin for finite sample)
  expect_true(coverage_rate >= 0.85 && coverage_rate <= 1.0)
})