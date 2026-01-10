test_that("clamp01 clamps into (eps, 1-eps)", {
  x <- c(-1, 0, 0.2, 1, 2, NA_real_)
  out <- debiasLLMReporting:::clamp01(x, eps = 1e-3)
  expect_true(all(out[!is.na(out)] >= 1e-3))
  expect_true(all(out[!is.na(out)] <= 1 - 1e-3))
  expect_true(is.na(out[is.na(x)]))
})

test_that("generate_dgp_data produces correct length and marginal prevalence ~ theta", {
  set.seed(1)
  N <- 20000
  theta <- 0.3
  dgp <- generate_dgp_data(N, theta, signal = 1)
  expect_equal(length(dgp$Y), N)
  # marginal mean should be close to target (Monte Carlo tolerance)
  expect_lt(abs(mean(dgp$Y) - theta), 0.01)
})

test_that("logit_ci returns bounds in [0,1] and lower<=upper", {
  ci <- logit_ci(theta_hat = 0.2, var_hat = 0.01, alpha = 0.05)
  expect_true(ci$lower >= 0 && ci$lower <= 1)
  expect_true(ci$upper >= 0 && ci$upper <= 1)
  expect_true(ci$lower <= ci$upper)
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

test_that("ppi_point_and_ci handles empty and length-1 inputs", {
  # Empty inputs lead to NaN outputs (current behavior)
  out_empty <- ppi_point_and_ci(numeric(0), numeric(0), numeric(0), alpha = 0.05)
  expect_true(is.nan(out_empty$theta))
  expect_true(is.nan(out_empty$var))
  expect_true(is.nan(out_empty$ci_lower))
  expect_true(is.nan(out_empty$ci_upper))

  # Length-1 inputs should be finite with zero variance
  out_one <- ppi_point_and_ci(Y_L = 1, f_L = 1, f_U = 0, alpha = 0.05)
  expect_true(is.finite(out_one$theta))
  expect_true(is.finite(out_one$var))
  expect_true(is.finite(out_one$ci_lower))
  expect_true(is.finite(out_one$ci_upper))
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
  lam_hat <- debiasLLMReporting:::lambda_hat_ppi_pp_optimal(Y_L, f_L, f_U)

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

test_that("lambda_hat_ppi_pp_optimal falls back when f_L/f_U are constant", {
  Y_L <- c(0, 1, 1, 0, 1)
  f_L <- rep(0.5, length(Y_L))
  f_U <- rep(0.5, 20)

  lam_hat <- lambda_hat_ppi_pp_optimal(Y_L, f_L, f_U)
  expect_equal(lam_hat, 1)

  est <- ppi_pp_point_and_ci(Y_L, f_L, f_U, alpha = 0.05)
  expect_true(is.finite(est$theta))
  expect_true(is.finite(est$var))
  expect_true(is.finite(est$ci_lower))
  expect_true(is.finite(est$ci_upper))
})

test_that("ppi_pp_var_hat handles empty and length-1 inputs", {
  # Empty inputs lead to NaN outputs (current behavior)
  var_empty <- ppi_pp_var_hat(1, numeric(0), numeric(0), numeric(0))
  expect_true(is.nan(var_empty))

  # Length-1 inputs yield zero variance
  var_one <- ppi_pp_var_hat(1, Y_L = 1, f_L = 1, f_U = 0)
  expect_equal(var_one, 0)
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

test_that("eif_point_and_ci filters NA and non-binary inputs", {
  Y_cal <- c(1, 0, 1, NA, 2)
  Yhat_cal <- c(1, 0, 1, NA, 2)
  Yhat_test <- c(1, 0, 1, NA, 2, 0)

  result <- eif_point_and_ci(Y_cal, Yhat_cal, Yhat_test)

  expect_true(is.finite(result$theta))
  expect_true(is.finite(result$var))
  expect_true(result$theta >= 0 && result$theta <= 1)
  expect_true(result$ci_lower <= result$ci_upper)
  expect_true(result$ci_lower >= 0 && result$ci_upper <= 1)
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

test_that("fit_misclass_mle recovers parameters in a simple simulation", {
  set.seed(7)
  N <- 5000
  m_cal <- 500
  theta_true <- 0.4
  q0_true <- 0.85
  q1_true <- 0.8

  Y_all <- rbinom(N, 1, theta_true)
  Yhat_all <- ifelse(
    Y_all == 1,
    rbinom(N, 1, q1_true),
    rbinom(N, 1, 1 - q0_true)
  )

  idx_cal <- sample(seq_len(N), m_cal)
  y_cal <- Y_all[idx_cal]
  yhat_cal <- Yhat_all[idx_cal]
  yhat_test <- Yhat_all[-idx_cal]

  fit <- fit_misclass_mle(y_cal, yhat_cal, yhat_test, level = 0.90)

  expect_true(is.finite(fit$logLik))
  expect_true(abs(fit$theta_hat - theta_true) <= 0.1)
  expect_true(abs(fit$q0_hat - q0_true) <= 0.1)
  expect_true(abs(fit$q1_hat - q1_true) <= 0.1)
})
