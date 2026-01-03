# Continuous Score Simulation: Core Functions
# DGP for latent scores, Model Training, Calibration Methods, Estimators
#
# DGP IN WORDS:
# - We have 10 features (X1, ..., X10)
# - True score Y depends on X1-X5 with nonlinearity: Y = f(X1,...,X5) + noise
# - We train ML models on 20% of data to predict Y from X
# - Models have varying quality:
#   - "good": knows true features (X1-X5), flexible form (GAM)
#   - "medium": uses all features, random forest (some overfitting to noise features)
#   - "poor": linear model on all features (misspecified + noise features)
#
# DGP IN MATH:
# - X_i ~ N(0, I_10) for i = 1, ..., N
# - mu(X) = beta0 + beta1*X1 + beta2*X2^2 + beta3*sin(X3) + beta4*X4*X5
# - Y = mu(X) + epsilon, epsilon ~ N(0, sigma_y^2)
# - Yhat = g(X; trained model)
#
# WHY BIAS OCCURS:
# - "poor" model uses X6-X10 which are noise → predictions regress toward mean
# - "poor" model is linear → can't capture X2^2 or sin(X3) terms
# - This creates E[Yhat|X] ≠ E[Y|X], and when aggregated, E[Yhat] may ≠ E[Y]

# ============================================================================
# 1. DATA GENERATING PROCESS
# ============================================================================

#' Generate continuous outcome data with nonlinear X dependence
#'
#' True model: Y = f(X1, ..., X5) + noise
#' X6-X10 are noise features not in the true model
#'
#' KEY DESIGN: X ~ N(mu_x, 1). Train on mu_x=1, test on mu_x=1,3,5,7,9.
#' This creates distribution shift that exposes bias in misspecified models.
#'
#' @param N Number of observations
#' @param mu_x Mean of X features (controls distribution shift)
#' @param sigma_y Noise standard deviation
#' @param seed Random seed
#' @return List with X matrix, Y vector, and true conditional mean
#' @keywords internal
generate_continuous_dgp <- function(N, mu_x = 1, sigma_y = 2.0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Generate 10 features from N(mu_x, 1)
  X <- matrix(rnorm(N * 10, mean = mu_x, sd = 1), nrow = N, ncol = 10)
  colnames(X) <- paste0("X", 1:10)

  # True conditional mean: nonlinear function of X1-X5 only
  # No intercept shift needed - E[Y] changes naturally with mu_x
  mu_X <- 0.5 * X[,1] +
    1.0 * X[,2]^2 +
    2.0 * sin(X[,3]) +
    0.8 * X[,4] * X[,5] +
    0.5 * exp(-0.3 * X[,1]^2)

  # Add noise
  Y <- mu_X + rnorm(N, 0, sigma_y)

  list(
    X = X,
    Y = Y,
    mu_X = mu_X,
    mu_x = mu_x,
    sigma_y = sigma_y
  )
}

# ============================================================================
# 2. MODEL TRAINING FUNCTIONS
# ============================================================================

#' Train "good" model: GAM on correct features X1-X5
#'
#' Knows which features matter, uses flexible functional form
#' @keywords internal
train_model_good <- function(X, Y) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' required")
  }

  df <- data.frame(Y = Y, X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4], X5 = X[,5])

  # GAM with smooth terms - can capture nonlinearity
  # Note: No interaction term, so won't perfectly capture X4*X5
  fit <- mgcv::gam(Y ~ s(X1, k = 20) + s(X2, k = 20) + s(X3, k = 20) +
                       s(X4, k = 10) + s(X5, k = 10),
                   data = df)

  function(X_new) {
    df_new <- data.frame(X1 = X_new[,1], X2 = X_new[,2], X3 = X_new[,3],
                         X4 = X_new[,4], X5 = X_new[,5])
    as.numeric(predict(fit, newdata = df_new))
  }
}

#' Train "medium" model: Random Forest on ALL features
#'
#' Uses all 10 features including noise features X6-X10
#' RF is flexible but will partially fit to noise
#' @keywords internal
train_model_medium <- function(X, Y) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' required")
  }

  df <- data.frame(Y = Y, X)

  # RF on all features - will use some noise features
  fit <- randomForest::randomForest(Y ~ ., data = df, ntree = 100)

  function(X_new) {
    df_new <- data.frame(X_new)
    colnames(df_new) <- paste0("X", 1:10)
    as.numeric(predict(fit, newdata = df_new))
  }
}

#' Train "poor-OLS" model: Linear regression on ALL features
#'
#' Uses all 10 features, linear form only
#' Can't capture nonlinearity + fits to noise features
#' NOTE: OLS still guarantees E[Yhat] = E[Y] on training distribution!
#' @keywords internal
train_model_poor_ols <- function(X, Y) {
  df <- data.frame(Y = Y, X)
  colnames(df) <- c("Y", paste0("X", 1:10))

  # Linear model - misspecified (missing nonlinear terms)
  fit <- lm(Y ~ ., data = df)

  function(X_new) {
    df_new <- data.frame(X_new)
    colnames(df_new) <- paste0("X", 1:10)
    as.numeric(predict(fit, newdata = df_new))
  }
}

#' Train "poor-tree" model: Regression tree on X1 and X6 ONLY
#'
#' Uses only X1 (relevant) and X6 (noise) - severely misspecified
#' KEY: NOT OLS, so no guarantee that E[Yhat] = E[Y] even on training data!
#' This should show bias even without distribution shift.
#' @keywords internal
train_model_poor_tree <- function(X, Y) {
  if (!requireNamespace("rpart", quietly = TRUE)) {
    stop("Package 'rpart' required")
  }

  # Only use X1 (relevant) and X6 (noise)
  df <- data.frame(Y = Y, X1 = X[,1], X6 = X[,6])

  # Regression tree - NOT OLS, won't calibrate mean
  fit <- rpart::rpart(Y ~ X1 + X6, data = df,
                      control = rpart::rpart.control(maxdepth = 5, minsplit = 20))

  function(X_new) {
    df_new <- data.frame(X1 = X_new[,1], X6 = X_new[,6])
    as.numeric(predict(fit, newdata = df_new))
  }
}

# ============================================================================
# 3. CALIBRATION METHODS FOR E[Y | Yhat]
# ============================================================================

#' Fit calibration model using linear regression
#' @export
fit_calibration_linear <- function(Y_cal, Yhat_cal) {
  df <- data.frame(Y = Y_cal, Yhat = Yhat_cal)
  fit <- lm(Y ~ Yhat, data = df)

  function(Yhat_new) {
    predict(fit, newdata = data.frame(Yhat = Yhat_new))
  }
}

#' Fit calibration model using GAM
#' @export
fit_calibration_gam <- function(Y_cal, Yhat_cal, k = 10) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' required")
  }

  df <- data.frame(Y = Y_cal, Yhat = Yhat_cal)
  k_use <- min(k, floor(length(Y_cal) / 2))
  k_use <- max(k_use, 3)

  fit <- tryCatch(
    mgcv::gam(Y ~ s(Yhat, k = k_use), data = df),
    error = function(e) lm(Y ~ Yhat, data = df)
  )

  function(Yhat_new) {
    predict(fit, newdata = data.frame(Yhat = Yhat_new))
  }
}

#' Fit calibration model using cubic spline
#' @export
fit_calibration_spline <- function(Y_cal, Yhat_cal, df = 4) {
  data_df <- data.frame(Y = Y_cal, Yhat = Yhat_cal)
  df_use <- min(df, floor(length(Y_cal) / 3))
  df_use <- max(df_use, 1)

  fit <- tryCatch({
    if (df_use >= 2) {
      lm(Y ~ splines::ns(Yhat, df = df_use), data = data_df)
    } else {
      lm(Y ~ Yhat, data = data_df)
    }
  }, error = function(e) lm(Y ~ Yhat, data = data_df))

  yhat_range <- range(Yhat_cal)

  function(Yhat_new) {
    Yhat_clamped <- pmin(pmax(Yhat_new, yhat_range[1]), yhat_range[2])
    predict(fit, newdata = data.frame(Yhat = Yhat_clamped))
  }
}

#' Fit calibration model using isotonic regression
#' @export
fit_calibration_isotonic <- function(Y_cal, Yhat_cal) {
  ord <- order(Yhat_cal)
  Yhat_sorted <- Yhat_cal[ord]
  Y_sorted <- Y_cal[ord]

  iso_fit <- isoreg(Yhat_sorted, Y_sorted)
  fitted_vals <- iso_fit$yf

  function(Yhat_new) {
    approx(x = Yhat_sorted, y = fitted_vals, xout = Yhat_new, rule = 2)$y
  }
}

# ============================================================================
# 4. ESTIMATORS FOR CONTINUOUS SCORES
# ============================================================================

#' EIF estimator for continuous surrogates with calibration
#' @export
eif_continuous_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test,
                                         calibration_method = "linear",
                                         alpha = 0.10) {
  g <- switch(calibration_method,
    "linear" = fit_calibration_linear(Y_cal, Yhat_cal),
    "gam" = fit_calibration_gam(Y_cal, Yhat_cal),
    "spline" = fit_calibration_spline(Y_cal, Yhat_cal),
    "isotonic" = fit_calibration_isotonic(Y_cal, Yhat_cal),
    stop("Unknown calibration method: ", calibration_method)
  )

  g_cal <- g(Yhat_cal)
  g_test <- g(Yhat_test)

  m <- length(Y_cal)
  n <- length(Yhat_test)

  theta_hat <- mean(g_test)

  var_g_test <- var(g_test) / n
  resid_cal <- Y_cal - g_cal
  var_cal <- var(resid_cal) / m

  var_hat <- var_g_test + var_cal
  var_hat <- max(var_hat, 1e-10)

  z <- qnorm(1 - alpha / 2)
  se <- sqrt(var_hat)

  list(
    theta = theta_hat,
    var = var_hat,
    ci_lower = theta_hat - z * se,
    ci_upper = theta_hat + z * se,
    calibration_method = calibration_method
  )
}

#' PPI estimator for continuous scores
#' @export
ppi_continuous_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test, alpha = 0.10) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  bias_hat <- mean(Y_cal - Yhat_cal)
  theta_hat <- mean(Yhat_test) + bias_hat

  eps_cal <- Y_cal - Yhat_cal
  var_eps <- if (m > 1) var(eps_cal) else 0
  var_test <- if (n > 1) var(Yhat_test) else 0
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

#' PPI++ estimator for continuous scores
#' @export
ppi_pp_continuous_point_and_ci <- function(Y_cal, Yhat_cal, Yhat_test,
                                            alpha = 0.10,
                                            lambda_range = c(0, 1)) {
  m <- length(Y_cal)
  n <- length(Yhat_test)

  var_fn <- function(lambda) {
    Z_cal <- Y_cal - lambda * Yhat_cal
    var_Z <- if (m > 1) var(Z_cal) else 0
    var_test <- if (n > 1) var(Yhat_test) else 0
    var_Z / m + lambda^2 * var_test / n
  }

  opt <- optimize(var_fn, interval = lambda_range)
  lambda_hat <- opt$minimum

  theta_hat <- mean(Y_cal) + lambda_hat * (mean(Yhat_test) - mean(Yhat_cal))
  var_hat <- max(var_fn(lambda_hat), 1e-10)

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

#' Naive estimator (just use model predictions)
#' @export
naive_continuous_point_and_ci <- function(Yhat_test, alpha = 0.10) {
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
