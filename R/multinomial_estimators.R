###############################################################
# Multinomial (K-class) LLM-as-a-judge calibration utilities
#   - Multiclass RG (known M / plug-in M)
#   - Multiclass PPI (g = onehot) and EIF/RePPI-style (g = mu_hat)
#   - Optional: Joint MLE (pi, M) in the misclassification model via EM
#
# Conventions (IMPORTANT):
#   * Labels are integers in {1,2,...,K}.
#   * N[a, k] = count{ y_hat = a, y_true = k }  (rows = judge, cols = true)
#   * Confusion matrix M[a, k] = P(y_hat = a | y_true = k) (columns sum to 1)
#   * mu[a, k] = P(y_true = k | y_hat = a) (rows sum to 1)
###############################################################

# -----------------------------
# Basic checks + encodings
# -----------------------------

as_labels <- function(x, K = NULL, allow_na = FALSE, name = "x") {
  if (is.factor(x)) x <- as.integer(x)
  if (!is.numeric(x)) stop(sprintf("%s must be numeric/integer/factor.", name))

  if (!allow_na && any(is.na(x))) stop(sprintf("%s contains NA.", name))
  x <- as.integer(x)
  if (!is.null(K)) {
    bad <- (!is.na(x)) & (x < 1L | x > K)
    if (any(bad)) stop(sprintf("%s has values outside 1..K.", name))
  } else {
    if (any(x < 1L, na.rm = TRUE)) stop(sprintf("%s has values < 1.", name))
  }
  x
}

onehot <- function(labels, K) {
  labels <- as_labels(labels, K = K, name = "labels")
  n <- length(labels)
  M <- matrix(0.0, nrow = n, ncol = K)
  M[cbind(seq_len(n), labels)] <- 1.0
  M
}

count_test <- function(y_hat, K) {
  y_hat <- as_labels(y_hat, K = K, name = "y_hat")
  tab <- tabulate(y_hat, nbins = K)
  as.integer(tab)
}

count_calibration_table <- function(y_true, y_hat, K) {
  y_true <- as_labels(y_true, K = K, name = "y_true")
  y_hat  <- as_labels(y_hat,  K = K, name = "y_hat")
  if (length(y_true) != length(y_hat)) stop("y_true and y_hat must have same length.")
  N <- matrix(0L, nrow = K, ncol = K)
  # rows = y_hat (judge), cols = y_true
  for (i in seq_along(y_true)) {
    N[y_hat[i], y_true[i]] <- N[y_hat[i], y_true[i]] + 1L
  }
  N
}

prep_counts <- function(y_hat_test = NULL, T = NULL,
                        y_true_cal = NULL, y_hat_cal = NULL, N = NULL,
                        K = NULL) {
  # Infer K if not provided
  if (is.null(K)) {
    candidates <- integer(0)
    if (!is.null(y_hat_test)) candidates <- c(candidates, as.integer(max(y_hat_test, na.rm = TRUE)))
    if (!is.null(y_true_cal)) candidates <- c(candidates, as.integer(max(y_true_cal, na.rm = TRUE)))
    if (!is.null(y_hat_cal))  candidates <- c(candidates,  as.integer(max(y_hat_cal, na.rm = TRUE)))
    if (!is.null(T))          candidates <- c(candidates, length(T))
    if (!is.null(N))          candidates <- c(candidates, nrow(N), ncol(N))
    if (length(candidates) == 0) stop("Cannot infer K; please provide K.")
    K <- max(candidates)
  }

  # Test counts
  if (is.null(T)) {
    if (is.null(y_hat_test)) stop("Provide either y_hat_test or T.")
    T <- count_test(y_hat_test, K)
  } else {
    if (length(T) != K) stop("T must have length K.")
    if (any(T < 0)) stop("T must be nonnegative.")
    T <- as.integer(T)
  }

  # Calibration table
  if (is.null(N)) {
    if (is.null(y_true_cal) || is.null(y_hat_cal)) stop("Provide either (y_true_cal,y_hat_cal) or N.")
    N <- count_calibration_table(y_true_cal, y_hat_cal, K)
  } else {
    if (!is.matrix(N) || nrow(N) != K || ncol(N) != K) stop("N must be a KxK matrix.")
    if (any(N < 0)) stop("N must be nonnegative.")
    N <- matrix(as.integer(N), nrow = K, ncol = K)
  }

  list(K = K, T = T, N = N, n = sum(T), m = sum(N))
}

# -----------------------------
# Simplex projection (for stability / post-processing)
# -----------------------------
proj_simplex <- function(v) {
  # Euclidean projection onto {x >= 0, sum x = 1}
  v <- as.numeric(v)
  if (any(!is.finite(v))) stop("proj_simplex: v must be finite.")
  K <- length(v)
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u - (cssv - 1) / seq_len(K) > 0))
  theta <- (cssv[rho] - 1) / rho
  w <- pmax(v - theta, 0)
  # numerical guard
  s <- sum(w)
  if (s <= 0) rep(1 / K, K) else w / s
}

# -----------------------------
# Estimating M and mu from calibration table
# -----------------------------

estimate_M_from_N <- function(N, alpha = 0.5) {
  # M[a,k] = P(y_hat=a | y_true=k), columns sum to 1
  if (!is.matrix(N)) stop("N must be a matrix.")
  K <- nrow(N)
  if (ncol(N) != K) stop("N must be KxK.")
  if (alpha < 0) stop("alpha must be >= 0.")
  col_tot <- colSums(N)
  denom <- col_tot + K * alpha
  M <- (N + alpha)
  M <- sweep(M, 2, denom, "/")
  list(M = M, m_k = as.integer(col_tot))
}

estimate_mu_from_N <- function(N, alpha = 0.5) {
  # mu[a,k] = P(y_true=k | y_hat=a), rows sum to 1
  if (!is.matrix(N)) stop("N must be a matrix.")
  K <- nrow(N)
  if (ncol(N) != K) stop("N must be KxK.")
  if (alpha < 0) stop("alpha must be >= 0.")
  row_tot <- rowSums(N)
  denom <- row_tot + K * alpha
  mu <- (N + alpha)
  mu <- sweep(mu, 1, denom, "/")
  list(mu = mu, r_a = as.integer(row_tot))
}

# -----------------------------
# Covariances for multinomial proportions
# -----------------------------

cov_multinom_phat <- function(p_hat, n) {
  # Cov( p_hat ) for p_hat = empirical proportions of a K-category variable
  p_hat <- as.numeric(p_hat)
  (diag(p_hat) - tcrossprod(p_hat)) / n
}

# -----------------------------
# Multiclass RG (known M)
# -----------------------------

multiclass_rg_known <- function(y_hat_test = NULL, T = NULL, M,
                                project = FALSE) {
  if (!is.matrix(M)) stop("M must be a KxK matrix.")
  K <- nrow(M)
  if (ncol(M) != K) stop("M must be KxK.")
  pre <- prep_counts(y_hat_test = y_hat_test, T = T,
                     y_true_cal = rep(1L, K), y_hat_cal = rep(1L, K), # dummy
                     N = diag(0L, K), K = K)
  T <- pre$T
  n <- pre$n
  p_hat <- as.numeric(T) / n

  # Solve M pi = p
  # NOTE: M must be invertible for RG-known
  if (kappa(M) > 1e10) warning("M is ill-conditioned; inversion may be unstable.")
  pi_hat <- as.numeric(qr.solve(M, p_hat))

  if (project) pi_hat <- proj_simplex(pi_hat)

  # Asymptotic covariance under known M:
  A <- qr.solve(M, diag(1, K))  # M^{-1}
  cov_p <- cov_multinom_phat(p_hat, n)
  cov_pi <- A %*% cov_p %*% t(A)

  list(
    method = "RG_known",
    K = K,
    n = n,
    pi_hat = pi_hat,
    p_hat = p_hat,
    cov = cov_pi,
    se = sqrt(pmax(diag(cov_pi), 0))
  )
}

# -----------------------------
# Multiclass RG plug-in (estimate M from calibration, then invert)
#   Includes delta-method covariance ignoring randomness of m_k
# -----------------------------

multiclass_rg_plugin <- function(y_hat_test = NULL, T = NULL,
                                 y_true_cal = NULL, y_hat_cal = NULL, N = NULL,
                                 K = NULL,
                                 alpha = 0.5,
                                 project = FALSE,
                                 warn_kappa = 1e10) {

  pre <- prep_counts(y_hat_test = y_hat_test, T = T,
                     y_true_cal = y_true_cal, y_hat_cal = y_hat_cal, N = N,
                     K = K)
  K <- pre$K; T <- pre$T; N <- pre$N; n <- pre$n; m <- pre$m

  p_hat <- as.numeric(T) / n

  estM <- estimate_M_from_N(N, alpha = alpha)
  M_hat <- estM$M
  m_k <- estM$m_k

  # Invert M_hat
  kap <- kappa(M_hat)
  if (kap > warn_kappa) warning(sprintf("M_hat is ill-conditioned (kappa=%.2e). Consider smoothing (alpha) or ridge/CLS.", kap))
  A_hat <- qr.solve(M_hat, diag(1, K))  # M_hat^{-1}

  pi_hat <- as.numeric(A_hat %*% p_hat)
  if (project) pi_hat <- proj_simplex(pi_hat)

  # Delta-method covariance:
  # Cov(p_hat) + contribution from estimating columns of M_hat
  cov_p <- cov_multinom_phat(p_hat, n)

  # Column-wise contribution from M_hat uncertainty:
  # If column k is estimated from m_k samples:
  #   Cov(M_hat[,k]) ≈ (1/m_k)(diag(M_hat[,k]) - M_hat[,k] M_hat[,k]^T)
  # Perturbation: d pi ≈ - M^{-1} (dM) pi, so column k contributes pi_k^2 * A Cov(col_k) A^T
  S <- matrix(0.0, nrow = K, ncol = K)
  for (k in seq_len(K)) {
    if (m_k[k] <= 0) next
    Mk <- as.numeric(M_hat[, k])
    cov_Mk <- (diag(Mk) - tcrossprod(Mk)) / m_k[k]
    S <- S + (pi_hat[k]^2) * cov_Mk
  }

  cov_pi <- A_hat %*% (cov_p + S) %*% t(A_hat)

  list(
    method = "RG_plugin",
    K = K,
    n = n,
    m = m,
    pi_hat = pi_hat,
    p_hat = p_hat,
    M_hat = M_hat,
    m_k = m_k,
    cov = cov_pi,
    se = sqrt(pmax(diag(cov_pi), 0)),
    kappa_Mhat = kap
  )
}

# -----------------------------
# Multiclass PPI / EIF-style estimator for pi
#   General form:
#     pi_hat = mean_test g(y_hat) + mean_cal (Y_onehot - g(y_hat))
#   Choices:
#     g = "onehot" -> g(a)=e(a)   (vanilla PPI)
#     g = "mu"     -> g(a)=mu_hat(a)=P(Y=· | y_hat=a)  (EIF / RePPI-style)
# -----------------------------

multiclass_ppi <- function(y_hat_test = NULL, T = NULL,
                           y_true_cal = NULL, y_hat_cal = NULL, N = NULL,
                           K = NULL,
                           g = c("onehot", "mu"),
                           alpha = 0.5,
                           project = FALSE) {

  g <- match.arg(g)
  pre <- prep_counts(y_hat_test = y_hat_test, T = T,
                     y_true_cal = y_true_cal, y_hat_cal = y_hat_cal, N = N,
                     K = K)
  K <- pre$K; T <- pre$T; N <- pre$N; n <- pre$n; m <- pre$m

  # Reconstruct label vectors ONLY if user provided them.
  # If user passes only counts, we can still compute the estimator via expectations:
  # - For g="onehot", mean_test g = p_hat, and mean_cal (Y - g) can be computed from N.
  # - For g="mu", mean_test mu(y_hat) = sum_a p_hat[a] mu[a,], and mean_cal (Y - mu(y_hat)) from N too.
  p_hat <- as.numeric(T) / n

  if (g == "onehot") {
    # mean_test e(y_hat) = p_hat
    mean_g_test <- p_hat

    # mean_cal Y_onehot is true class proportions in calibration
    pi_cal <- as.numeric(colSums(N)) / m

    # mean_cal e(y_hat) is judge proportions within calibration
    p_cal <- as.numeric(rowSums(N)) / m  # NOTE: row sums are y_hat counts in calibration

    pi_hat <- mean_g_test + (pi_cal - p_cal)

    # Covariance estimate:
    # Cov(g_test)/n where g_test is onehot -> analytic
    cov_test <- cov_multinom_phat(p_hat, n) * n  # careful: cov_multinom_phat gives Cov(mean); multiply by n to get Cov(g)
    # Better: directly Cov(mean_g_test)= (diag(p_hat)-p_hat p_hat^T)/n
    cov_mean_g_test <- cov_multinom_phat(mean_g_test, n)

    # For calibration residual R = Y_onehot - e(y_hat),
    # estimate Cov(R) from counts N (exact sample covariance from cell probabilities).
    # Compute empirical residual vectors for each (a,k) cell:
    # r(a,k) = e(k) - e(a), with probability N[a,k]/m
    # Then Cov(mean residual) = Cov(R)/m.
    # We'll compute Cov(R) from the distribution over cells.
    probs <- as.numeric(N) / m
    # Build all residual vectors efficiently:
    # For each cell (a,k), residual = e(k) - e(a)
    # E[residual] already equals pi_cal - p_cal
    mean_res <- (pi_cal - p_cal)
    CovR <- matrix(0.0, K, K)
    for (a in seq_len(K)) {
      ea <- rep(0.0, K); ea[a] <- 1.0
      for (k in seq_len(K)) {
        pk <- probs[a + (k - 1L) * K]
        if (pk <= 0) next
        ek <- rep(0.0, K); ek[k] <- 1.0
        r <- ek - ea
        d <- r - mean_res
        CovR <- CovR + pk * tcrossprod(d)
      }
    }
    cov_mean_res <- CovR / m

    cov_pi <- cov_mean_g_test + cov_mean_res

    if (project) pi_hat <- proj_simplex(pi_hat)

    return(list(
      method = "PPI_onehot",
      K = K, n = n, m = m,
      pi_hat = pi_hat,
      cov = cov_pi,
      se = sqrt(pmax(diag(cov_pi), 0)),
      p_hat = p_hat
    ))
  }

  # g == "mu": estimate mu[a,k] = P(Y=k | y_hat=a) from calibration N
  estmu <- estimate_mu_from_N(N, alpha = alpha)
  mu_hat <- estmu$mu

  # mean_test mu(y_hat) = sum_a p_hat[a] * mu_hat[a,]
  mean_g_test <- as.numeric(t(p_hat) %*% mu_hat)

  # mean_cal of Y_onehot = pi_cal
  pi_cal <- as.numeric(colSums(N)) / m

  # mean_cal mu(y_hat) = sum_a p_cal[a] * mu_hat[a,]
  p_cal <- as.numeric(rowSums(N)) / m
  mean_g_cal <- as.numeric(t(p_cal) %*% mu_hat)

  pi_hat <- mean_g_test + (pi_cal - mean_g_cal)

  # Covariance estimate via plug-in sample covariances from cell distribution:
  # Test: g(y_hat)=mu_hat[row=y_hat] takes only K distinct values -> exact covariance from p_hat
  #   Cov(g) = sum_a p_hat[a] (mu_a - mean)^T(mu_a - mean)
  #   Cov(mean_g_test) = Cov(g)/n
  CovG_test <- matrix(0.0, K, K)
  for (a in seq_len(K)) {
    da <- mu_hat[a, ] - mean_g_test
    CovG_test <- CovG_test + p_hat[a] * tcrossprod(da)
  }
  cov_mean_g_test <- CovG_test / n

  # Calibration residual: R = Y_onehot - mu_hat[y_hat]
  # Distribution over cells (a,k) with prob N[a,k]/m:
  #   r(a,k) = e(k) - mu_hat[a,]
  # Compute Cov(R) and divide by m
  mean_res <- (pi_cal - mean_g_cal)
  CovR <- matrix(0.0, K, K)
  probs <- as.numeric(N) / m
  for (a in seq_len(K)) {
    mu_a <- mu_hat[a, ]
    for (k in seq_len(K)) {
      pk <- probs[a + (k - 1L) * K]
      if (pk <= 0) next
      ek <- rep(0.0, K); ek[k] <- 1.0
      r <- ek - mu_a
      d <- r - mean_res
      CovR <- CovR + pk * tcrossprod(d)
    }
  }
  cov_mean_res <- CovR / m

  cov_pi <- cov_mean_g_test + cov_mean_res

  if (project) pi_hat <- proj_simplex(pi_hat)

  list(
    method = "PPI_mu(EIF/RePPI-style)",
    K = K, n = n, m = m,
    pi_hat = pi_hat,
    cov = cov_pi,
    se = sqrt(pmax(diag(cov_pi), 0)),
    p_hat = p_hat,
    mu_hat = mu_hat
  )
}

# -----------------------------
# Optional: Joint MLE in Model B via EM
#   Misclassification model: P(y_hat=a, y_true=k) = pi_k * M[a,k]
#   Test data provides only y_hat counts T; calibration provides N.
# -----------------------------

multiclass_mle_em <- function(T, N, alpha = 0.5,
                              maxit = 2000, tol = 1e-10,
                              verbose = FALSE) {
  if (!is.matrix(N)) stop("N must be a KxK matrix.")
  K <- nrow(N)
  if (ncol(N) != K) stop("N must be KxK.")
  if (length(T) != K) stop("T must have length K.")
  if (any(T < 0) || any(N < 0)) stop("Counts must be nonnegative.")
  if (alpha < 0) stop("alpha must be >= 0.")

  T <- as.numeric(T)
  N <- matrix(as.numeric(N), nrow = K, ncol = K)
  n <- sum(T); m <- sum(N); Ntot <- n + m

  # Initialize from calibration with smoothing
  estM <- estimate_M_from_N(N, alpha = alpha)
  M <- estM$M
  pi <- as.numeric(colSums(N) + alpha)  # mild smoothing
  pi <- pi / sum(pi)

  loglik <- function(pi, M) {
    # loglik = sum_{a,k} N_{a,k} log(pi_k M_{a,k}) + sum_a T_a log(p_a)
    p <- as.numeric(M %*% pi)
    # guard
    eps <- 1e-15
    p <- pmax(p, eps)
    val1 <- 0
    for (a in seq_len(K)) for (k in seq_len(K)) {
      if (N[a,k] > 0) val1 <- val1 + N[a,k] * (log(pmax(pi[k], eps)) + log(pmax(M[a,k], eps)))
    }
    val2 <- sum(T * log(p))
    val1 + val2
  }

  ll_old <- loglik(pi, M)

  for (it in seq_len(maxit)) {
    # E-step on test counts:
    # w[a,k] = P(Y=k | y_hat=a) under current params = pi_k M[a,k] / p[a]
    p <- as.numeric(M %*% pi)
    eps <- 1e-15
    p <- pmax(p, eps)

    W <- sweep(M, 2, pi, "*")   # multiply each column k by pi_k
    W <- sweep(W, 1, p, "/")    # divide each row a by p_a
    # expected test cell counts:
    Ntest <- sweep(W, 1, T, "*")  # row a multiplied by T_a

    # M-step: combine observed calibration counts with expected test counts
    Nall <- N + Ntest

    # update pi
    pi_new <- colSums(Nall)
    pi_new <- pi_new / sum(pi_new)

    # update M columns
    col_tot <- colSums(Nall)
    # avoid division by zero
    col_tot <- pmax(col_tot, eps)
    M_new <- sweep(Nall, 2, col_tot, "/")

    # check convergence
    ll_new <- loglik(pi_new, M_new)
    if (verbose && (it %% 50 == 0 || it == 1)) {
      cat(sprintf("iter=%d  loglik=%.6f  dloglik=%.3e\n", it, ll_new, ll_new - ll_old))
    }
    if (abs(ll_new - ll_old) <= tol * (1 + abs(ll_old))) {
      pi <- pi_new; M <- M_new; ll_old <- ll_new
      break
    }

    pi <- pi_new
    M <- M_new
    ll_old <- ll_new
  }

  list(
    method = "MLE_EM(Model B)",
    K = K, n = n, m = m,
    pi_hat = as.numeric(pi),
    M_hat = M,
    loglik = ll_old
  )
}

# -----------------------------
# Naive estimator (no correction)
# -----------------------------

naive_multinomial <- function(y_hat_test = NULL, T = NULL, K = NULL) {
  # Handle both raw vectors and count vectors
  if (is.null(T)) {
    if (is.null(y_hat_test)) stop("Provide either y_hat_test or T.")
    if (is.null(K)) K <- max(y_hat_test, na.rm = TRUE)
    T <- count_test(y_hat_test, K)
  } else {
    if (is.null(K)) K <- length(T)
  }

  n <- sum(T)
  pi_hat <- as.numeric(T) / n

  # Multinomial variance
  cov_pi <- cov_multinom_phat(pi_hat, n)

  list(
    method = "Naive",
    K = K,
    n = n,
    pi_hat = pi_hat,
    cov = cov_pi,
    se = sqrt(pmax(diag(cov_pi), 0))
  )
}
