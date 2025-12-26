suppressPackageStartupMessages({
  library(dplyr)
})

source("R/joint_mle.R")
source("R/sim_estimators.R")

run_probe <- function(seed,
                      theta = 0.3,
                      q0 = 0.6,
                      q1 = 0.6,
                      sample_size = 1000L,
                      label_ratio = 0.1,
                      m0_prop = 0.5,
                      signal = 1) {
  set.seed(seed)
  n_labels <- max(4L, round(label_ratio * sample_size))
  m0 <- max(2L, round(n_labels * m0_prop))
  m1 <- max(2L, n_labels - m0)
  counts <- list(
    n_ppi = n_labels,
    N_ppi = sample_size,
    n_llm_test = sample_size,
    m0 = m0,
    m1 = m1
  )
  N_required <- max(counts$n_ppi + counts$N_ppi,
                    counts$n_llm_test + counts$m0 + counts$m1)
  dgp <- generate_dgp_data(N_required, theta, signal)
  Y <- dgp$Y
  Z <- ifelse(Y == 1,
              rbinom(N_required, 1, q1),
              rbinom(N_required, 1, 1 - q0))
  idx_neg <- which(Y == 0)
  idx_pos <- which(Y == 1)
  if (length(idx_neg) < counts$m0 || length(idx_pos) < counts$m1) {
    stop("Not enough calibration positives/negatives; rerun with larger sample.")
  }
  idx_cal0 <- sample(idx_neg, counts$m0)
  idx_cal1 <- sample(idx_pos, counts$m1)
  y_cal_idx <- c(idx_cal0, idx_cal1)
  y_cal <- Y[y_cal_idx]
  yhat_cal <- Z[y_cal_idx]
  yhat_test <- Z[seq_len(counts$n_llm_test)]
  fit <- fit_misclass_mle(y_cal, yhat_cal, yhat_test)
  tibble(
    seed = seed,
    theta_hat = fit$theta_hat,
    se_obs = fit$se_obs["theta"],
    se_exp = fit$se_exp["theta"],
    ci_obs_lo = fit$ci_theta_obs[1],
    ci_obs_hi = fit$ci_theta_obs[2]
  )
}

keep <- list()
seed <- 1
while (length(keep) < 5 && seed <= 50) {
  res <- run_probe(seed)
  if (!is.na(res$se_obs)) {
    keep[[length(keep) + 1]] <- res
  }
  seed <- seed + 1
}

results <- bind_rows(keep)
print(results)
