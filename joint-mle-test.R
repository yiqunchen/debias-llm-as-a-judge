## ============================================================
## Joint MLE diagnostics demo
## ============================================================

## When working inside the package repository without installing it,
##   source the package code so these helpers are available.
if (!exists("fit_misclass_mle")) {
  source("R/joint_mle.R")
}

## Align diagnostics with the smaller MLE simulation settings used in
## simulation_llm_vs_ppi_mle.R (lighter grid, lower q0/q1, smaller m_cal).
set.seed(123)
res <- sim_compare_info(
  B = 300,
  N = 3000,
  theta_true = 0.3,
  q0_true = 0.7,
  q1_true = 0.7,
  m_cal = round(0.01 * 3000),  # match label_ratio = 0.10
  level = 0.95
)

print(res$summary)
print(res$coverage)
print(res$theta_summary)

message("Access plots via res$plots, e.g., print(res$plots$theta_hist)")
