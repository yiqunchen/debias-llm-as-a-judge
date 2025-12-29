## debiasLLMReporting

R package for comparing methods that correct for measurement error when using LLM-as-a-judge (or any ML/AI/... classifiers) classifiers to estimate population prevalence.

### Overview

When using large language models (LLMs) as classifiers to estimate the prevalence of some outcome in a population, the LLM's imperfect accuracy introduces bias. This package provides simulation tools and estimators to study and correct for this measurement error.

We compare three main approaches:

1. **Measurement Error Correction**: Classical approach using the Rogan-Gladen estimator, which corrects for misclassification using estimated sensitivity (q1) and specificity (q0) from a calibration sample.

2. **Prediction-Powered Inference (PPI)**: Uses a small labeled dataset to debias predictions from a larger unlabeled dataset, combining human labels with LLM predictions.

3. **PPI++**: An extension of PPI that optimizes a tuning parameter (lambda) to minimize variance while maintaining valid confidence intervals.

4. **EIF**: An estimator based on the semi-parametric efficiency theory, which coincides with PPI++ with optimal tuning parameter in the binary case

### Installation

```r
# Install from GitHub
# install.packages("remotes")  # if needed
remotes::install_github("yiqunchen/debias-llm-as-a-judge")
```

### Tutorials and Use

Visit https://yiqunchen.github.io/debias-llm-as-a-judge/ for tutorials and examples. Please file an [issue](https://github.com/yiqunchen/debias-llm-as-a-judge/issues) if you have a request for a tutorial that is not currently included.

### Citation

If you use `debias-llm-as-a-judge` for your analysis, please cite our manuscript:
TBA
 
### Quick Start

```r
library(debiasLLMReporting)

# Simulated example
set.seed(42)
N <- 500; m <- 50
theta_true <- 0.6; q0 <- 0.85; q1 <- 0.80

# Generate data
Y_all <- rbinom(N, 1, theta_true)
Yhat_all <- ifelse(Y_all == 1, rbinom(N, 1, q1), 1 - rbinom(N, 1, q0))

# Split into calibration and test
idx_cal <- sample(N, m)
Y_cal <- Y_all[idx_cal]
Yhat_cal <- Yhat_all[idx_cal]
Yhat_test <- Yhat_all[-idx_cal]

# Apply PPI++
result <- ppi_pp_point_and_ci_general(
  Y_L = Y_cal, f_L = Yhat_cal, f_U = Yhat_test, alpha = 0.10
)
cat("Estimate:", round(result$theta, 3),
    "95% CI: [", round(result$ci_lower, 3), ",", round(result$ci_upper, 3), "]\n")
```


### Package Functions

Core estimators:

- `llm_point_and_ci()` - Rogan-Gladen estimator with delta method CI
- `ppi_point_and_ci()` - Standard PPI estimator
- `ppi_pp_point_and_ci_general()` - PPI++ with optimized lambda
- `eif_point_and_ci()` - Efficient influence function estimator
- `fit_misclass_mle()` - Joint MLE for misclassification model

### Reproducibility

The `reproduce/` directory contains scripts to replicate our analyses:

- `simulation_llm_vs_ppi.R` - Monte Carlo simulation comparing estimators across parameter grid
- `apply_estimators.R` - Apply estimators to real LLM judge evaluation data

Results are saved to timestamped directories under `results/`.

### License

MIT
